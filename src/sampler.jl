"""
    ODE_MH(u_t::Function, t_m::AbstractVector{<:AbstractFloat}, y::AbstractMatrix{<:AbstractFloat}, t_span::Tuple{Float64,Float64}, K::Int, K_b::Int, k_d::Int, f_theta!::Function, g_theta!::Function, log_pdf_w_theta::Function, log_pdf_theta::Function, theta_0::AbstractVector{<:AbstractFloat}, log_pdf_x_init::Function, x_init_0::AbstractVector{<:AbstractFloat}, propose_z::Function, log_proposal_ratio_z::Function; ODE_solver=RK4(), ODE_solver_opts=(dt=0.1, adaptive=false), print_progress=true)

Run a Metropolis-Hastings (MH) sampler with ODE-integrated latent trajectory to obtain samples ``\\{\\theta, x(t=0)\\}^{[1:K]}`` from the joint parameter and state posterior distribution ``p(\\theta, x(t=0) \\mid \\mathbb{D}=\\{u(t), y_{1:M}\\})``.

# Arguments
- `u_t`: input trajectory; function of time ``t``
- `t_m`: time points of the output measurements, must lie within `t_span`
- `y`: output measurements
- `t_span`: timespan of the training trajectory
- `K`: number of models/scenarios to be sampled
- `K_b`: length of the burn in period (in MH iterations)
- `k_d`: number of MH iterations to be skipped between retained samples to decrease correlation (thinning)
- `f_theta!`: dynamics function parametrized by theta (mutating); has inputs ``(\\dot{x}, \\theta, x, u, t)``
- `g_theta!`: measurement function parametrized by theta (mutating); has inputs ``(g, \\theta, x, u, t)``
- `log_pdf_w_theta`: function that returns the logarithm of the probability density function of the measurement noise parametrized by theta; has inputs ``(\\theta, w)``
- `log_pdf_theta`: function that returns the logarithm of the probability density function of theta (prior); has input ``(\\theta)``
- `theta_0`: ``\\theta`` used to initialize the MH sampler
- `log_pdf_x_init`: function that returns the logarithm of probability density function of ``x(t=-T)`` (prior); has input ``x(t=-T)``
- `x_init_0`: ``x(t=-T)`` used to initialize the MH sampler
- `propose_z`: function that proposes new parameters ``z' = [\\theta'; x'(t=-T)]`` (proposal distribution); has input (``z = [\\theta; x(t=-T)``])
- `log_proposal_ratio_z`: function that computes the logarithm of the ratio of proposal densities for ``z``; has input arguments ``(z, z')``
- `ODE_solver`: ODE solver algorithm to use (we use RK4 as default as this is also used in the optimal control formulation)
- `ODE_solver_opts`: ODE solver setting
- `print_progress`: if true, the progress is printed (default: `true`)

# Returns
- `MH_samples`: MH samples
- `acceptance_ratio`: acceptance ratio of the MH sampler
- `runtime`: runtime of the sampling process
"""
function ODE_MH(u_t::Function, t_m::AbstractVector{<:AbstractFloat}, y::AbstractMatrix{<:AbstractFloat}, t_span::Tuple{Float64,Float64}, K::Int, K_b::Int, k_d::Int, f_theta!::Function, g_theta!::Function, log_pdf_w_theta::Function, log_pdf_theta::Function, theta_0::AbstractVector{<:AbstractFloat}, log_pdf_x_init::Function, x_init_0::AbstractVector{<:AbstractFloat}, propose_z::Function, log_proposal_ratio_z::Function; ODE_solver=RK4(), ODE_solver_opts=(dt=0.1, adaptive=false), print_progress=true)
    # Total number of MH iterations (length of the Markov chain including the initial state)
    N = K_b + 1 + (K - 1) * (k_d + 1)

    # Get number of parameters, etc.
    n_theta = length(theta_0)
    n_x = length(x_init_0)
    M = size(y, 2)

    # Initialize and pre-allocate.
    MH_samples = Vector{MH_sample}(undef, K)
    for k in 1:K
        MH_samples[k] = MH_sample(Array{Float64}(undef, n_theta), Array{Float64}(undef, n_x), Array{Float64}(undef, n_x))
    end
    n_accepted = 0
    current_sample = 1
    t_save = vcat(t_m, t_span[2])

    # Simulate system forward using a numerical ODE solver and compute the log-likelihood of the measurements.
    function simulate_log_likelihood(theta, x_init)
        # Update model, i.e., update right hand side of ODE, measurement function, and noise distribution.
        ode_rhs(dx, x, p, t) = f_theta!(dx, theta, x, u_t(t), t)
        g!(g, x, u, t) = g_theta!(g, theta, x, u, t)
        log_pdf_w(w) = log_pdf_w_theta(theta, w)

        prob = ODEProblem(ode_rhs, x_init, t_span, theta)
        sol = solve(prob, ODE_solver; ODE_solver_opts..., saveat=t_save)
        x_sim = Array(sol)

        # Likelihood based on measurement model (logarithms are used for numerical reasons).
        log_likelihood_sim = 0.0
        g_sim = Array{Float64}(undef, size(y, 1))
        w = Array{Float64}(undef, size(y, 1))
        for m in 1:M
            g!(g_sim, x_sim[:, m], u_t(t_m[m]), t_m[m])
            w .= y[:, m] - g_sim
            log_likelihood_sim += log_pdf_w(w)
        end

        return log_likelihood_sim, x_sim[:, end]
    end

    # Store the k-th state of the chain if the burn-in period is reached and the sample is not removed by thinning.
    function store_sample!(k, theta, x_init, x_t)
        if (K_b < k) && (mod(k - (K_b + 1), k_d + 1) == 0)
            MH_samples[current_sample].theta .= theta
            MH_samples[current_sample].x_t .= x_t
            MH_samples[current_sample].x_init .= x_init
            current_sample += 1
        end
    end

    # Time MH sampler.
    sampling_timer = time()

    if print_progress
        println("### Started MH sampling")
    end

    # Initialize the chain with z[1] = [theta_0; x_init_0].
    theta = copy(theta_0)
    x_init = copy(x_init_0)
    z = vcat(theta, x_init)
    log_p_theta = log_pdf_theta(theta)
    log_p_x_init = log_pdf_x_init(x_init)
    if !isfinite(log_p_theta) || !isfinite(log_p_x_init)
        error("Initial value z[1] = [theta_0; x_init_0] has zero prior density.")
    end
    log_likelihood, x_t = simulate_log_likelihood(theta, x_init)
    store_sample!(1, theta, x_init, x_t)

    for k in 1:N-1
        # Propose new parameters and a new initial state.
        z_prop = propose_z(z)
        theta_prop = z_prop[1:n_theta]
        x_init_prop = z_prop[n_theta+1:end]
        log_p_theta_prop = log_pdf_theta(theta_prop)
        log_p_x_init_prop = log_pdf_x_init(x_init_prop)

        # Proposals outside the support of the prior are rejected without simulating the system.
        accepted = false
        if isfinite(log_p_theta_prop) && isfinite(log_p_x_init_prop)
            log_likelihood_prop, x_t_prop = simulate_log_likelihood(theta_prop, x_init_prop)

            # Compute acceptance probability.
            log_acceptance_ratio = log_likelihood_prop - log_likelihood + log_p_theta_prop - log_p_theta + log_p_x_init_prop - log_p_x_init + log_proposal_ratio_z(z, z_prop)

            # Accept or reject the proposal.
            if log(rand()) < log_acceptance_ratio
                accepted = true
                n_accepted += 1
                theta = theta_prop
                x_init = x_init_prop
                z = z_prop
                x_t = x_t_prop
                log_likelihood = log_likelihood_prop
                log_p_theta = log_p_theta_prop
                log_p_x_init = log_p_x_init_prop
            end
        end

        # z[k+1] is the proposal if accepted, otherwise the current state is repeated.
        store_sample!(k + 1, theta, x_init, x_t)

        # Print progress.
        if print_progress
            @printf("%s%i/%i iterations (%i accepted)\e[0m\n", accepted ? "\e[32m" : "\e[31m", k + 1, N, n_accepted)
        end
    end

    runtime = time() - sampling_timer
    acceptance_ratio = (n_accepted / (N - 1)) * 100

    # Print runtime and acceptance ratio.
    if print_progress
        @printf("### MH sampling complete\nRuntime: %.2f s\nAcceptance ratio: %.2f %%\n", runtime, acceptance_ratio)
    end

    return MH_samples, acceptance_ratio, runtime
end

"""
    staged_ODE_MH(u_t::Function, t_m::AbstractVector{<:AbstractFloat}, y::AbstractMatrix{<:AbstractFloat}, t_span::Tuple{Float64,Float64}, K::Int, K_b::Int, k_d::Int, f_theta!::Function, g_theta!::Function, log_pdf_w_theta::Function, log_pdf_theta::Function, theta_0::AbstractVector{<:AbstractFloat}, log_pdf_x_init::Function, x_init_0::AbstractVector{<:AbstractFloat}, proposal_z_cov_0::AbstractMatrix{<:AbstractFloat}, M_chunk::Int, K_stage::Int, alpha::Union{AbstractFloat,AbstractVector{<:AbstractFloat}}; regularizer::Float64=1e-8, ODE_solver=RK4(), ODE_solver_opts=(dt=0.1, adaptive=false), print_progress=true)

Run a Metropolis-Hastings (MH) sampler with ODE-integrated latent trajectory with incremental data and adaptive proposal to obtain samples ``\\{\\theta, x(t=0)\\}^{[1:K]}`` from the joint parameter and state posterior distribution ``p(\\theta, x(t=0) \\mid \\mathbb{D}=\\{u(t), y_{1:M}\\})``.
The number of data points used in the likelihood computation is gradually increased by a fixed chunk size. At each stage, the MH sampler is run on the current data subset, and the proposal distribution is adapted based on the empirical covariance of the collected samples.

# Arguments
- `u_t`: input trajectory; function of time ``t``
- `t_m`: time points of the output measurements, must lie within `t_span`
- `y`: output measurements
- `t_span`: timespan of the training trajectory
- `K`: number of models/scenarios to be sampled
- `K_b`: length of the burn in period
- `k_d`: number of models/scenarios to be skipped to decrease correlation (thinning)
- `f_theta!`: dynamics function parametrized by theta (mutating); has inputs ``(\\dot{x}, \\theta, x, u, t)``
- `g_theta!`: measurement function parametrized by theta (mutating); has inputs ``(g, \\theta, x, u, t)``
- `log_pdf_w_theta`: function that returns the logarithm of the probability density function of the measurement noise parametrized by theta; has inputs ``(\\theta, w)``
- `log_pdf_theta`: function that returns the logarithm of the probability density function of theta (prior); has input ``(\\theta)``
- `theta_0`: ``\\theta`` used to initialize the MH sampler
- `log_pdf_x_init`: function that returns the logarithm of probability density function of ``x(t=-T)`` (prior); has input ``x(t=-T)``
- `x_init_0`: ``x(t=-T)`` used to initialize the MH sampler
- `proposal_z_cov_0`: initial covariance of the multivariate normal proposal distribution for ``z = [\\theta; x(t=-T)]`` (e.g., covariance of the prior)
- `M_chunk`: number of data points added at each stage
- `K_stage`: number of samples per stage
- `alpha`: proposal scaling factor (scalar or vector; if a vector its i‐th element is used at stage i)
- `regularizer`: small constant added to the diagonal of the proposal covariance matrix to ensure positive definiteness (default: `1e-8`)
- `ODE_solver`: ODE solver algorithm to use (we use RK4 as default as this is also used in the optimal control formulation)
- `ODE_solver_opts`: ODE solver setting
- `print_progress`: if true, the progress is printed (default: `true`)

# Returns
- `MH_samples`: final samples from full-data posterior
- `acceptance_ratio`: vector containing the acceptance ratio of each stage
- `runtime`: total runtime of the staged sampling process
"""
function staged_ODE_MH(u_t::Function, t_m::AbstractVector{<:AbstractFloat}, y::AbstractMatrix{<:AbstractFloat}, t_span::Tuple{Float64,Float64}, K::Int, K_b::Int, k_d::Int, f_theta!::Function, g_theta!::Function, log_pdf_w_theta::Function, log_pdf_theta::Function, theta_0::AbstractVector{<:AbstractFloat}, log_pdf_x_init::Function, x_init_0::AbstractVector{<:AbstractFloat}, proposal_z_cov_0::AbstractMatrix{<:AbstractFloat}, M_chunk::Int, K_stage::Int, alpha::Union{AbstractFloat,AbstractVector{<:AbstractFloat}}; regularizer::Float64=1e-8, ODE_solver=RK4(), ODE_solver_opts=(dt=0.1, adaptive=false), print_progress=true)
    # Get number of parameters, etc.
    n_theta = length(theta_0)
    n_x = length(x_init_0)
    n_z = n_theta + n_x # total number of variables
    M = size(y, 2)

    # Determine the number of stages (each stage adds M_chunk data points)
    N_stages = ceil(Int, M / M_chunk)

    # Initialize current parameter vector
    theta = theta_0
    x_init = x_init_0
    proposal_cov = proposal_z_cov_0

    # Allocate an array to store samples from each stage.
    acceptance_ratio = zeros(N_stages)
    MH_samples = Vector{MH_sample}(undef, K)

    sampling_timer = time()

    if print_progress
        println("### Started staged MH sampling")
    end

    for i in 1:N_stages
        # Select current data chunk - use data from 1 to T_i.
        M_i = min(i * M_chunk, M)
        y_i = y[:, 1:M_i]
        t_m_i = t_m[1:M_i]

        # Define the proposal function as sampling from a multivariate normal.
        propose_z(z) = rand(MvNormal(z, proposal_cov))
        log_ratio_proposal_pdf(z_accepted, z_prop) = 0

        # Call the base MH sampler.
        if i < N_stages
            # For intermediate stages, sample K_stage samples without thinning.
            MH_samples_stage, acceptance_ratio_stage = ODE_MH(u_t, t_m_i, y_i, t_span, K_stage, K_b, 0, f_theta!, g_theta!, log_pdf_w_theta, log_pdf_theta, theta, log_pdf_x_init, x_init, propose_z, log_ratio_proposal_pdf; ODE_solver=ODE_solver, ODE_solver_opts=ODE_solver_opts, print_progress=false)
        else
            # In the final stage, sample K samples with thinning parameter k_d.
            MH_samples_stage, acceptance_ratio_stage = ODE_MH(u_t, t_m_i, y_i, t_span, K, K_b, k_d, f_theta!, g_theta!, log_pdf_w_theta, log_pdf_theta, theta, log_pdf_x_init, x_init, propose_z, log_ratio_proposal_pdf; ODE_solver=ODE_solver, ODE_solver_opts=ODE_solver_opts, print_progress=false)
        end

        # Save stage samples and acceptance ratio.
        acceptance_ratio[i] = acceptance_ratio_stage

        if i < N_stages
            # Update proposal covariance based on the empirical covariance
            Z = hcat([vcat(s.theta, s.x_init) for s in MH_samples_stage]...)

            post_cov_theta = cov(transpose(Z))
            if isa(alpha, Number)
                proposal_cov = alpha * post_cov_theta + regularizer * Matrix(I, n_z, n_z)
            else
                proposal_cov = alpha[i] * post_cov_theta + regularizer * Matrix(I, n_z, n_z)
            end

            # Update the current state to the last sample from the current stage.
            theta = MH_samples_stage[end].theta
            x_init = MH_samples_stage[end].x_init
        else
            MH_samples = MH_samples_stage
        end
        if print_progress
            @printf("Stage %i/%i complete\nAcceptance ratio: %.2f %%\n", i, N_stages, acceptance_ratio_stage)
        end
    end

    # Print runtime and average acceptance ratio.
    average_acceptance_ratio = mean(acceptance_ratio)
    runtime = time() - sampling_timer
    if print_progress
        @printf("### Staged MH sampling complete\nRuntime: %.2f s\nAverage acceptance ratio: %.2f %%\n",
            runtime, average_acceptance_ratio)
    end

    return MH_samples, acceptance_ratio, runtime
end
