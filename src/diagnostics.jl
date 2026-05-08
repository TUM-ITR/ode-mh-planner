"""
    compute_autocorrelation(MH_samples::Vector{MH_sample}; max_lag::Int=100)

Compute the autocorrelation function (ACF) of the MH samples.

# Arguments
- `MH_samples`: MH samples
- `max_lag`: maximum lag at which to calculate the ACF

# Returns
- `autocorrelation`: matrix containing the ACF for each variable
"""
function compute_autocorrelation(MH_samples::Vector{MH_sample}; max_lag::Int=100)
    # Get number of models.
    K = size(MH_samples, 1)

    # Get number of parameters of the MH samples.
    n_variables = length(MH_samples[1].theta) + length(MH_samples[1].x_init)

    # Fill matrix with the series of the parameters of the MH samples.
    sample_matrix = Array{Float64}(undef, K, n_variables)
    for i in 1:K
        sample_matrix[i, :] .= vcat(MH_samples[i].theta, MH_samples[i].x_init)
    end

    # Calculate the autocorrelation.
    autocorrelation = autocor(sample_matrix, Array(0:max_lag); demean=true)

    return autocorrelation
end

"""
    compute_ess(MH_samples::Vector{MH_sample}; max_lag::Int=100)

Compute the effective sample size (ESS) for each parameter and initial state.

# Arguments
- `MH_samples`: MH samples
- `max_lag`: maximum lag for autocorrelation estimation

# Returns
- `ess`: vector of ESS estimates for all variables
"""
function compute_ess(MH_samples::Vector{MH_sample}; max_lag::Int=100)
    # Get number of models.
    K = size(MH_samples, 1)

    # Get number of parameters of the MH samples.
    n_variables = length(MH_samples[1].theta) + length(MH_samples[1].x_init)

    # Compute autocorrelation.
    autocorrelation = compute_autocorrelation(MH_samples; max_lag=max_lag)

    ess = zeros(n_variables)

    for i in 1:n_variables
        # Sum autocorrelation of i-th variable until first negative or max_lag.
        autocorrelation_sum = 0.0
        for lag in 1:max_lag
            if autocorrelation[lag+1, i] < 0
                break
            end
            autocorrelation_sum += autocorrelation[lag+1, i]
        end

        # Compute the effective sample size.
        ess[i] = K / (1 + 2 * autocorrelation_sum)
    end

    return ess
end

"""
    compute_gelman_rubin(MH_chains::Vector{Vector{MH_sample}})

Compute the Gelman–Rubin statistic R̂ for each parameter and initial state from a vector of PMCMC chains.
The Gelman–Rubin statistic R̂ quantifies convergence by comparing within-chain to between-chain variance.
R̂ close to 1 (typically R̂ < 1.05) indicates good convergence across chains.

# Arguments
- `MH_chains`: vector of chains, where each chain is a vector of MH samples

# Returns
- `R_hat`: vector of R̂ values, one for each variable
"""
function compute_gelman_rubin(MH_chains::Vector{Vector{MH_sample}})
    N = length(MH_chains) # Number of chains
    K = length(MH_chains[1]) # Number of samples per chain

    # Get number of parameters of the MH samples.
    n_variables = length(MH_chains[1][1].theta) + length(MH_chains[1][1].x_init)

    # Extract samples from each chain
    sample_matrices_chains = Array{Float64}[]
    for chain in MH_chains
        sample_matrix = Array{Float64}(undef, K, n_variables)
        for i in 1:K
            sample_matrix[i, :] .= vcat(chain[i].theta, chain[i].x_init)
        end

        push!(sample_matrices_chains, sample_matrix)
    end

    # Compute means and variances
    means = zeros(N, n_variables)
    variances = zeros(N, n_variables)
    for n in 1:N
        means[n, :] .= vec(mean(sample_matrices_chains[n], dims=1))
        variances[n, :] .= vec(var(sample_matrices_chains[n], dims=1, corrected=true))
    end

    # Between-chain and within-chain variance
    mean_overall = mean(means, dims=1)
    B = K / (N - 1) .* sum((means .- mean_overall) .^ 2, dims=1)
    W = mean(variances, dims=1)

    # Estimated marginal posterior variance and R̂
    V_hat = (K - 1) / K .* W .+ B / K
    R_hat = sqrt.(V_hat ./ W)

    # Clip from below at 1.0 for numerical consistency
    return max.(R_hat, 1.0)
end