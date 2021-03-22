# Tauchen (1986)'s method for approximating AR(1) process with finite markov chain
# Inspirations: Miao (2020) and QuantEcon

using Distributions

#=  Structure
    z' = ρz + σ*ϵ', where ϵ ~ N(0,1) iid

    Inputs
    - N: number of states in discrete state space
    - ρ: autocorrelation of AR(1) process
    - σ: standard deviation of random component of AR(1) process 
    - n_std: span of discrete state space relative to standard deviation

    Outputs
    - transition matrix Π
    - discretised state space z

=#

function tauchen(N::Integer, ρ::Real, σ::Real = 1, n_std::Integer = 4)

    # Step 1: discretise the state space
    σ_z = sqrt( σ / (1 - ρ^2) )
    zmax = n_std * σ_z
    zmin = -zmax

    dist = (zmax - zmin)/(N - 1)
    z = zmin:dist:zmax

    # Step 2: choose N intervals 
    # this is done on the fly in the third step

    # Step 3: calculate transition probabilities
    # note: we evaluate the probabilities using a standard normal distribution
    Π = zeros(Float64, N, N)

    Distr = Normal(0,1)

    for row = 1:N

        Π[row, 1] = cdf(Distr, σ^(-1) * (z[1] + dist/2 - ρ*z[row]) )
        Π[row, N] = 1 - cdf(Distr, σ^(-1) * (z[N] - dist/2 - ρ*z[row]) )
        for col = 2:N-1
            Π[row, col] = (cdf(Distr, (z[col] + dist/2 - ρ*z[row]) / σ ) 
                            - cdf(Distr, (z[col] - dist/2 - ρ*z[row]) / σ))
    
        end

    end

    if sum(abs.(sum(Π, dims = 2) .- 1)) > 2e-15
        error("Rows do not sum to one.")
    end

    return (Π, Array(z))

end 

# Example

using Plots
include("simul_mc.jl");

ρ = 0.5; σ = 2; N = 7; T = 100;

(Π, S) = tauchen(N, ρ, σ);
z_mc = simul_mc(Π, S, 4, T);

plot(0:T, z_mc)
