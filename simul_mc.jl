# Simulates a discrete Markov chain
# Inspirations: Miao (2020) and QuantEcon

using Distributions

#=  Structure

    Inputs
    - Π: transition matrix
    - S: state space
    - S0: initial state
    - T: simulation periods

    Outputs
    - simulated Markov chain

=#

function simul_mc(Π, S, S0::Int64 = 1, T::Int64 = 100)

    N = length(S)
    dists = [Categorical(Π[i, :]) for i in 1:N]

    z_mc = zeros(Int64, T + 1);
    z_mc[1] = S0;

    for t in 1:T
        dist = dists[z_mc[t]] 
        z_mc[t+1] = rand(dist)
    end

    z_mc = S[z_mc]

end

# Example
using Plots

Π = [0.9 0.1 0.0; 0.4 0.4 0.2; 0.1 0.1 0.8];
S = [-1, 0, 1];

z_mc = simul_mc(Π, S, 2, 100);

plot(z_mc)
