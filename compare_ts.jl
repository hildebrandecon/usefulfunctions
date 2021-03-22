# Simulates an AR(1) process and compares it to a Tauchen discretisation

include("tauchen.jl")
include("simul_ar.jl")

#=  Structure
    z' = ρz + σ*ϵ', where ϵ ~ N(0,1) iid

    Inputs
    - ρ: autocorrelation of AR(1) process
    - σ: standard deviation of random component of AR(1) process 
    - S0: initial state
    - T: simulation periods

    Outputs
    - simulated AR(1) process

=#

function compare_ts(N::Integer, ρ::Real, σ::Real = 1, z0::Real = 0, T::Int64 = 100)

    z_ar = simul_ar(ρ, σ, z0, T);

    (Π, S) = tauchen(N, ρ, σ);
    dist = (S[N]-S[1])/(N-1);
    
    z_mc = similar(z_ar)
    
    ub = Array(S) .+ dist/2
    ub[N] = Inf

    for t = 1:length(z_mc)

        flag = true
        count = 0
        while flag 

            count += 1

            if z_ar[t] <= ub[count]
                z_mc[t] = S[count]
                flag = false
            end

        end

    end

    return (z_ar, z_mc, ub)

end

# Comparison
using Plots

N = 5; ρ = 0.5; σ = 2; z0 = 0; T = 100;

(z_ar, z_mc, ub) = compare_ts(N, ρ, σ, z0, T);

plot(0:T, [z_ar z_mc], legend=:false)
hline!(ub, linestyle = :dot, linecolor = :black)
