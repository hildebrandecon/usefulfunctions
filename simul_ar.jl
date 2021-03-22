# Simulates an AR(1) process

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

function simul_ar(ρ::Real, σ::Real = 1, z0::Real = 0, T::Int64 = 100)

    z_ar = zeros(Real, T + 1);
    z_ar[1] = z0;

    for t in 1:T
        z_ar[t+1] = ρ * z_ar[t] + σ * randn()
    end

    return z_ar

end

# Example
using Plots

ρ = 0.5; σ = 2; z0 = 0; T = 100;
z_ar = simul_ar(ρ, σ, z0, T);

plot(0:T, z_ar)

