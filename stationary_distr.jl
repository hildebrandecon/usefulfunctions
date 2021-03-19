# Simulates a discrete Markov chain
# Inspirations: Miao (2020), Stachurski (2009) and QuantEcon

using LinearAlgebra

#=  Stationary distribution satisfies:
        (I - Π' + ones(N,N)) * π = ones(N,1)
                           A * π = b
                               π = inv(A) * b

    Inputs
    - Π: transition matrix

    Outputs
    - invariant distribution

=#

function invariant(Π, check::Bool = false)

    N = size(Π,1)

    A = (I - Π' + ones(N,N))
    b = ones(N,1)
    
    πstar = inv(A)*b
    
    if sum(abs.(sum(Π'*πstar - πstar, dims = 2))) > 2e-15
        error("Something went wrong here.")
    end

    if abs(sum(πstar) - 1) > 2e-15
        error("Something went wrong here.")
    end

    if check

        π0 = repeat([1/N], N)
        dist = Inf
        count = 0
        
        while dist > 2e-15

            count += 1
            π1 = Π' * π0
            dist = norm(π1 - π0)
            π0 = π1

        end

        if norm(π1 - πstar) > 2e-10
            error("Something went wrong here.")
        end

    end

    return πstar

end

# Example
Π = [0.9 0.1 0.0; 0.4 0.4 0.2; 0.1 0.1 0.8];

πstar = invariant(Π);

norm(πstar - Π' * πstar)
