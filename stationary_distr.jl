using LinearAlgebra

Π = [0.9 0.1 0.0; 0.4 0.4 0.2; 0.1 0.1 0.8];


#=  Stationary distribution satisfies:
        ones(N,1) = (I - Π' + ones(N,N)) π
                b = A π
                π = inv(A) * b
=#

function stationary_distr(Π)

    N = size(Π,1)

    A = (I - Π' + ones(N,N))
    b = ones(N,1)
    
    π = inv(A)*b
    
    if sum(abs.(sum(Π'*π - π, dims = 2))) > 2e-15
        error("Something went wrong here.")
    end
    
    return π

end







