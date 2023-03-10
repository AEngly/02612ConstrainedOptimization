using LinearAlgebra;
using SparseArrays;

# The documentation can be found here: https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg-factorizations

function LUDecomposition(A)

    L,U,p = lu(A)
    return(L, U, p)

end

function QRDecomposition(A)

    Q,R = qr(A)
    return(Q, R)

end

function CholeskyFactorization(A)

    F = cholesky(A)
    L = F.L
    return(L)

end

# The following works for a sparse matrix

function LDLDecomposition(A)

    A = sparse(A)
    F = ldlt(A)
    p = F.p
    D = F.D
    L = F.L

    return(L, D, p)

end

# Cholesky

function CholeskyDecomposition(A)

    rows, columns = size(A)
    eigenvalueTolerance = 10^-7

    # Checks before proceeding
    if rows != columns
        display("Matrix is not rectangular.")
        return
    elseif A != transpose(A)
        display("Matrix is not symmetric.")
        return
    elseif !all(map(ev -> ev > eigenvalueTolerance, eigvals(A)))
        display("Matrix is not positive definite.")
        return
    end

    # Create zero matrix for L
    L = zeros(rows, columns)

    # Perform the Cholesky decomposition
    for i in 1:rows
        for k in 1:i
            tmp_sum = sum(L[i,j] * L[k,j] for j in 1:k)

            # See formulas: https://en.wikipedia.org/wiki/Cholesky_decomposition
            if (i == k) # Diagonal elements
                # LaTeX: l_{kk} = \sqrt{ a_{kk} - \sum^{k-1}_{j=1} l^2_{kj}}
                L[i,k] = sqrt(A[i,i] - tmp_sum)
            else
                # LaTeX: l_{ik} = \frac{1}{l_{kk}} \left( a_{ik} - \sum^{k-1}_{j=1} l_{ij} l_{kj} \right)
                L[i,k] = (1.0 / L[k,k] * (A[i,k] - tmp_sum))
            end

        end
    end
    return L

end

function LinearSolve(A, b, decomposition = "Cholesky")

    if decomposition == "Cholesky"

        L = CholeskyDecomposition(A)
        L_star = transpose(L)
        n,_ = size(A)

        # Step 1: Solve Ly = b with forward substitution
        y = zeros(n)

        y[1] = b[1] / L[1,1]

        for i in 2:n
            y[i] = (b[i] - sum([L[i,j]*y[j] for j in 1:(i-1)]))/L[i,i]
        end

        # Step 2: Solve L^{*}x = y by backward substitution
        x = zeros(n)

        x[n] = y[n] / L_star[n,n]

        for i in reverse([k for k in 1:(n-1)])
            x[i] = (y[i] - sum([L_star[i,j]*x[j] for j in (i+1):n]))/L_star[i,i]
        end

        return x

    end

end