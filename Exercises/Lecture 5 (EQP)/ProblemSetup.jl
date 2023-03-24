using LinearAlgebra;
using SparseArrays;
using CPUTime;
using Profile;
using Printf;
using Plots;

include("/Users/andreasengly/Documents/Danmarks Tekniske Universitet/Portfolio for GitHub/02612 Constrained Optimization/Factorizations/LUFactorization.jl");

# EXERCISE 1

function ProblemSetup(n,u,d0)

    H = diagm([1 for i in 1:(n+1)])
    A = zeros((n,n+1))

    # First row
    A[1,1] = -1
    A[1,n] = 1

    for j in 1:(n-2)
        i = j+1
        A[i,i-1] = 1
        A[i,i] = -1
    end

    # Last row
    A[n, n-1] = 1
    A[n, n] = -1
    A[n, n+1] = -1

    # Then we construct g
    g = [-u for i in 1:(n+1)]

    # Then we construct b
    b = [0 for i in 1:(n)]
    b[1] = -d0

    return H, g, A, b

end

function KKT_matrix(n,u,d0)

    H, g, A, b = ProblemSetup(n, u, d0)
    A = A'

    # Then we need to transpose A (as it is effectively transposed already)

    zero_rows = size(-A')[1]
    zero_columns = size(-A)[2]

    row1 = [H -A];
    row2 = [-A' zeros(zero_rows, zero_columns)]

    LHS = [row1;
           row2]

    #-A' zeros(size(A))

    RHS = -[g; b]

    return LHS, RHS, H, g, A, b

end

# Then we test it for n = 10.

LHS, RHS, H, g, A, b = KKT_matrix(3, 0.2, 1)

# 1.5) LU factorization

function LUSolver(LHS, RHS)
    L,U,p = LUDecomposition(LHS)
    u1 = U \ (L \ RHS[p]);
    return u1;
end

# 1.6) LDL factorization

function LDLSolver(LHS, RHS)
    L, D, p = LDLDecomposition(LHS)
    u2 = zeros((size(LHS)[2], 1))
    u2[p] = L' \ (D \ (L \ RHS[p]));
    return u2
end

# 1.7) Solving with Null-Space Method

function NullSpace(H, g, A, b)
    Q,R = QRDecomposition(A);
    m = size(A)[2]
    n = size(A)[1]
    Q1 = Q[:,1:m]
    Q2 = Q[:,m+1:n]
    R = R[1:m,1:m]
    xy = R' \ b
    xz = Q2'*H*Q2 \ (-Q2'*(H * Q1 * xy + g))
    x = Q1 * xy + Q2 * xz
    lambda = R \ Q1'*(H*x + g)
    return x, lambda
end

# 1.8) Solving with Range-Space Method (Schur-Complement Method)

# See page

function RangeSpace(H, g, A, b)
    L = CholeskyFactorization(H);
    HG = L \ (L' \ g);
    HA = L \ (L' \ A);
    lambda = inv(A' * HA) * (b + A'*HG);
    x = HA*lambda - HG;
    return x, lambda
end

# 1.9) Plotting

n = zeros(25)
LU_time = zeros(25)
LDL_time = zeros(25)
NullSpace_time = zeros(25)
RangeSpace_time = zeros(25)

counter = 1

for i in [k for k in 20:20:500]

    @printf "Iteration: %.f\n" i

    n[counter] = i;
    LHS, RHS, H, g, A, b = KKT_matrix(i, 0.2, 1);

    CPUtic();
    u1 = LUSolver(LHS, RHS);
    LU_time[counter] = CPUtoc();

    CPUtic();
    u2 = LDLSolver(LHS, RHS);
    LDL_time[counter] = CPUtoc();

    CPUtic();
    x, lambda = NullSpace(H, g, A, b);
    NullSpace_time[counter] = CPUtoc();

    CPUtic();
    x, lambda = RangeSpace(H, g, A, b);
    RangeSpace_time[counter] = CPUtoc();

    counter = counter + 1;

end

plot(n, [LU_time LDL_time NullSpace_time RangeSpace_time], title="Benchmarking Methods", label=["LU" "LDL" "Null-Space" "Range-Space"], linewidth=2)
savefig("BenchmarkingECQP.png")
plot(n, [LU_time LDL_time NullSpace_time RangeSpace_time], title="Benchmarking Methods", label=["LU" "LDL" "Null-Space" "Range-Space"], linewidth=2, ylims=(0,2))
savefig("BenchmarkingECQP_zoomed.png")

# EXERCISE 2

# 2.1) Make a contour plot

H = [2 0;
     0 2];

g = [-2; -5];

A = ([1 -2;
      -1 -2;
      -1 2])';

b = [-2; -6; -2];

function q(x1, x2)
    x = [x1; x2]
    return 1/2 * x' * H * x + x' * g
end

function c(x1, x2)
    x = [x1; x2]
    return 1/2 * x' * H * x + x' * g
end

x = -15:15
y = -15:15

contour(x, y, q, fill=true, plot_title="")
plot(X, f, fill = (the_max, 0.5, :auto))
#savefig("output/plot.png")

# Implementation of Active-Set Method for Convex QP.

function ActiveSetConvexQP(H, g, A, b, max_iter)

    Ak = A;
    bk = b;

    # First we have to compute the LHS and RHS
    LHS = [H -Ak;
        -Ak' zeros(size(Ak,2), size(Ak,2))];

    RHS = -[g; bk]

    # Then we can solve it using any of the four method we implemented earlier on.
    x0, lambda = RangeSpace(H, g, Ak, bk)

    xk = x0
    W = [i for i in 1:size(A, 2)]
    Wk = [i for i in 1:size(A, 2)]

    for k in 1:max_iter

        println("(iteration, value): ", k, xk)

        # First we have to compute the LHS and RHS
        LHS = [H -Ak;
            -Ak' zeros(size(Ak,2), size(Ak,2))];

        RHS = -[g; bk]

        # Then we can solve it using any of the four method we implemented earlier on.
        u = RangeSpace(H, g, Ak, bk)
        pk = u[1:size(H,1)]
        lambda = u[size(H,1)+1:end]

        if !all(pk .> 0) && !all(pk .< 0)

            if all(lambda .> 0)
                return xk
            else
                j = findmin(lambda)[2]
                Wk = setdiff(Wk, j)
                xk = xk
            end

        else

            constraints_to_check = setdiff(W, Wk)
            filtered_constraints = setdiff(W, Wk)

            for i in constraints_to_check

                if !(A[:,i]' * xk < 0)
                    filtered_constraints = setdiff(constraints_to_check, i)
                end

            end

            ak = min(1, [(b[i] - A[:,i]'*xk)/(A[:,i]'*pk) for i in filtered_constraints])
            xk = xk + ak * pk;

            if any(A' * xk .< b)
                constraintsToConsider = W[A' * xk .< b]
                randomIndex = rand(constraintsToConsider)
                Wk = union(Wk, W[randomIndex])
            else
                Wk = Wk
            end

        end

        # Then we have to change the constraints of the problem
        Ak = A[:,Wk]
        bk = b[Wk]

    end

    return xk

end


test = ActiveSetConvexQP(H, g, A, b, 100)