function [x, lambda] = EqualityQPSolverLUsparse(H, g, A, b)

% ---------------- DESCRIPTION --------------
%
% Name: EqualityQPSolverLUsparse 
% Type: Solver for Equality Constrained Quadratic Programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%
% Syntax: x = EqualityQPSolverLUsparse(LHS, RHS)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------

    [LHS, RHS, n, m] = KKT_matrix(H, g, A, b);
    x_lambda = zeros(n+m,1);
    LHS = sparse(LHS); % Changing the KKT-matrix to sparse representation
    [L, U, p] = lu(LHS,'vector');
    x_lambda = U\(L\RHS(p)); % Solve with LU
    x = x_lambda(1:n);
    lambda = x_lambda(n+1:end);
    
end