function [x, lambda] = EqualityQPSolverLDLdense(H, g, A, b)

% ---------------- DESCRIPTION --------------
%
% Name: EqualityQPSolverLDLdense 
% Type: Solver for Equality Constrained Quadratic Programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%
% Syntax: [x, lambda] = EqualityQPSolverLDLdense(H, g, A, b)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------

    [LHS, RHS,n,m] = KKT_matrix(H, g, A, b);
    x_lambda = zeros(n+m,1);
    [L,D,p] = ldl(LHS, 'lower','vector');
    x_lambda(p) = L'\(D\(L\RHS(p)));
    x = x_lambda(1:n);
    lambda = x_lambda(n+1:end);

end