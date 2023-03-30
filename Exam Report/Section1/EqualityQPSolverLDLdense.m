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

    [LHS, RHS, size_x] = KKT_matrix(H, g, A, b);
    [L,D,~] = ldl(LHS);
    z = L \ RHS; % Forward substitution
    y = D \ z; % Forward substitution
    x_lambda = L' \ y;   % Backward substitution
    x = x_lambda(1:size_x);
    lambda = x_lambda(size_x+1:end);

end