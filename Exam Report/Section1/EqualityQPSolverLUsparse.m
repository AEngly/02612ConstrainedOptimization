
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

    [LHS, RHS, size_x] = KKT_matrix(H, g, A, b);
    LHS = sparse(LHS);
    RHS = sparse(RHS);
    [L, U] = lu(LHS);
    y = L \ RHS; % Forward substitution
    x_lambda = U \ y;   % Backward substitution
    x = x_lambda(1:size_x);
    lambda = x_lambda(size_x+1:end);

end