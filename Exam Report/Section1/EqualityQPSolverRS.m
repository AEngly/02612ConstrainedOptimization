function [x, lambda] = EqualityQPSolverRS(H, g, A, b)

% ---------------- DESCRIPTION --------------
%
% Name: Range Space Solver for EQP
% Type: Solver for Equality Constrained Quadratic Programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%
% Syntax: [x, lambda] = EqualityQPSolverRS(H, g, A, b)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    L = chol(H);
    HG = L \ (L' \ g);
    HA = L \ (L' \ A);
    lambda = (A' * HA) \ (-b + A'*HG);
    x = HA*lambda - HG;
    
end