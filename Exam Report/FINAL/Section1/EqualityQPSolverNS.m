function [x, lambda] = EqualityQPSolverNS(H, g, A, b)

% ---------------- DESCRIPTION --------------
%
% Name: Null Space Solver for EQP
% Type: Solver for Equality Constrained Quadratic Programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%
% Syntax: [x, lambda] = EqualityQPSolverNS(H, g, A, b)
%
% Assumptions:
%             1) A has full row rank
%             2) Z'HZ is positive definite (Z is null-space of A)
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    [Q, R] = qr(A);
    [n, m] = size(A);
    Y = Q(:,1:m);
    Z = Q(:,m+1:n);
    R = R(1:m,1:m);
    L = chol(Z'*H*Z);
    xy = R'\-b;
    xz = L\(L'\(-Z'*(H*Y*xy + g)));
    x = Y*xy + Z*xz;
    lambda = R\Y'*(H*x + g);

end