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
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    [Q, R] = qr(A);
    [n, m] = size(A);
    Q1 = Q(:,1:m);
    Q2 = Q(:,m+1:n);
    R = R(1:m,1:m);
    xy = R' \ -b;
    xz = Q2'*H*Q2 \ (-Q2'*(H * Q1 * xy + g));
    x = Q1 * xy + Q2 * xz;
    lambda = R \ Q1'*(H*x + g);
    
end