function [H,g,A,b,x,lambda] = GeneratorECQP(n,alpha,beta,density)

% ---------------- DESCRIPTION --------------
%
% Name: Generator for ECQP
% Type: Generates random programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x input \
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

    m = round(beta*n);
    A = sprandn(n,m,density);
    A = full(A);
    M = sprandn(n,n,density);
    H = M*M' + alpha*eye(n,n);

    x = randn(n,1);
    lambda = randn(m,1);

    y = [H -A; -A' zeros(m,m)]*[x;lambda];

    g = -y(1:n);
    b = y(n+1:n+m);

end