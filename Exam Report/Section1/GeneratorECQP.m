function [H,g,A,b] = GeneratorECQP(n,alpha,beta,density)

% ---------------- DESCRIPTION --------------
%
% Name: Generator for ECQP
% Type: Generates random programs
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

    m = round(beta*n);
    A = sprandn(n,m,density);
    b = rand(m,1);
    M = sprandn(n,n,density);
    H = M*M' + alpha*eye(n,n);
    % g should be calculated from random x and lambda
    g = randn(n,1);

end