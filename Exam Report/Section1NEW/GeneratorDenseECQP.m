function [H,g,A,b,x,lambda] = GeneratorDenseECQP(n,beta,density)

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
% Created: 17.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    m = round(beta*n);
    A = sprandn(n,m,density);
    A = full(A);
    H = sprandsym(n,density);
    Ii = eye(n);
    while ~all(eig(H) > 0)
        H = H + Ii;
    end
    x = randn(n,1);
    lambda = randn(m,1);

    y = [H -A; -A' zeros(m,m)]*[x;lambda];

    g = -y(1:n);
    b = y(n+1:n+m);

end