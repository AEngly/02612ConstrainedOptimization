function [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density)

% ---------------- DESCRIPTION --------------
%
% Name: Generator for Random QP
% Type: Generates random programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%                   dl <= C'x <= du
%                   l <= x <= u
%
% Syntax: [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    m = round(beta*n);
    b = [];
    A = [];
    C = sprandn(n,m,density);
    dl = -rand(m,1);
    du = rand(m,1);
    l = [];
    u = [];
    M = sprandn(n,n,density);
    H = M*M' + alpha*eye(n,n);
    g = randn(n,1);

end