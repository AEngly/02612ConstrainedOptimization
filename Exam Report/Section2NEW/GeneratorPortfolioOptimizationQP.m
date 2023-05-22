function [H,g,A,b,C,dl,du,l,u] = GeneratorPortfolioOptimizationQP(k,gamma,density)

% ---------------- DESCRIPTION --------------
%
% Name: Generator for Optimal Control QP
% Type: Generates random programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%                   dl <= C'x <= du
%                   l <= x <= u
%
% Syntax: [H,g,A,b,C,dl,du,l,u] = GeneratorOptimalControlQP(n,alpha,beta,density)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------
    
    % Define risk-aversion
    m = k;
    n = 100*k;
    % Generate problem
    D = diag(rand(1,n))*sqrt(k);
    F = sprandn(n,k,density);
    mu = randn(n,1);


    
    % Create empty A and b, as no equality constraints are considered
    b = [zeros(k,1); 1];
  
    A1 = [full(F); -eye(k);];
    A2 = [ones(n,1); zeros(k,1)];
    A = [A1 A2];

    C = zeros(n,0);
    dl = zeros(m,0);
    du = zeros(m,0);
    % Create empty upper and lower bounds for x, as x is unconstrained
    l = [zeros(n,1); -inf*ones(k,1)];
    u = zeros(n,0);

    H = [D zeros(n,k); zeros(k,n) eye(k)];

    g = -[gamma^(-1)*mu; zeros(k,1)];

end