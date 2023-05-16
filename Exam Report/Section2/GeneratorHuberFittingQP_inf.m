function [H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP_inf(n,beta,density)

% ---------------- DESCRIPTION --------------
%
% Name: Generator for Huber Fitting QP
% Type: Generates random programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%                   dl <= C'x <= du
%                   l <= x <= u
%
% Syntax: [H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP(n,alpha,beta,density)
%
%
% Created: 05.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------
    
    bn = 1e15;
    m = round(beta*n);

    Ax = zeros(n,m);
    while n ~= rank(Ax)
    Ax = full(sprandn(n,m,density));
    end

    rank(Ax)
    v = randn(n,1)/n;
    e = randn(m,1)/4;
    e(randsample(m,floor(m/20))) = rand(1,floor(m/20))*10;
    
    b = -(Ax'*v+e);

    Au = -eye(m);
    As = eye(m);
    Ar = -eye(m);

    A = [Ax; Au; As; Ar];

    C = zeros(n+3*m,0);
    dl = zeros(0,1);
    du = zeros(0,1);

    % lower limit for x and u set to a small enough number.
    l = [-inf*ones(n,1); -inf*ones(m,1); zeros(m,1); zeros(m,1);];
    u = inf*ones(n+3*m,0);
    
    e = 1e-8;
    H = zeros(n+3*m,n+3*m);
    %H(1:n,1:n) = e*eye(n);
    H(n+1:n+m,n+1:n+m) = eye(m);
    %H(n+m+1:n+2*m,n+m+1:n+2*m) = e*eye(m);
    %H(n+2*m+1:n+3*m,n+2*m+1:n+3*m) = e*eye(m);
    

    g = [zeros(n,1); zeros(m,1); ones(m,1)*2; ones(m,1)*2];

end