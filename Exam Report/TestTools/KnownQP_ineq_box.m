function [xstar,H,g,bl,A,bu,l,u] = KnowQP_ineq_box(n,m,alpha,density)
% KnowQP_ineq_box  Generates data for a known inequality and box constrained Convex QP
%
%   min     0.5 x' H x + g' x
%    x
%   s.t.    bl <= A' x <= bu
%            l <=    x <= u
%
% Syntax: [xstar,H,g,bl,A,bu,l,u] = KnowQP_ineq_box(n,m,alpha,density)
%
%   Inputs:
%       n       number of variabels
%       m       number of constraints
%
%   Outputs: QP data (xstar,H,g,bl,A,bu,l,u)

%% Must be implemented
A = sprandn(n,m,density);

x = zeros(n,1);
x(1:m,1) = abs(rand(m,1));

bl = -rand(m,1);
bu = rand(m,1);
M = sprandn(n,n,density);
H = M*M' + alpha*eye(n,n);
g = randn(n,1);
l = -ones(n,1);
u = ones(n,1);