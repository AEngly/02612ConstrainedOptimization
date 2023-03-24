% ---------------- DESCRIPTION --------------
%
% Name: driver.m  
% Type: Testing Script for Primal-Dual Active-Set QP Solver
%
% Problem structure:
%          min  1/2*x'Hx + g'x
%           
%          s.t.   A'x = b      (Lagrange multiplier: lambda1)
%                 dl <= C'x <= du      (Lagrange multiplier: lambda2 & lambda3)
%                 l <= x <= u      (Lagrange multiplier: lambda4 & lambda5)
%
%
% Created: 24.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          IMM, Technical University of Denmark

% We start by generating random matrices 

alpha = 0.10;
density = 0.15;
n = 100;

[H,g,bl,A,bu,l,u] = RandomQP(n,alpha,density);

[all_xk, mu_star, active_constraints] = QP_ineqBox_primalActiveSet(H, g, [], [], A, bl, bu, l, u);