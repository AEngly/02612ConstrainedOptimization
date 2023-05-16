function [x, k, x_k] = QP_PrimalActiveSet_linprog(H, g, A, b, C, dl, du, l, u, maxiter, varargin)

% ---------------- DESCRIPTION --------------
%
% Name: QP_PrimalActiveSet_linprog   
% Type: Primal Active-Set QP Solver
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'*x + b = 0
%                   dl <= C' x <= du    (Lagrange multiplier: mu)
%                   l <=    x <= u      (Lagrange multiplier: lamba)  
%
% Syntax: [x,info,mu,lambda,iter] = QP_dualActiveSet(g,A,b,x)
%
%         info = true   : Converged
%              = false  : Not Converged
%
% Created: 24.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          IMM, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------
%
% ---------------- Check input arguments --------------



args = [];

if ~exist('verbose','var')
  args = [args, "verbose"];
end

if ~exist('x_k','var')
  args = [args, "x_k"];
end

% % ---------------- Settings --------------
num_tol = 1e-9;
bn = 1e9;

% ---------------- Initial point --------------
n = length(g);
mb = length(b);
mdl = length(dl);
mdu = length(du);
ml = length(l);
mu = length(u);

gi = [zeros(n,1); ones(2*mb,1); ones(mdl,1); ones(mdu,1)];
bi = [-b];
di = [ -dl; du;];
Ai = [A; eye(mb,mb); -eye(mb,mb); zeros(mdl,mb); zeros(mdu,mb) ];

if mdl > 0
    C1 = [-C; zeros(2*mb,mdl); -eye(mdl,mdl); zeros(mdu,mdl) ];
else
    C1 = zeros(n+mb*2+mdl+mdu,0);
end
if mdu > 0
    C2 = [C; zeros(2*mb,mdu); zeros(mdl,mdu); -eye(mdu,mdu) ];
else
    C2 = zeros(n+mb*2+mdl+mdu,0);
end

Ci = [C1 C2];

if ml == 0
    lb = [-bn*ones(n,1); zeros(2*mb,1); zeros(mdl,1); zeros(mdu,1)];
else
    lb = [l; zeros(2*mb,1); zeros(mdl,1); zeros(mdu,1);];
end

if mu == 0
    ub = zeros(n+mb*2+mdl+mdu,0);
else
    ub = [u; bn*ones(2*mb,1); bn*ones(mdl,1); bn*ones(mdu,1);];
end

x0 = linprog(gi',Ci',di,Ai',bi,lb,ub);

% ---------------- optimal point --------------
[n,ma] = size(A);

A = [A C -C eye(length(l)) -eye(length(u))];
b = [-b; dl; -du; l; -u];

x0 = x0(1:n);

[x, k, x_k] = QP_primalActiveSet_core(H, g, A, b, x0, ma, maxiter, num_tol,args);

end