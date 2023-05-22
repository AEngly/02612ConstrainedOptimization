function [x, k, x_k] = QP_PrimalActiveSet_PAS(H, g, A, b, C, dl, du, l, u, maxiter,varargin)
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

% ---------------- Settings --------------

num_tol = 1e-9;
bn = 1e9;

% ---------------- Setup for initial point --------------
[n,ma] = size(A);
n = length(g);
mb = length(b);
mdl = length(dl);
mdu = length(du);
ml = length(l);
mu = length(u);

gbar = [zeros(n,1); ones(2*mb+mdl+mdu+ml+mu,1)];
bbar = [-b];

A1 = [A; eye(mb,mb); -eye(mb,mb); zeros(mdl+mdu+ml+mu,mb)];

if mdl > 0
    C1 = [C; zeros(2*mb,mdl); eye(mdl,mdl); zeros(mdu+ml+mu,mdl)];
    t1 = max(zeros(mdl,1),dl);
else
    C1 = zeros(n+mb*2+mdl+mdu+ml+mu,0);
    t1 = zeros(0,1);
end
if mdu > 0
    C2 = [-C; zeros(2*mb,mdu); zeros(mdl,mdu); -eye(mdu,mdu); zeros(ml+mu,mdu)];
    t2 = max(zeros(mdu,1),-du);
else
    C2 = zeros(n+mb*2+mdl+mdu+ml+mu,0);
    t2 = zeros(0,1);
end


if ml == 0
    ly = [-bn*ones(n,1); zeros(2*mb+mdl+mdu+ml+mu,1)];
else
    ly = [l; zeros(2*mb+mdl+mdu+ml+mu,1)];
end

C3 = [eye(n) zeros(n,2*mb+mdl+mdu+ml+mu); [zeros(2*mb+mdl+mdu,n); eye(ml,n); zeros(mu,n)] eye(2*mb+mdl+mdu+ml+mu)];
s1 = max(zeros(n,1),l);

if mu == 0
    C4 = zeros(n,0);
    s2 = zeros(0,1);
else
    C4 = [-eye(mu); zeros(mdl+mdu+ml,mu); eye(mu)];
    s2 = max(zeros(mu,1),-u);
end

dbar = [ dl; -du; ly; -u];

Abar = [A1 C1 C2 C3 C4];
bbar = [-b; dl; -du; ly; -u];

Hi = 1e-3*eye(n+mb*2+mdl+mdu+ml+mu);

v = zeros(mb,1);
w = zeros(mb,1);

v(b >= 0) = -b(b >= 0);
w(b < 0) = b(b < 0);

x0 = [ zeros(n,1); v; w; t1; t2; s1; s2];

gi = [ zeros(n,1); ones(2*mb+mdl+mdu+ml+mu,1)];

% ---------------- Find initial point --------------
[x0, k] = QP_primalActiveSet_core(Hi, gi, Abar, bbar, x0, ma, maxiter, num_tol);

% ---------------- Setup for solution --------------
Abar = [A C -C eye(length(l)) -eye(length(u))];
bbar = [-b; dl; -du; l; -u];
x0 = x0(1:n);
% ---------------- Find solution --------------
[x, k, x_k]= QP_primalActiveSet_core(H, g, Abar, bbar, x0, ma, maxiter, num_tol,args);


end