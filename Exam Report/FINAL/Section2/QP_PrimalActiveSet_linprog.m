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
% Note: Initial point is found with linprog.
%
% Syntax: [x, k, x_k] = QP_PrimalActiveSet_linprog(H, g, A, b, C, dl, du, l, u, maxiter, varargin)
%
%         info = true   : Converged
%              = false  : Not Converged
%
% Created: 24.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          IMM, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

args = [];

if ~exist('verbose','var')
  args = [args, "verbose"];
end

if ~exist('x_k','var')
  args = [args, "x_k"];
end

% % ---------------- Settings --------------
num_tol = 1e-8;
bn = 1e9;

% ---------------- Initial point --------------
n = length(g);
mb = length(b);
mdl = length(dl);
mdu = length(du);
ml = length(l);
mu = length(u);

gbar = [zeros(n,1); ones(2*mb+mdl+mdu+ml+mu,1)];
bbar = [-b];

Abar = [A; eye(mb,mb); -eye(mb,mb); zeros(mdl+mdu+ml+mu,mb)];

if mdl > 0
    C1 = [C; zeros(2*mb,mdl); eye(mdl,mdl); zeros(mdu+ml+mu,mdl)];
else
    C1 = zeros(n+mb*2+mdl+mdu+ml+mu,0);
end
if mdu > 0
    C2 = [-C; zeros(2*mb,mdu); zeros(mdl,mdu); -eye(mdu,mdu); zeros(ml+mu,mdl)];
else
    C2 = zeros(n+mb*2+mdl+mdu+ml+mu,0);
end


if ml == 0
    ly = [-bn*ones(n,1); zeros(2*mb+mdl+mdu+ml+mu,1)];
else
    ly = [l; zeros(2*mb+mdl+mdu+ml+mu,1)];
end

C3 = [eye(n) zeros(n,2*mb+mdl+mdu+ml+mu); [zeros(2*mb+mdl+mdu,n); eye(ml,n); zeros(mu,n)] eye(2*mb+mdl+mdu+ml+mu)];

if mu == 0
    C4 = zeros(n,0);
else
    C4 = [-eye(mu); zeros(mdl+mdu+ml,mu); eye(mu)];
end

dbar = [ dl; -du; ly; -u];

Cbar = [C1 C2 C3 C4];
x0 = linprog(gbar',-Cbar',-dbar,Abar',bbar);
%fprintf("Starting point:\n")
%disp(x0);
%fprintf("Starting point end.\n")

% ---------------- optimal point --------------
[n,ma] = size(A);

Abar = [A C -C eye(length(l)) -eye(length(u))];
bbar = [-b; dl; -du; l; -u];

x0 = x0(1:n);

[x, k, x_k] = QP_primalActiveSet_core(H, g, Abar, bbar, x0, ma, maxiter, num_tol,args);

end