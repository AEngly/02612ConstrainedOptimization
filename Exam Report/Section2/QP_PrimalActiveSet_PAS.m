function [xk, k] = test_new(H, g, A, b, C, dl, du, l, u, maxiter)

% settings
num_tol = 1e-9;
bn = 1e9;

n = length(g);
mb = length(b);
mdl = length(dl);
mdu = length(du);
ml = length(l);
mu = length(u);

if ml == 0
    lb = [-bn*ones(n,1); zeros(2*mb,1); zeros(mdl,1); zeros(mdu,1)];
else
    lb = [l; zeros(2*mb,1); zeros(mdl,1); zeros(mdu,1);];
end

if mu == 0
    ub = zeros(n+mb*2+mdl+mdu,0);
    A5 = zeros(n+mb*2+mdl+mdu,0);
else
    ub = [u];
    A5 = [ -eye(mu); zeros(mb*2+mdl+mdu, mu)];
end

[n,ma] = size(A);

A1 = [ A; eye(mb); -eye(mb); zeros(mdl, mb); zeros(mdu, mb)];

if mdl > 0
    A2 = [ C; zeros(mb,mdl); zeros(mb,mdl); eye(mdl, mdl); zeros(mdu, mdl)];
else
    A2 = zeros(n+mb*2+mdl+mdu,0);
end

if mdu > 0
    A3 = [ -C; zeros(mb,mdu); zeros(mb,mdu); zeros(mdl, mdu); eye(mdu, mdu)];
else
    A3 = zeros(n+mb*2+mdl+mdu,0);
end

A4 = eye(n+mb*2+mdl+mdu);

Ai = [A1 A2 A3 A4 A5];
bi = [-b; dl; -du; lb; -ub];

Hi = 1e-3*eye(n+mb*2+mdl+mdu);

v = zeros(mb,1);
w = zeros(mb,1);

v(b < 0) = -b(b < 0);
w(b >= 0) = b(b >= 0);

t1 = max(zeros(mdl,1),dl);
t2 = max(zeros(mdu,1),-du);
x0 = [ zeros(n,1); v; w; t1; t2];

gi = [ zeros(n,1); ones(mb*2+mdl+mdu,1)];

[x0, k] = QP_primalActiveSet_core(Hi, gi, Ai, bi, x0, ma, maxiter, num_tol);

A = [A C -C eye(length(l)) -eye(length(u))];
b = [-b; dl; -du; l; -u];

x0 = x0(1:n);
[xk, k] = QP_primalActiveSet_core(H, g, A, b, x0, ma, maxiter, num_tol);


end