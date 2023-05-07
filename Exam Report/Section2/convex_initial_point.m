n = length(g);
mb = length(b);
mdl = length(dl);
mdu = length(du);
ml = length(l);
mu = length(u);

Hi = [1e-8*eye(n) zeros(n,mdl+mdu+ml+mu); zeros(mdl+mdu+ml+mu, n) eye(mdl+mdu+ml+mu)];
gi = [zeros(n,1); ones(mdl,1); ones(mdu,1); ones(ml,1); ones(mu,1)];
bi = [b; -dl; du; -l; u];
A1 = [A; zeros(mb,mdl); zeros(mb,mdu); zeros(mb,ml); zeros(mb,mu) ];
A2 = [C; -eye(mdl,mdl); zeros(mdu,mdl); zeros(ml,mdl); zeros(mu,mdl) ];
A3 = [-C; zeros(mdl,mdu); -eye(mdu,mdu); zeros(ml,mdu); zeros(mu,mdu) ];
A4 = [eye(n,ml); zeros(mdl,ml); zeros(mdu,ml); -eye(ml,ml); zeros(mu,ml) ];
A5 = [eye(n,mu); zeros(mdl,mu); zeros(mdu,mu); -zeros(ml,mu); eye(mu,mu) ];
Ai = [A1 A2 A3 A4 A5];

[x, lambda] = EqualityQPSolverLUsparse(Hi, gi, Ai, bi);