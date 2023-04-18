function [x,z,s] = QP_ineq_box_InteriorPointPDPC_initial_point(H,g,dl,C,du,l,u,x0,z0,s0)

%Check problem size
n = length(x0);

% Permute matrices
Cbar = [full(C) full(-C) eye(n,n) -eye(n,n)]; 
dbar = [-dl; du; -l; u];

% Calculate residuals
rL = H*x0+g-Cbar*z0;
rC = s0-dbar-Cbar'*x0;
rsz = (s0.*z0);

%Compute LDL factorization of modified KKT system
Czs = Cbar*diag(z0./s0);
Hbar = H + Czs*Cbar';
[L,D,p] = ldl(Hbar,'lower','vector');

%Affine direction
rbarL = rL - Czs * (rC-s0);
rtemp = [ -rbarL];
deltaxaff(p) = L'\( D \(L\rtemp(p)));
deltazaff = -diag(z0./s0)*Cbar'*deltaxaff' + diag(z0./s0)*(rC-s0);
deltasaff = -s0 -diag( s0./z0)*deltazaff;

% Compute initial point
x = x0;
z = max(ones(length(z0),1),abs(z0+deltazaff));
s = max(ones(length(s0),1),abs(s0+deltasaff));
end