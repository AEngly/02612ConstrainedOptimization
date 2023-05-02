function [x,z,s] = QP_general_InteriorPointPDPC_initial_point(H,g,A,b,C,dl,du,l,u)

%Check problem size
n = length(g);
m = length(b);
% Set initial guesses for X0, y0, s0 and z0

if isempty(l)
    x0 = zeros(n,1);
else
    x0 = mean(l,u);
end

if isempty(b)
    y0 = [];
else
    y0 = zeros(m,1);
end

% Permute matrices
Cbar = [full(C) full(-C) eye(n,n) -eye(n,n)]; 
dbar = [-dl; du; -l; u];

s0 = 

% Calculate residuals
rL = H*x0+g-A*y0-Cbar*z0;
rA = b-A'*x0;
rC = s0+dbar-Cbar'*x0;
rsz = (s0.*z0);

%Compute LDL factorization of modified KKT system
Czs = Cbar*diag(z0./s0);
Hbar = H + Czs*Cbar';
KKT = [ Hbar -A; -A' zeros(m,m)];
[L,D,p] = ldl(KKT,'lower','vector');

%Affine direction
rbarL = rL - Czs * (rC-s0);
rtemp = [-rbarL; -rA];
deltaxyaff(p) = L'\( D \(L\rtemp(p)));
deltaxaff = deltaxyaff(1:n)';
deltazaff = -diag(z0./s0)*Cbar'*deltaxaff + diag(z0./s0)*(rC-s0);
deltasaff = -s0 -diag( s0./z0)*deltazaff;

% Compute initial point
x = x0;
y = y0;
z = max(ones(length(z0),1),abs(z0+deltazaff));
s = max(ones(length(s0),1),abs(s0+deltasaff));
end