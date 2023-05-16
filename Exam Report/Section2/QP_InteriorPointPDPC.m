function [x,k,x_k,y,z,s,info] = QP_InteriorPointPDPC(H,g,A,b,C,dl,du,l,u,maxiter,x0, varargin)
% ---------------- DESCRIPTION --------------
%
% Name: QP_ineq_box_InteriorPointPDPC   
% Type: Primal-Dual Predictor-Corrector Interior-Point QP Solver
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'*x + b = 0
%                   bl <= C' x <= bu
%                   l <=    x <= u
%
% Syntax: [x, y, info, z, s, iter] = QP_InteriorPointPDPC(H,g,A,b,C,d,x0,y0,z0,s0)
%
%         info = true   : Converged
%              = false  : Not Converged
%
% Created: 24.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------


% ---------------- Check input arguments --------------

n = length(g);

args = [];

if exist('verbose','var')
  verbose=true;
  args = [args, 'verbose'];
end

if exist('x_k','var')
    all = true;
    x_k = zeros(n, maxiter);
else
    x_k = [];
    all = false;
end

% ---------------- Setup tolerances --------------
tol_L = 1e-8;
tol_C = 1e-8;
tol_A = 1e-8;
tol_mu = 1e-8;

% ---------------- Find starting point --------------

%Check problem size
m = length(b);
% Set initial guesses for X0, y0, s0 and z0

x0 = zeros(n,1);
y0 = zeros(m,1);

% Permute matrices
Cbar = [full(C) full(-C) eye(length(l)) -eye(length(u))]; 
dbar = [dl; -du; l; -u];

s0 = ones(length(dbar),1);
z0 = ones(length(dbar),1);

% Calculate residuals
rL = H*x0+g-A*y0-Cbar*z0;
rA = -b-A'*x0;
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


% ---------------- iterate --------------


% Calculate residuals
rL = H*x+g-A*y-Cbar*z;
rA = -b-A'*x;
rC = s+dbar-Cbar'*x;
rsz = (s.*z);

% Setup constants
mc = length(dbar);

% Complementarity measure
mu = z'*s/mc;

% Setup loop
k=0;
terminate = (k > maxiter | norm(rL)<=tol_L & norm(rA)<=tol_A & norm(rC)<=tol_C & abs(mu)<=tol_mu );

% preallocation
%deltaxyaff = zeros(n+m,1);

while ~terminate
    k= k+1;
    %disp(k);
    % Compute LDL factorization of modified KKT system
    Czs = Cbar*diag(z./s);
    Hbar = H + Czs*Cbar';
    % Book mentions modified Cholesky here??
    KKT = [ Hbar -A; -A' zeros(m,m)];
    [L,D,p] = ldl(KKT,'lower','vector');

    %Affine direction
    rbarL = rL - Czs * (rC-s);
    rtemp = [-rbarL; -rA];
    deltaxyaff(p) = L'\( D \(L\rtemp(p)));
    deltaxaff = deltaxyaff(1:n)';
    %deltayaff = deltaxyaff(n+1:end)';
    deltazaff = -diag(z./s)*Cbar'*deltaxaff + diag(z./s)*(rC-s);
    deltasaff = -s -diag( s./z)*deltazaff;

    %compute alpha affine
    idxdeltazaff = find(deltazaff < 0.0);
    idxdeltasaff = find(deltasaff < 0.0);
    alphaaff = min([1.0 (-z(idxdeltazaff)./deltazaff(idxdeltazaff))' (-s(idxdeltasaff)./deltasaff(idxdeltasaff))']);

    % Duality gap and centering parameter
    muaff = (z+alphaaff*deltazaff)'*(s+ alphaaff*deltasaff)/mc;
    sigma = (muaff/mu)^3;

    % Affine-Centering-Correction Direction 
    rbarsz = rsz + deltasaff.*deltazaff-sigma*mu;
    rbarL = rL - Czs * (rC-diag(z)\rbarsz);
    rtemp =  [-rbarL; -rA];

    deltaxy(p) = L'\( D \(L\rtemp(p))); 
    deltax = deltaxy(1:n)';
    deltay = deltaxy(n+1:end)';
    deltaz = -diag(z./s)*Cbar'*deltax + diag(z./s)*(rC-diag(z)\rbarsz);
    deltas = -diag(z)\rbarsz-diag( s./z )*deltaz;

    %compute alpha 
    idxdeltaz = find(deltaz < 0.0);
    idxdeltas = find(deltas < 0.0);
    alpha = min([1.0 (-z(idxdeltaz)./deltaz(idxdeltaz))' (-s(idxdeltas)./deltas(idxdeltas))']);

    %Update iteration
    nabla = 0.995;
    alphabar = alpha*nabla;
    x = x + alphabar*deltax;
    y = y + alphabar*deltay;
    z = z + alphabar*deltaz;
    s = s + alphabar*deltas;
    
    if all
        x_k(:,k) = x;
    end
    %calculate residuals
    rL = H*x+g-A*y-Cbar*z;
    rA = -b-A'*x;
    rC = s+dbar-Cbar'*x;
    rsz = (s.*z);
    mu = z'*s/mc;
    
    %disp(x(1:2));
    % Check convergence
    terminate = (k >= maxiter | norm(rL)<=tol_L & norm(rA)<=tol_A & norm(rC)<=tol_C & abs(mu)<=tol_mu );
end

info = k <= maxiter;
end