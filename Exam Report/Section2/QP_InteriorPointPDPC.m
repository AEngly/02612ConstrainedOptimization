function [x,y,z,s,info,iter] = QP_general_InteriorPointPDPC(H,g,A,b,C,dl,du,l,u)
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

itermax = 20;

% Setup tolerances
tol_L = 1e-8;
tol_C = 1e-8;
tol_mu = 1e-8;

%Find starting point
[x,y,z,s] = QP_general_InteriorPointPDPC_initial_point(H,g,A,b,C,dl,du,l,u);

%Check problem size
n = length(x);
m = length(b);

% Permute matrices
Cbar = [full(C) full(-C) eye(length(l),length(l)) -eye(length(l),length(l))]; 
dbar = [-dl; du; -l; u];

% Calculate residuals
rL = H*x+g-A*y-Cbar*z;
rA = -b-A'*x;
rC = s-dbar-Cbar'*x;
rsz = (s.*z);

% Setup constants
mc = length(dbar);

% Complementarity measure
mu = z'*s/mc;

% Setup loop
k=0;
terminate = (k > itermax | norm(rL)<=tol_L | norm(rC)<=tol_C | abs(mu)<=tol_mu );

% preallocation
%deltaxyaff = zeros(n+m,1);

while ~terminate
    k= k+1;
    disp(k);
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
    
    %calculate residuals
    rL = H*x+g-A*y-Cbar*z;
    rA = -b-A'*x;
    rC = s-dbar-Cbar'*x;
    rsz = (s.*z);
    mu = z'*s/mc;

    % Check convergence
    terminate = (k > itermax | norm(rL) <= tol_L | norm(rC) <= tol_C | abs(mu) <= tol_mu );
end
iter = k;
info = k <= itermax;
end