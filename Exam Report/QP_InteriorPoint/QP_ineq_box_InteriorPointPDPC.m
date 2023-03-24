function [x, info,z,s, iter] = QP_ineq_box_InteriorPointPDPC(H,g,dl,C,du,l,u,x0,z0,s0)
% ---------------- DESCRIPTION --------------
%
% Name: QP_ineq_box_InteriorPointPDPC   
% Type: Primal-Dual Predictor-Corrector Interior-Point QP Solver
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    bl <= A' x <= bu
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

% TODO
% Test algorithmn.
% Implement starting point heuristic
% Setup for general problem formulation
% Improve performance
% Add info
% Return iter

itermax = 20;

% Setup tolerances
tol_L = 1e-8;
tol_C = 1e-8;
tol_mu = 1e-8;

n = length(x0);

Cbar = [full(C) full(-C) eye(n,n) -eye(n,n)]; 
dbar = [-dl; du; -l; u];

% Calculate residuals
rL = H*x0+g-Cbar*z0;
rC = s0-dbar-Cbar'*x0;
rsz = (s0.*z0);

% Setup constants
mc = length(dbar);
mu = z0'*s0/mc;

x = x0;
z = z0;
s = s0;

% Setup loop
k=0;
terminate = (k > itermax | norm(rL)<=tol_L | norm(rC)<=tol_C | abs(mu)<=tol_mu );


while ~terminate
    k= k+1;
    Czs = Cbar*diag(z./s);
    Hbar = H + Czs*Cbar';
    [L,D,p] = ldl(Hbar,'lower','vector');

    %Affine direction
    rbarL = rL - Czs * (rC-s);
    rtemp = [ -rbarL];
    deltaxaff(p) = L'\( D \(L\rtemp(p)));
    deltazaff = -diag(z./s)*Cbar'*deltaxaff' + diag(z./s)*(rC-s);
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
    rtemp = [ -rbarL];

    deltax = L'\( D \(L\rtemp(p))); %Might have error
    %deltax = deltax';
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
    z = z + alphabar*deltaz;
    s = s + alphabar*deltas;
    
    %calculate residuals
    rL = H*x+g-Cbar*z;
    rC = s-dbar-Cbar'*x;
    rsz = (s.*z);
    mu = z'*s/mc;

    % Check convergence
    terminate = (k > itermax | norm(rL) <= tol_L | norm(rC) <= tol_C | abs(mu) <= tol_mu );
end
iter = k;
info = k <= itermax;
end