function [x,z,s] = PDPC_InteriorPoint__ineq(H,g,C,d,x0,z0,s0)

% ---------------- DESCRIPTION --------------
%
% Name: QP_dualActiveSet   
% Type: Primal-Dual Active-Set QP Solver
%
% Problem structure:
%          min  g'*x
%           x
%          s.t. A x  = b      (Lagrange multiplier: mu)
%                 x >= 0      (Lagrange multiplier: lamba)
%
% Syntax: [x,info,mu,lambda,iter] = QP_dualActiveSet(g,A,b,x)
%
%         info = true   : Converged
%              = false  : Not Converged
%
% Created: 24.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

% Setup tolerances
tol_L = 1e-8;
tol_C = 1e-8;
tol_mu = 1e-8;

% implement starting point heuristic

%Calculate residuals
rL = H*x0+g-C*z0;
rC = s0+d-C'*x0;
rsz = (s0.*z0);

% Setup constants
mc = length(d);
mu = z0'*s0/mc;

x = x0;
s = s0;
z = z0;

% Setup loop
terminate = (k >20 | norm(rL)<=tol_L | norm(rC)<=tol_C | abs(mu)<=tol_mu );
k=0;

while ~terminate
    k= k+1;
    disp(k);
    Czs = C*diag(z./s);
    Hbar = H + Czs*C';
    [L,D,p] = ldl([Hbar],'lower','vector');

    %Affine direction
    rbarL = rL - Czs * (rC-s);
    rtemp = [ -rbarL];
    deltaxaff(p) = L'\( D \(L\rtemp(p)));
    deltazaff = -diag(z./s)*C'*deltaxaff' + diag(z./s)*(rC-s);
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

    deltax = L'\( D \(L\rtemp(p)));
    %deltax = deltax';
    deltaz = -diag(z./s)*C'*deltax + diag(z./s)*(rC-diag(z)\rbarsz);
    deltas = -diag(z)\rbarsz-diag( s./z )*deltaz;

    %compute alpha 
    idxdeltaz = find(deltaz < 0.0);
    idxdeltas = find(deltas < 0.0);
    alpha = min([1.0 (-z(idxdeltaz)./deltaz(idxdeltaz))' (-s(idxdeltas)./deltas(idxdeltas))']);

    %Update iteration
    nabla = 0.95;
    alphabar = alpha*nabla;
    x = x + alphabar*deltax;
    z = z + alphabar*deltaz;
    s = s + alphabar*deltas;
    
    disp(x);
    %calculate residuals
    rL = H*x+g-C*z;
    rC = s+d-C'*x;
    rsz = (s.*z);
    mu = z'*s/mc;
    terminate = (k >20 | norm(rL)<=tol_L | norm(rC)<=tol_C | abs(mu)<=tol_mu );

end
end