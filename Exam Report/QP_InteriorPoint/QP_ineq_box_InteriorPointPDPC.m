function [x,y,info,z,s, iter] = QP_ineq_box_InteriorPointPDPC(H,g,A,b,C,d,x0,y0,z0,s0)
% ---------------- DESCRIPTION --------------
%
% Name: QP_dualActiveSet   
% Type: Primal-Dual Active-Set QP Solver
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    bl <= A' x <= bu
%                   l <=    x <= u
%
% Syntax: [x,y,info,z,s, iter] = QP_InteriorPointPDPC(H,g,A,b,C,d,x0,y0,z0,s0)
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

% Setup tolerances
tol_L = 1e-8;
tol_C = 1e-8;
tol_mu = 1e-8;



% Calculate residuals
rL = H*x0+g-A*y0-C*z0;
rA = b-A'*x0;
rC = s0+d-C'*x0;
rsz = (s0.*z0)';

% Setup constants
mc = length(d);
mu = z0'*s0/mc;

x = x0;
s = y0;
z = z0;
s = s0;
deltaxyaff = zeros(length(g)+length(b),1);

% Setup loop
terminate = (k >20 | norm(rL)<=tol_L | norm(rC)<=tol_C | abs(mu)<=tol_mu );
k=0;

while ~terminate
    Czs = C*diag(z./s);
    Hbar = H + Czs*C';
    O = zeros(length(b),length(b));
    [L,D,p] = ldl([Hbar -A; -A' O],'lower','vector');

    % Affine direction
    rbarL = rL - Czs * (rC-s);
    rtemp = [ -rbarL;-rA];
    deltaxyaff(p) = L'\( D \(L\rtemp(p)));
    deltaxaff = deltaxyaff(1:n);
    deltayaff = deltaxyaff(n+1:end);
    deltazaff = -(z./s)*C'*deltaxaff + (z./s)*(rC-s');
    deltasaff = -s' -( s./z)*deltazaff;

    % Compute alpha affine
    idxdeltazaff = find(deltazaff < 0.0);
    idxdeltasaff = find(deltasaff < 0.0);
    alphaaff = min([1.0 -z(idxdeltazaff)/deltazaff(idxdeltazaff) -s(idxdeltasaff)/deltasaff(idxdeltasaff)]);

    % Duality gap and centering parameter
    muaff = (z+alphaaff*deltazaff)'*(s+ alphaaff*deltasaff)/mc;
    sigma = (muaff/mu)^3;

    % Affine-Centering-Correction Direction 
    rbarsz = rsz + deltasaff.*deltazaff-sigma*mu;
    rbarL = rL - Czs * (rC-inv(diag(z))*rbarsz);
    rtemp = [ -rbarL;-rA];

    deltaxy(p) = L'\( D \(L\rtemp(p)));
    deltax = deltaxy(1:n);
    deltay = deltaxy(n+1:end);
    deltaz = -(z./s)*C'*deltax + (z./s)*(rC-diag(z)\rbarsz);
    deltas = -diag(z)\rbarsz-( s./z )*deltaz;

    % Compute alpha 
    idxdeltaz = find(deltaz < 0.0);
    idxdeltas = find(deltas < 0.0);
    alpha = min([1.0 -z(idxdeltaz)/deltaz(idxdeltaz) -s(idxdeltas)/deltas(idxdeltas)]);

    % Update iteration
    nabla = 0.995;
    alphabar = alpha*nabla;
    x = x + alphabar*deltax;
    y = y + alphabar*deltay;
    z = z + alphabar*deltaz;
    s = s + alphabar*deltas;
    
    % Calculate residuals
    rL = H*x+g-A*y-C*z;
    rA = b-A'*x;
    rC = s+d-C'*x;
    rsz = (s.*z)';
    mu = z'*s/mc;

    % Check convergence
    terminate = (k >100 | norm(rL)<=tol_L | norm(rC)<=tol_C | abs(mu)<=tol_mu );
end
end