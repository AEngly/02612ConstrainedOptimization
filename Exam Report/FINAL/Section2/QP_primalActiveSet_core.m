function [xk, k, x_k] = QP_primalActiveSet_core(H, g, A, b, x0, ma, maxiter, num_tol, varargin)

% ---------------- DESCRIPTION --------------
%
% Name: QP_primalActiveSet_core  
% Type: Primal-Dual Active-Set QP Solver
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'*x + b = 0
%                   dl <= C' x <= du    (Lagrange multiplier: mu)
%                   l <=    x <= u      (Lagrange multiplier: lamba)  
%
% Syntax: [xk, k, x_k] = QP_primalActiveSet_core(H, g, A, b, x0, ma, maxiter, num_tol, varargin)
%
%         info = true   : Converged
%              = false  : Not Converged
%
% Created: 24.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          IMM, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

n = length(g);

args = [];

if ~exist('verbose','var')
  verbose=true;
  args = [args, 'verbose'];
end

if ~exist('x_k','var')
    full = true;
    x_k = zeros(n, maxiter);
else
    x_k = [];
    full = false;
end

%  ---------------- Pre-allocate --------------

xk = x0;
At = A';
w = abs(At*xk - b) < num_tol;

%  ---------------- Loop setup --------------
% 
terminate = false;
k = 1;

%  ---------------- Main loop --------------
while ~terminate

    
    mk = sum(w);
    %disp(k);
    

    [L, U, pe] = lu([H -A(:,w); -At(w,:) zeros(mk, mk)],'vector');
    RHS = -[H*xk + g;zeros(mk,1)];
    z = L \ RHS(pe); % Forward substitution
    ptemp = U \ z;   % Backward substitution

    pstar = ptemp(1:n);
    muk = ptemp(n+ma+1:n+mk);

    if norm(pstar) < num_tol
        if all(muk >= 0)
            k=k-1;
            return;
        else
            active = find(w == 1);
            w(active(ma + find(muk == min(muk(:))))) = 0;
        end
    else
        nom = (b(~w) - At(~w,:)*xk);
        denom = (At(~w,:)*pstar);
        idx =denom < 0;
        alphas = nom(idx)./denom(idx);
        j = alphas == min(alphas);
        alpha = alphas(j);

        if alpha < 1
            xk = xk + alpha.*pstar;
            m = find(idx == 1);
            l = find(w == 0);
            w(l(m(j))) = 1;

        else
            xk = xk + pstar;
        end
    end
    terminate = ~(k <maxiter); % terminate is max iterations is met.


    if full
        x_k(:,k) = xk;
    end
    k= k+1;
end

end