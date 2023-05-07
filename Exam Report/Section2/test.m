function [xk] = test(H, g, A, b, C, dl, du, l, u, maxiter)

% settings
num_tol = 1e-8;

x0 = linprog(g,[C'; -C'],[du; -dl],A',b,l,u);

A = [A C -C eye(length(l)) -eye(length(u))];
b = [-b; dl; -du; l; -u];

% preallocate
At = A';
w = At*x0 - b <= 1e-8;
xk = x0;

% loop setup
terminate = false;
n = length(xk);
k = 0;
while ~terminate
    k= k+1;
    
    mk = sum(w);
    disp(k);
    

    [L, U, pe] = lu([H -A(:,w); -At(w,:) zeros(mk, mk)],'vector');
    RHS = -[H*xk + g;zeros(mk,1)];
    z = L \ RHS(pe); % Forward substitution
    ptemp = U \ z;   % Backward substitution

    pstar = ptemp(1:n);
    muk = ptemp(n+1:n+mk);

    if norm(pstar) < num_tol
        if all(muk >= 0)
            return;
        else
            active = find(w == 1);
            w(active(find(muk == min(muk(:))))) = 0;
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
    terminate = ~(k <maxiter);
end
end