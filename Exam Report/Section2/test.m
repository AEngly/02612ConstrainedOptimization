function [xk] = test(H, g, A, b, C, dl, du, l, u, maxiter)

% settings
num_tol = 1e-9;

n = length(g);
mb = length(b);
mdl = length(dl);
mdu = length(du);
ml = length(l);
mu = length(u);

gi = [zeros(n,1); ones(2*mb,1); ones(mdl,1); ones(mdu,1)];
bi = [-b; dl; du;];
A1 = [A; eye(mb,mb); -eye(mb,mb); zeros(mdl,mb); zeros(mdu,mb) ];
if mdl > 0
    A2 = [C; zeros(2*mb,mdl); -eye(mdl,mdl); zeros(mdu,mdl) ];
else
    A2 = zeros(n+mb*2+mdl+mdu,0);
end
if mdu > 0
    A3 = [C; zeros(2*mb,mdu); zeros(mdl,mdu); eye(mdu,mdu) ];
else
    A3 = zeros(n+mb*2+mdl+mdu,0);
end

Ai = [A1 A2 A3];

if ml == 0
    lb = [-inf*ones(n,1); zeros(2*mb,1); zeros(mdl,1); zeros(mdu,1)];
else
    lb = [l; zeros(2*mb,1); zeros(mdl,1); zeros(mdu,1);];
end

if mu == 0
    ub = zeros(n+mb*2+mdl+mdu,0);
else
    ub = [u; inf*ones(2*mb,1); inf*ones(mdl,1); inf*ones(mdu,1);];
end

x0 = linprog(gi,[],[],Ai',bi,lb,ub);

[n,ma] = size(A);

A = [A C -C eye(length(l)) -eye(length(u))];
b = [-b; dl; -du; l; -u];

% preallocate
xk = x0(1:n);
At = A';
w = At*xk - b <= 1e-8;


% loop setup
terminate = false;
n = length(xk);


k = 0;
while ~terminate
    k= k+1;
    
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
    terminate = ~(k <maxiter);
end
end