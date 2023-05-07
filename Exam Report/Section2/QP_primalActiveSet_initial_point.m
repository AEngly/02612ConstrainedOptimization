
function [all_xk] = QP_primalActiveSet_initial_point(H, g, A, b, C, dl, du, l, u)

% ---------------- DESCRIPTION --------------
%
% Name: QP_dualActiveSet   
% Type: Primal-Dual Active-Set QP Solver
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'*x + b = 0
%                   dl <= C' x <= du    (Lagrange multiplier: mu)
%                   l <=    x <= u      (Lagrange multiplier: lamba)  
%
% Syntax: [x,info,mu,lambda,iter] = QP_dualActiveSet(g,A,b,x)
%
%         info = true   : Converged
%              = false  : Not Converged
%
% Created: 24.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          IMM, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

% ---------------- Initial point --------------

% Add slack variable to convert to equality constrained problem
n = length(g);

x = linprog(g,[C'; -C'],[du; -dl],A',b,l,u);

% We start by defining the KKT-matrix

A_bar = [A C -C eye(length(l)) -eye(length(u))];
b_bar = [-b; dl; -du; l; -u];

optimal = false;
counter = 1;
alpha = 0;
num_tol = 1.0e-6;
iter = 0;
verbose = false;

% Then we define the initial working set
xk = x;
A_bar_t = A_bar';
Wk = A_bar_t*xk - b_bar <= 1e-8;
[n, ma] = size(A);
[n, m] = size(A_bar);
all_constraints = (1:m)';

% All xk
all_xk = xk;

while ~optimal && iter < 100

    iter = iter + 1;

    % We start with a naive implementation - and later extend to
    % null-space method.

    Ak = A_bar_t(Wk,:)';
    bk = b_bar(Wk);
    gk = H * xk + g;   
    mk = sum(Wk);

    row1 = [H -Ak];
    row2 = [-Ak' zeros(mk, mk)];
    LHS = [row1; row2];
    RHS = -[gk; zeros(mk,1)];

    [L, U, pe] = lu(LHS,'vector');
    z = L \ RHS(pe); % Forward substitution
    p_mu = U \ z;   % Backward substitution
    disp(iter);
    p = p_mu(1:n);
    muk = p_mu(n+1:n+mk);

    % Then we use the null-space method (see slides EqualityConstrainedQP
    % from week 5).

    % If step size is 0, then we need to find the multipliers
    if norm(p) < num_tol
        if all(muk > 0)
            x_star = xk;
            mu_star = zeros(m,1);
            mu_star(Wk,1) = muk;
            all_constraints = (1:m)';
            active_constraints = all_constraints(Wk);
            optimal = true;
        else
            % Then identify index of most negative multipliers
            remove_index = muk == min(muk);
            active_indices = find(Wk == 1);
            % Update working set
            Wk(active_indices(remove_index),:) = 0;
            % Keep the current solution
            xk = xk;
        end
    else

        % Compute the distance to the nearest inactive constraint in the search direction
        nominator = b_bar(~Wk) - A_bar_t(~Wk,:)*xk;
        denominator = A_bar_t(~Wk, :)*p;
        idx =denominator < 0;
        nominator = nominator(idx);
        denominator = denominator(idx);
        alphas = nominator ./ denominator;
        j = alphas == min(alphas);
        alpha = alphas(j);

        if alpha < 1
            xk = xk + alpha*p;
            Wk(j) = 1; 
        else
            xk = xk + p;
            Wk = Wk;
        end
        
    end 

    all_xk = [all_xk xk];

    if verbose

        fprintf("ITERATION: %0.0f\n\n", counter);
        fprintf("Working set (Wk) in iteration: %0.0f\n", counter);
        disp(all_constraints(Wk));  
    
        fprintf("Step (p) in iteration: %0.0f\n", counter);
        disp(p);
    
        fprintf("x (xk) in iteration: %0.0f\n", counter);
        disp(xk);
    
        fprintf("Alpha in iteration: %0.0f\n", counter);
        disp(alpha);

    end

    counter = counter + 1;

end

end