function [x_final, solverInformation] = SQPLineSearchDampedBFGS(fun,x0,lb,ub,clb,cub,nonlcon,options)

    % Auxiliary variables
    iter = 0;
    maxit = options.maxit;
    BFGS = options.BFGS;
    stepSolver = options.stepSolver;
    lineSearch = options.lineSearch;
    converged = false;
    sigma = 1;
    rho = 0.9; % (has to be between 0 and 1)
    l1_penalty = options.l1Penalty;
    mu_margin = 1e-3;

    % Create a struct to store information
    solverInformation = struct();
    solverInformation.stepSequence = x0;

    % Variables for Line Search
    tau = 0.9; % 
    armijoC = 0.1;

    % Algorithm used for quadprog
    if stepSolver == "quadprog"
        %optionsQP = optimoptions('quadprog','Algorithm','active-set','Display','off');
        optionsQP = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
    end

    % Require starting condition
    xk = x0;

    % Compute necessary values in first step
    [c,ceq,GC,GCeq] = nonlcon(xk);
    [f, fGrad, H] = fun(xk);

    % Save number of variables
    n = length(x0);
    m = length(c);
    meq = length(ceq);

    % Require initial dual variables
    lk = zeros(meq+2*m,1);

    % Initialize positive definite matrix for damped BFGS
    Bk = eye(n);

    % Compute steps until convergence test is satisfied
    while ~converged && (iter < maxit)

        % 1) Solve Equation 5.18 (page X in N&W)
        % We start by formulating the nonlinear constraints (but without the x's for quadprog).
        % For quadprog, we need the form Aeq*x = beq and A*x <= b.
        % We formulate it initially as Ceq + deq = 0 and C + d >= 0.
        % Then clearly Aeq = Ceq, beq = -deq, A = -C, and b = d.
        % In addition, we add infeasible handling.

        if stepSolver == "quadprog"  
            % Create matrices
            Ceq = [];
            deq = [];
            if ~isempty(ceq)
                Ceq = [Ceq; GCeq' -eye(meq) eye(meq) zeros(meq,2*m)];
                deq = [deq; ceq];
            end
            C = [];
            d = [];
            if ~isempty(cub)
                C = [C; -GC' zeros(m,2*meq+2*m)];
                d = [d; cub-c];
            end
            if ~isempty(clb)
                C = [C; GC' zeros(m,2*meq+2*m)];
                d = [d; c-clb];
            end

            % Extend lower and upper bound to quadprog
            lb_new = [lb; zeros(2*meq+2*m,1)];
            ub_new = [ub; inf*ones(2*meq+2*m,1)];

        end

        % This section should be uncommented when it is run with another
        % solver than quadprog
        %if ~isempty(ub)
        %    C = [C; -eye(length(xk)) zeros(m,2*meq+2*m+2*n)];
        %    d = [d; ub-xk];
        %end
        %if ~isempty(lb)
        %    C = [C; eye(length(xk)) zeros(m,2*meq+2*m+2*n)];
        %    d = [d; xk-lb];
        %end
        %C = [C; zeros(n) zeros(m,2*meq+2*m+2*n)];
        %d = [d; zeros(n,1)];

        % Then we create F and G as stated in the report
        F = [fGrad; l1_penalty*ones(2*meq + 2*m,1)];
        G = [Bk zeros(n,2*meq+2*m); zeros(2*meq+2*m,n) zeros(2*meq+2*m,2*meq+2*m)];

        % Translate to quadprog syntax
        Aeq = Ceq;
        beq = -deq;
        A = -C;
        b = d;

        % 2) Solving local QP subproblem
        [steps,fval,exitflag,output,dual] = quadprog(G,F,A,b,Aeq,beq,lb_new,ub_new,x0,optionsQP);
        pk = steps(1:n);
        yk1 = dual.eqlin;
        zk1 = dual.ineqlin;
        lk1 = [yk1;zk1];

        % Compute step size of lambda
        p_lambda = lk1 - lk;

        if lineSearch == "second-order"

            % 3) Use Line Search from slide 6/11 ("Lecture 09B_SQP")
            alpha = 1;
            fun_alpha = fun(xk + alpha*pk);
            [c_alpha,ceq_alpha] = nonlcon(xk + alpha*pk);
    
            % Powell's update of penalty parameters
            lambda = abs(yk1);
            mu = abs(zk1);
    
            % Compute merit functions and its directional derivative
            phi_alpha = fun_alpha + lambda'*abs(ceq_alpha) + mu'*abs(min(0,c_alpha));
            phi_0 = f + lambda'*abs(ceq) + mu'*abs(c);
            d_phi_0 = fGrad'*pk - lambda'*abs(ceq) + mu'*abs(min(0,c_alpha));
    
            while ~(phi_alpha <= phi_0 + armijoC*d_phi_0*alpha)
    
                % Update alpha
                a = (phi_alpha - (phi_0 + d_phi_0*alpha))/(alpha^2);
                alpha_min = -phi_0/(2*a);
                alpha = min(tau*alpha, max(alpha_min,armijoC*alpha));
    
                % Repeat first step
                fun_alpha = fun(xk + alpha*pk);
                [c_alpha,ceq_alpha] = nonlcon(xk + alpha*pk);
    
                % Compute merit functions and its directional derivative
                phi_alpha = fun_alpha + lambda'*abs(ceq_alpha) + mu'*abs(min(0,c_alpha));
                phi_0 = f + lambda'*abs(ceq) + mu'*abs(c);
                d_phi_0 = fGrad'*pk - lambda'*abs(ceq) + mu'*abs(min(0,c_alpha));
    
            end
    
            % Accept step
            xk = xk + alpha*pk;

        elseif lineSearch == "slides"

            % Initialize alpha to 1
            alpha = 1;

            % Compute constraints and function with full step
            fun_alpha = fun(xk + alpha*pk);
            [c,ceq] = nonlcon(xk);
            [c_alpha,ceq_alpha] = nonlcon(xk + alpha*pk);

            % Take care of bound constraints
            c_alpha_ineq = [];
            c_ineq = [];
            ck = ceq;
            if ~isempty(cub)
                c_alpha_ineq = [c_alpha_ineq; cub-c_alpha];
                c_ineq = [c_ineq; cub-c];
            end
            if ~isempty(clb)
                c_alpha_ineq = [c_alpha_ineq; c_alpha-clb];
                c_ineq = [c_ineq; c-clb];
            end

            % Compute mu in iteration k
            ck = [ck; c_ineq];
            muk = (fGrad'*pk + (1/2)*pk'*Bk*pk)/((1-rho)*norm(ck,1)) + mu_margin;

            % Compute merit functions and its directional derivative
            term_ceq_alpha = muk*sum(abs(ceq_alpha));
            term_c_alpha = muk*sum(abs(min(0,c_alpha_ineq))); 
            phi_alpha = fun_alpha;
            phi_0 = f;
            d_phi_0 = fGrad'*pk;

            if ~isempty(term_ceq_alpha)
                phi_alpha = phi_alpha + term_ceq_alpha;
                phi_0 = phi_0 + muk*sum(abs(ceq));
                d_phi_0 = d_phi_0 - muk*sum(abs(ceq));
            end
            if ~isempty(term_c_alpha)
                phi_alpha = phi_alpha + term_c_alpha;
                phi_0 = phi_0 + muk*sum(abs(min(0,c_ineq)));
                d_phi_0 = d_phi_0 - muk*sum(abs(min(0,c_ineq)));
            end

            while phi_alpha > phi_0 + armijoC*alpha*d_phi_0

                % Decrease alpha accordingly
                alpha = alpha*tau;

                % Recompute constraints and objective with new alpha
                fun_alpha = fun(xk + alpha*pk);
                [c_alpha,ceq_alpha] = nonlcon(xk + alpha*pk);
    
                % Take care of bound constraints
                c_alpha_ineq = [];
                if ~isempty(cub)
                    c_alpha_ineq = [c_alpha_ineq; cub-c_alpha];
                end
                if ~isempty(clb)
                    c_alpha_ineq = [c_alpha_ineq; c_alpha-clb];
                end
    
                % Compute merit functions and its directional derivative
                term_ceq_alpha = muk*sum(abs(ceq_alpha));
                term_c_alpha = muk*sum(abs(min(0,c_alpha_ineq))); 
                phi_alpha = fun_alpha;
    
                if ~isempty(term_ceq_alpha)
                    phi_alpha = phi_alpha + term_ceq_alpha;
                end
                if ~isempty(term_c_alpha)
                    phi_alpha = phi_alpha + term_c_alpha;
                end
            
            end

            % Accept step size
            xk = xk + alpha*pk;
            lk = lk + alpha*p_lambda;

            % Store step
            solverInformation.stepSequence = [solverInformation.stepSequence xk];

        end

        % In order to compute the lagrangian, we need the full setup
        [c_xk,ceq_xk,GC_xk,GCeq_xk] = nonlcon(xk-alpha*pk);
        [c_xk1,ceq_xk1,GC_xk1,GCeq_xk1] = nonlcon(xk);
        [fun_xk, fun_xkGrad] = fun(xk - alpha*pk);
        [fun_xk1, fun_xk1Grad] = fun(xk);

        % Then we construct the appropriate row vector
        gGrad = [-GC_xk GC_xk];
        LGrad = fun_xkGrad - gGrad*lk;
        gGrad1 = [-GC_xk1 GC_xk1];             
        LGrad1 = fun_xk1Grad - gGrad1*lk;
        
        % Then get pk and qk
        sk = alpha*pk;
        qk = LGrad1 - LGrad;

        % Then be can use the modified BFGS procedure
        if pk'*qk >= 0.2*pk'*Bk*pk
            theta = 1;
        else
            theta = (0.8*pk'*Bk*pk)/(pk'*Bk*pk - pk'*qk);
        end

        % Then we can find rk
        rk = theta*qk + (1 - theta)*Bk*pk;

        % Then we can update the matrix
        Bk = Bk + (rk*rk')/(pk'*rk) - ((Bk*pk)*(Bk*pk)')/(pk'*Bk*pk);

        % Iteration completed
        iter = iter + 1;

        % Test if converged
        if norm(pk*alpha,2) < options.convergenceRequirement

            converged = true;
            fprintf("Algorithm converged succesfully.\n");

            % Save important information
            solverInformation.iterations = iter;

            % Set results
            x_final = xk;

        end

    end

end