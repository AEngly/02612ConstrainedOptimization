function [primal_final, dual_final, solverInformation] = SQPLineSearchDampedBFGS(fun,x0,xlb,xub,clb,cub,nonlcon,options)

    % Auxiliary variables
    iter = 0;
    maxit = options.maxit;
    BFGS = options.BFGS;
    stepSolver = options.stepSolver;
    lineSearch = options.lineSearch;
    converged = false;
    rho = 0.9; % (has to be between 0 and 1)
    muk = options.l1Penalty;
    l1Penalty = options.l1Penalty;
    mu_margin = 1;

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
    [f, fGrad] = fun(xk);

    % Save number of variables
    n = length(x0);
    miq = length(c);
    meq = length(ceq);

    % Save number of bounds
    m_cub = length(cub);
    m_clb = length(clb);
    m_xub = length(xub);
    m_xlb = length(xlb);

    % Specify large number and small number to cope with missing bounds on
    % constraints (only if partially specified)
    mBIG = 10^6;
    mSMALL = -10^6;

    % Append to ensure that all constraints get bounds if only partially
    % specified

    if m_cub > 0 && miq > m_cub
        cub = [cub; mBIG*ones(miq-m_cub,1)];
        m_cub = miq;
    end
    if m_clb > 0 && miq > m_clb
        clb = [clb; mSMALL*ones(miq-m_clb,1)];
        m_clb = miq;
    end
    if m_xub > 0 && n > m_xub
        xub = [xub; mBIG*ones(n-m_xub,1)];
        m_xub = n;
    end
    if m_xlb > 0 && n > m_xlb
        xlb = [xlb; mSMALL*ones(n-m_xlb,1)];
        m_xlb = n;
    end

    % Require initial dual variables
    lk = ones(4*meq+2*m_cub+2*m_clb+2*m_xub+2*m_xlb,1);

    % Initialize positive definite matrix for damped BFGS
    Bk = eye(n);

    % Compute steps until convergence test is satisfied
    while ~converged && (iter < maxit)

        % 1) Solve Equation 18.11 on page 533 in N&W
        % We start by formulating the nonlinear constraints (but without the x's for quadprog).
        % For quadprog, we need the form Aeq*x = beq and A*x <= b.
        % We formulate it initially as Ceq + deq = 0 and C + d >= 0.
        % Then clearly Aeq = Ceq, beq = -deq, A = -C, and b = d.
        % In addition, we add infeasible handling.

        if stepSolver == "quadprog"  
            
            Ceq = [];
            deq = [];
            if ~isempty(ceq)
                Ceq = [Ceq; GCeq' -eye(meq) eye(meq) zeros(meq,m_cub+m_clb+m_xub+m_xlb)];
                deq = [deq; ceq];
            end

            C = [];
            d = [];
            if ~isempty(cub)
                blockRow = [-GC'  zeros(miq,meq) zeros(miq,meq) eye(m_cub) zeros(miq,m_clb) zeros(miq,m_xub) zeros(miq,m_xlb)];
                fprintf("Display blockrow in 1 if:\n")
                disp(blockRow);
                C = [C; blockRow];
                d = [d; cub-c];
            end
            if ~isempty(clb)
                blockRow = [GC'  zeros(miq,meq) zeros(miq,meq) zeros(miq,m_cub) eye(m_clb) zeros(miq,m_xub) zeros(miq,m_xlb)];
                fprintf("Display blockrow in 2 if:\n")
                disp(blockRow);
                C = [C; blockRow];
                d = [d; c-clb];
            end
            if ~isempty(xub)
                blockRow = [-eye(m_xub)  zeros(m_xub,meq) zeros(m_xub,meq) zeros(m_xub,miq) zeros(m_xub) eye(m_xub) zeros(m_xub)];
                fprintf("Display blockrow in 3 if:\n")
                disp(blockRow);
                C = [C; blockRow];
                d = [d; xub - xk];
            end
            if ~isempty(xlb)
                fprintf("Display blockrow in 4 if:\n")
                blockRow = [eye(m_xlb)  zeros(m_xlb,meq) zeros(m_xlb,meq) zeros(m_xlb,miq) zeros(m_xlb) zeros(m_xlb) eye(m_xlb)];
                disp(blockRow);
                C = [C; blockRow];
                d = [d; xk - xlb];
            end
            C = [C; zeros(2*meq+m_cub+m_clb+m_xub+m_xlb,n) eye(2*meq+m_cub+m_clb+m_xub+m_xlb)];
            d = [d; zeros(2*meq+m_cub+m_clb+m_xub+m_xlb,1)];

        end

        % Then we create F and G as stated in the report
        F = [fGrad; l1Penalty*ones(2*meq+m_cub+m_clb+m_xub+m_xlb,1)];
        mG = 2*meq+m_cub+m_clb+m_xub+m_xlb;
        G = [Bk zeros(n,mG); zeros(mG,n) zeros(mG,mG)];

        fprintf("Display Aeq:\n")
        disp(Ceq);
        fprintf("Display beq:\n")
        disp(-deq);
        fprintf("Display A:\n")
        disp(-C);
        fprintf("Display b:\n")
        disp(d);

        % Translate to quadprog syntax
        Aeq = Ceq;
        beq = -deq;
        A = -C;
        b = d;

        % 2) Solving local QP subproblem
        [steps,fval,exitflag,output,dual] = quadprog(G,F,A,b,Aeq,beq,[],[],[],optionsQP);
        if ~(exitflag == 1)
            fprintf("Step could not be computed!");
        end
        pk = steps(1:n);
        yk1 = dual.eqlin;
        zk1 = dual.ineqlin;
        lk1 = [yk1;zk1];

        % Compute step size of lambda
        disp(size(A));
        disp(size(lk));
        disp(size(lk1));
        p_lambda = lk1 - lk;

        if lineSearch == "slides"

            % Initialize alpha to 1
            alpha = 1;

            % Compute constraints and function with full step (part 1)
            fun_alpha = fun(xk + alpha*pk);
            [c,ceq] = nonlcon(xk);
            [c_alpha,ceq_alpha] = nonlcon(xk + alpha*pk);

            % Compute constraints and function with full step (part 2)
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
            if ~isempty(xub)
                c_alpha_ineq = [c_alpha_ineq; xub-(xk + alpha*pk)];
                c_ineq = [c_ineq; xub-xk];
            end
            if ~isempty(xlb)
                c_alpha_ineq = [c_alpha_ineq; (xk + alpha*pk) - xlb];
                c_ineq = [c_ineq; xk-xlb];
            end

            % Compute constraints and function with full step (part 2)
            ck = [ck; c_ineq];

            % Compute penalty in each iteration
            %muk = (fGrad'*pk + (1/2)*pk'*Bk*pk)/((1-rho)*norm(ck,1)) + mu_margin;
            muk = max(max(abs(lk),1/2*(muk+abs(lk))));

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
                if ~isempty(xub)
                    c_alpha_ineq = [c_alpha_ineq; xub-xk];
                end
                if ~isempty(xlb)
                    c_alpha_ineq = [c_alpha_ineq; xk-xlb];
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

            fprintf("Printing value of alpha:\n")
            disp(alpha);

            % Accept step size
            xk = xk + alpha*pk;
            lk = lk + alpha*p_lambda;

            % Store step
            solverInformation.stepSequence = [solverInformation.stepSequence xk];

        else
            % Accept step size
            alpha = 1;
            xk = xk + alpha*pk;
            lk = lk + alpha*p_lambda;
        end

        disp(xk);

        % In order to compute the lagrangian, we need original xk and new
        % iterate
        [c_xk,ceq_xk,GC_xk,GCeq_xk] = nonlcon(xk-alpha*pk);
        [c_xk1,ceq_xk1,GC_xk1,GCeq_xk1] = nonlcon(xk);
        [fun_xk, fun_xkGrad] = fun(xk - alpha*pk);
        [fun_xk1, fun_xk1Grad] = fun(xk);

        % Then we construct the appropriate row vector
        mBC = m_cub + m_clb + m_xlb + m_xub;
        gGrad = [-GC_xk GC_xk -eye(n) eye(n) zeros(n,mBC)];
        LGrad = fun_xkGrad - gGrad*lk1;
        gGrad1 = [-GC_xk1 GC_xk1 -eye(m_xub) eye(m_xlb) zeros(n,mBC)];             
        LGrad1 = fun_xk1Grad - gGrad1*lk1;
        
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
        fprintf("\nFinished interation %d!\n", iter);
        iter = iter + 1;

        % Test if converged
        if norm(pk*alpha,2) < options.convergenceRequirement

            converged = true;
            fprintf("Algorithm converged succesfully.\n");

            % Save important information
            solverInformation.iterations = iter;

            % Set results
            primal_final = xk;
            dual_final = lk1;

        end

    end

end