function [primal_final, dual_final, solverInformation] = SQPSimpleDampedBFGS(fun,x0,xlb,xub,clb,cub,nonlcon,options)

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
    epsilon = options.convergenceRequirement;

    % Variables for Line Search
    tau = 0.9; % 
    armijoC = 0.1;

    % Options for quadprog
    optionsQP = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');

    % Compute necessary values in first step
    [c,ceq,GC,GCeq] = nonlcon(x0);
    [f, fGrad] = fun(x0);

    % Save number of variables
    n = length(x0);
    m = length(c);

    % Require starting condition
    xk = x0;
    lk = ones(m,1);

    % Create a struct to store information
    solverInformation = struct();
    solverInformation.primalSequence = xk;
    solverInformation.dualSequence = lk;
    solverInformation.stepSequence = [];
    solverInformation.stepLambdaSequence = [];

    % Specify large number and small number to cope with missing bounds on
    % constraints (only if partially specified)
    mBIG = 10^6;
    mSMALL = -10^6;

    % Initialize positive definite matrix for damped BFGS
    Bk = 5*eye(n);

    % Compute steps until convergence test is satisfied
    while ~converged && (iter < maxit)

        fprintf("Starting iteration %d:\n\n", iter)

        % Start by constructing the problem for quadprog
        C = GC';
        d = c;

        % Then construct H,f,A,b for quadprog
        H = Bk;
        f = fGrad;
        A = -C;
        b = d;

        % Quadprog
        [pk,fval,exitflag,output,dual] = quadprog(H,f,A,b,[],[],[],[],[],optionsQP);

        % Check if it succeded
        if ~(exitflag == 1)
            fprintf("Step could not be computed!");
        end

        % Then increase xk (solution to next multipliers is given directly from the above)
        lk1 = dual.ineqlin;
        p_lambda = lk1 - lk;
        xk = xk + pk;
        lk = lk1;

        % Store new step
        solverInformation.primalSequence = [solverInformation.primalSequence xk];
        solverInformation.dualSequence = [solverInformation.dualSequence lk];
        solverInformation.stepSequence = [solverInformation.stepSequence pk];
        solverInformation.stepLambdaSequence = [solverInformation.stepLambdaSequence p_lambda];

        disp(pk);
        disp(xk - pk);
        disp(xk);
        disp(lk);
        disp(p_lambda);
        disp(lk1);

        % Update Hessian with damped BFGS
                % In order to compute the lagrangian, we need original xk and new
        % iterate
        [c_step,ceq_step,GC_step,GCeq_step] = nonlcon(xk);
        [fStep, fStepGrad] = fun(xk);

        % Then we construct the appropriate row vector
        LGrad = fGrad - GC*lk1;           
        LGrad1 = fStepGrad - GC_step*lk1;
        
        % Then get pk and qk
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

        % Recompute function values
        [c,ceq,GC,GCeq] = nonlcon(xk);
        [f, fGrad] = fun(xk);

        % Check for convergence
        if norm(fGrad - GC*lk, 'Inf') < epsilon
            converged = true;
        end

        iter = iter + 1;

    end

    primal_final = xk;
    dual_final = lk;