function [primal_final, dual_final, solverInformation] = SQPLineSearchBFGS(fun,x0,xlb,xub,clb,cub,nonlcon,options)

    % ---------------- DESCRIPTION --------------
    %
    % Name: SQPLineSearchBFGS
    % Type: SQP procedure with Line Search, BFGS and infeasibility handling
    %
    % Note: This implementation used 'quadprog' to solve subproblems.
    %
    % Problem structure:
    %           min     f(x)
    %            x
    %           s.t.    h(x) = 0
    %                   gl  <= g(x) <= gu
    %                   xl  <= x    <= xu
    %
    % Syntax: [primal_final, dual_final, solverInformation] = SQPLineSearchBFGS(fun,x0,xlb,xub,clb,cub,nonlcon,options)
    %
    %
    % Created: 21.05.2023
    % Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
    %          Compute, Technical University of Denmark
    %
    
    % ---------------- IMPLEMENTATION --------------

    % Auxiliary variables
    iter = 0;
    maxit = options.maxit;
    converged = false;
    epsilon = options.convergenceRequirement;

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

    % Initialize penalties for Powell exact merit
    eqPenalty = 0; % NOT IMPLEMENTED CURRENTLY
    ineqPenalty = lk+1;
    slackPenalty = options.infeasibilityPenalty;

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
        H = [Bk zeros(n,m); zeros(n,n) zeros(n,m)];
        f = [fGrad; slackPenalty*ones(m,1)];
        A = [-C -eye(m); zeros(m,m) -eye(m)];
        b = [d; zeros(m,1)];

        % Quadprog
        disp(H);
        [primal,fval,exitflag,output,dual] = quadprog(H,f,A,b,[],[],[],[],[],optionsQP);

        % Extract the right variables
        pk = primal(1:n);
        t = primal(n+1:end);
        disp(t);

        % Check if it succeded
        if ~(exitflag == 1)
            fprintf("Step could not be computed!");
        end

        % Then we complete a line search before proceeding with the step
        eqPenalty = 0; % NOT IMPLEMENTED CURRENTLY
        ineqPenalty = max(abs(lk),0.5*(ineqPenalty + abs(lk)));
        alpha = LineSearch(xk,pk,eqPenalty,ineqPenalty,fun,nonlcon);
        disp(alpha);

        % Then increase xk (solution to next multipliers is given directly from the above)
        lk1 = dual.ineqlin(1:m); % Only take the duals corresponding to the actual constraints supplied to the program
        p_lambda = lk1 - lk;
        xk = xk + alpha*pk;
        lk = lk + alpha*p_lambda;

        % Store new step
        solverInformation.primalSequence = [solverInformation.primalSequence xk];
        solverInformation.dualSequence = [solverInformation.dualSequence lk];
        solverInformation.stepSequence = [solverInformation.stepSequence pk];
        solverInformation.stepLambdaSequence = [solverInformation.stepLambdaSequence p_lambda];

        % Update Hessian with damped BFGS
        % In order to compute the lagrangian, we need original xk and new
        % iterate
        [c_step,ceq_step,GC_step,GCeq_step] = nonlcon(xk);
        [fStep, fStepGrad] = fun(xk);

        % Then we construct the appropriate row vector
        LGrad = fGrad - GC*lk;           
        LGrad1 = fStepGrad - GC_step*lk;
        
        % Then get pk and qk
        pk = pk*alpha;
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
        fprintf("\nChecking convergence:\n")
        disp(norm(fGrad - GC*lk, 'Inf'));
        if norm(fGrad - GC*lk, 'Inf') < epsilon
            converged = true;
        end

        iter = iter + 1;

    end

    primal_final = xk;
    dual_final = lk;

end

function [val] = PowellExact(alpha,xk,pk,lambda,mu,fun,nonlcon)
    val = fun(xk + alpha*pk) + mu'*abs(min(0,nonlcon(xk + alpha*pk)));
end

function [val] = PowellExactDerivative(alpha,xk,pk,lambda,mu,fun,nonlcon)
    [f, fGrad] = fun(xk + alpha*pk);
    val = fGrad'*pk + mu'*abs(min(0,nonlcon(xk + alpha*pk)));
end

function [newAlpha] = LineSearch(xk,pk,eqPenalty,ineqPenalty,fun,nonlcon)

    [f, fGrad] = fun(xk);
    [c_step,ceq_step,GC_step,GCeq_step] = nonlcon(xk); 

    alpha = 1;
    acceptStep = false;

    % Compress functions
    phi = @(alpha) PowellExact(alpha,xk,pk,eqPenalty,ineqPenalty,fun,nonlcon);
    dphi = @(alpha) PowellExactDerivative(alpha,xk,pk,eqPenalty,ineqPenalty,fun,nonlcon);

    % Evaluate merit functions
    c = phi(0);
    b = dphi(0);

    while ~acceptStep

        % Take step
        x = xk + alpha*pk;

        % Evaluate functions
        [fx, fGradx] = fun(x);
        [c_x,ceq_x,GC_x,GCeq_x] = nonlcon(x);

        if phi(alpha) <= phi(0) + 0.1*dphi(0)*alpha
            acceptStep = true;
        else

            % Evaluate terms
            a = (phi(alpha) - (c+b*alpha))/(alpha^2);
            alpha_min = -b/(2*a);
            
            % Update alpha
            alpha = min(0.9*alpha,max(alpha_min,0.1*alpha));

        end

    end
    
    newAlpha = alpha;

end