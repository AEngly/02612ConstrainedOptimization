
% ---------------- DESCRIPTION --------------
%
% Name: SQPTrustRegion
% Type: Sequential l1 Quadratic Programming Algorithm
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
% Syntax: MISSING
%
%
% Created: 09.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------

function [primal_final, dual_final, solverInformation] = SQPTrustRegion(fun,x0,xlb,xub,clb,cub,nonlcon,options)

    % Auxiliary variables
    iter = 0;
    maxit = options.maxit;
    BFGS = options.BFGS;
    stepSolver = options.stepSolver;
    lineSearch = options.lineSearch;
    trustRegion = options.trustRegion;
    converged = false;
    muk = options.l1Penalty;
    acceptanceMargin = options.acceptanceMargin;
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

    % Initialize penalties for Powell exact merit
    eqPenalty = 0; % NOT IMPLEMENTED CURRENTLY
    ineqPenalty = lk;
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
    Bk = eye(n);

    % Compute steps until convergence test is satisfied
    while ~converged && (iter < maxit)

        fprintf("Starting iteration %d:\n\n", iter)

        % Start by constructing the problem for quadprog
        C = GC';
        d = c;

        % Then construct H,f,A,b for quadprog
        H = [Bk zeros(n,m+1); zeros(m+1,n) zeros(m+1,m+1)];
        f = [fGrad; slackPenalty*ones(m,1); zeros(1,1)];
        A = [-C -eye(m) zeros(m,1); % Lower bound
             eye(m) zeros(m,m) -ones(m,1); % Formalizing norm part 1
             -eye(m) zeros(m,m) -ones(m,1); % Formalizing norm part 2
             zeros(m,m) -eye(m) zeros(m,1); % Nonnegativity for all t's
             zeros(1,2*m) ones(1,1)]; % Nonnegativity for s
        b = [d; zeros(3*m,1); trustRegion];

        % Quadprog
        [primal,fval,exitflag,output,dual] = quadprog(H,f,A,b,[],[],[],[],[],optionsQP);

        % Extract primal variables
        pk = primal(1:n);
        t = primal(n+1:end-1);
        s = primal(end);

        % Check if it succeded
        if ~(exitflag == 1)
            fprintf("Step could not be computed!");
        end

        % Update penalty parameter muk (follow algorithm 18.5)
        % This part is not implemented yet.
        
        % Determine whether to accept the step
        ared_k = phi_1(xk,muk,fun,nonlcon) - phi_1(xk+pk,muk,fun,nonlcon);
        pred_k = phi_1(xk,muk,fun,nonlcon) - q_mu(xk,pk,Bk,muk,fun,nonlcon);
        rho_k = ared_k/pred_k;

        % Compute scalar to change trust region
        gamma = min(max((2*rho_k-1)^3+1,0.25),2);

        if rho_k > acceptanceMargin

            % Update step
            xk = xk + pk;
            lk1 = dual.ineqlin(1:m);
            p_lambda = lk1 - lk;
            lk = lk + p_lambda;

             % Update trust region
            trustRegion = gamma*trustRegion;

            % Update Hessian with damped BFGS
            % In order to compute the lagrangian, we need original xk and new
            % iterate
            [c_step,ceq_step,GC_step,GCeq_step] = nonlcon(xk);
            [fStep, fStepGrad] = fun(xk);
    
            % Then we construct the appropriate row vector
            LGrad = fGrad - GC*lk;           
            LGrad1 = fStepGrad - GC_step*lk;
            
            % Then get pk and qk
            pk = pk;
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

        else

            % Keep the same step
            xk = xk;
            p_lambda = lk - lk;
            lk = lk + p_lambda;

            % Update trust region
            trustRegion = gamma*vecnorm(pk,'Inf');

            % Store step
            pk = zeros(n,1);

        end

        % Store information from iteration
        solverInformation.primalSequence = [solverInformation.primalSequence xk];
        solverInformation.dualSequence = [solverInformation.dualSequence lk];
        solverInformation.stepLambdaSequence = [solverInformation.stepLambdaSequence p_lambda];
        solverInformation.stepSequence = [solverInformation.stepSequence pk];

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
    dual_final = lk(1:n);

end

function [val] = phi_1(xk,mu,fun,nonlcon)
    [ck,ceqk,GCk,GCeqk] = nonlcon(xk);
    eqTerm = abs(0);
    ineqTerm = max(0,-ck);
    val = fun(xk) + mu*sum(eqTerm) + mu*sum(ineqTerm);
end

function [val] = q_mu(xk,pk,Bk,mu,fun,nonlcon)
    [fk, fGradk] = fun(xk);
    [ck,ceqk,GCk,GCeqk] = nonlcon(xk);
    ineqTerm = max(0,-(ck + GCk'*pk));
    eqTerm = abs(0);
    val = fk + fGradk'*pk + 0.5*pk'*Bk*pk + mu*sum(eqTerm) + mu*sum(ineqTerm);
end
