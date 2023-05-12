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

function [x_final, solverInformation] = SQPTrustRegion(fun,x0,lb,ub,clb,cub,nonlcon,options)

    % Read option choices
    maxIter = options.maxIter;
    trustRegion = options.trustRegion;
    epsilon = 10^(-options.precision);
    eta = options.eta;
    mu = options.l1Penalty;

    % Save information variables
    k = 1;

    % Compute relevant function values
    xk = x0;
    [f,df] = fun(xk);
    [c,ceq,GC,GCeq] = nonlcon(xk);

    % Compute relevant dimension
    n = length(xk);
    miq = size(c,1);
    meq = size(ceq,1);

    % Define initial dual variables as vector  of 0s
    zk = zeros(6*n+4*miq+3*meq+2,1);

    % Prepare matrix for damped BFGS
    Bk = eye(n);

    % Define options for quadprog
    %options = optimset('Display', 'off');
    options = optimset();

    % Define additonals
    converged = false;
    
    % Start for loop
    while ~converged && k <= maxIter
        
        % PART 1: Solve subproblem
        % 1.1) Define matrices for objective
        H = zeros(3*n + 2*meq + 2*miq+1, 3*n + 2*meq + 2*miq+1);
        H(1:n,1:n) = Bk;
        d = [df; mu*ones(2*meq,1);mu*ones(2*miq+2*n,1);0];

        % 1.2) Define matrices for equality constraints
        Aeq = [GCeq' eye(meq) -eye(meq) zeros(meq,2*n + 2*miq)];
        beq = ceq;

        % 1.3) Define matrices for inequality constraints
        A = zeros(6*n+4*miq+2*meq+2, 3*n+2*meq+2*miq+1);
        A(1:2*miq+4*n,1:n) = [-GC'; GC; -eye(n); eye(n); -eye(n); eye(n)];
        A1 = blkdiag(eye(2*meq),eye(miq),eye(miq),eye(n),eye(n),[1;-1]);
        A2 = blkdiag(eye(miq),eye(miq),eye(n),eye(n), [ones(n,1); -ones(n,1)]);
        A(2*miq+4*n+1:end,n+1:end) = A1;
        A(1:2*miq+4*n,n+2*meq+1:end) = A2;

        b = zeros(6*n+4*miq+2*meq+2, 1);
        b(1:2*miq+2*n,1) = [cub - c; c - clb; ub - xk; xk - lb];
        b(1:2*miq+2*n,1) = [cub - c; c - clb; ub - xk; xk - lb];
        b(end-1:end,1) = [trustRegion;0];

        % The system is now defined as
        % Aeq + beq = 0
        % A   + b  >= 0
        % In quadprog, the notation is as follows:
        % Aeq = beq
        % A   => b

        % 1.4) Solve system
        [primal,~,~,~,dual] = quadprog(H,d,-A,b,Aeq,-beq,[],[],[],options);
        
        % 1.5) Separate primal variables
        pk = primal(1:n);
        disp(pk);
        w = primal(n+1:n+meq);
        v = primal(n+meq+1:n+meq*2);
        t = primal(n+meq*2+1:end-1);
        s = primal(end);

        % PART 2: Update penalty parameter (Powell update -> slides 9A page 22)
        lambdaHat = dual.ineqlin(1:2*n+2*miq);
        lambdaNorm = vecnorm(lambdaHat,'Inf');
        mu = max(1/2*(mu+lambdaNorm),lambdaNorm);
        
        % PART 3: Accept or reject step, then update trust region
        % Compute constraint and append to get right form
        [c,ceq,GC,GCeq] = nonlcon(xk);
        cFull = zeros(2*n+2*miq, 1);
        cFull(1:2*miq+2*n,1) = [cub - c; c - clb; ub - xk; xk - lb];
        ceqFull = zeros(meq, 1);
        ceqFull(1:meq,1) = ceq;  
        GCFull = [-GC GC -eye(n) eye(n)];
        GCeqFull = GCeq;
        
        % Repeat where we "take" step
        [c,ceq,GC,GCeq] = nonlcon(xk+pk);
        cFullp = zeros(2*n+2*miq, 1);
        cFullp(1:2*miq+2*n,1) = [cub - c; c - clb; ub - (xk+pk); (xk+pk) - lb];
        ceqFullp = zeros(meq, 1);
        ceqFullp(1:meq,1) = ceq;   
        GCFullp = [-GC GC -eye(n) eye(n)];
        GCeqFullp = GCeq;

        % Compute objective of l1 form (q(pk))
        qP1 = df'*pk;
        qP2 = 1/2*pk'*Bk*pk;
        if ~isempty(GCeqFull)
            qP3 = sum(mu*abs(ceqFull+GCeqFull'*pk));
        else
            qP3 = 0;
        end
        if ~isempty(GCFull)
            qP4 = sum(mu*max(0,-(cFull+GCFull'*pk)));
        else
            qP4 = 0;
        end
        qP = f+qP1+qP2+qP3+qP4;

        % Compute objective of l1 form (q(0))
        q01 = sum(mu*abs(ceqFull));
        q02 = sum(mu*max(0,-cFull));
        q0 = f+q01+q02;

        % Compute objective of l1 form (phi1(0) and phi1(pk))
        phi1 = q0;
        phi1P = fun(xk + pk) + sum(mu*abs(ceqFull)) + sum(mu*max(0,-cFullp));
        
        % The we can compute the acceptance ratio
        rho = (phi1-phi1P)/(q0-qP);

        % If trust region is accepted
        if rho>eta

            % Update trust region
            gamma = min(max((2*rho-1)^3+1,0.25),2);
            trustRegion = gamma*trustRegion;

            % Compute update of Bk
            [~,dfP] = fun(xk + pk);
            if ~isempty(GCeqFullp')
                LGradP1 = GCeqFullp*dual.eqlin;
            else
                LGradP1 = 0;
            end
            if ~isempty(GCFullp)
                LGradP2 = GCFullp*dual.ineqlin(1:2*n+2*miq);
            else
                LGradP2 = 0;
            end
            if ~isempty(GCeqFull)
                LGrad1 = GCeqFull*dual.eqlin;
            else
                LGrad1 = 0;
            end
            if ~isempty(GCFull)
                LGrad2 = GCFull*dual.ineqlin(1:2*n+2*miq);
            else
                LGrad2 = 0;
            end

            LGradP = dfP - LGradP1 - LGradP2;
            LGrad = df - LGrad1 - LGrad2;

            % Update the current iterate
            zk = lambdaHat;
            xk = xk + pk;

            % Update values for next iteration
            [f,df] = fun(xk);
            [c,~] = nonlcon(xk);

            % Then get qk
            qk = LGradP - LGrad;
            
            % Then be can use the modified BFGS procedure
            if pk'*qk >= 0.2*pk'*Bk*pk
                theta = 1;
            else
                denom = 0.8*pk'*Bk*pk;
                nom = pk'*Bk*pk - pk'*qk;
                theta = denom/nom;
            end
    
            % Then we can find rk
            rk = theta*qk + (1 - theta)*Bk*pk;
    
            % Then we can update the matrix
            Bk = Bk + (rk*rk')/(pk'*rk) - ((Bk*pk)*(Bk*pk)')/(pk'*Bk*pk);
        
        % If the trust region is not accepted
        else
            % Update the trust region
            trustRegion = gamma * vecnorm(pk,'Inf');
        end

        % PART 4: Finish iteration by checking for convergence
        fprintf("Finished iteration: %d", k);
        k = k + 1;
        
        if norm(pk,2) < epsilon
            converged = true;
            x_final = xk;
            solverInformation = struct();
        end

   end

end
