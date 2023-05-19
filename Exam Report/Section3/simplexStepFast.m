function [x,fval,exitflag,output,lambda,Abar,bbar,gbar] = simplexStepFast(gbar,Abar,bbar,information,options)

% ---------------- DESCRIPTION --------------
%
% Name: Implementation of Simplex
% Type: Solves linear programs
%
% Problem structure:
%           min     g' x
%            x
%           s.t.    A'x + b = 0
%                   dl <= C'x <= du
%                   l <= x <= u
%
% Syntax [x,fval,exitflag,output,mu,lambda,Abar,bbar,gbar] = simplexCore(g,Aeq,beq,A,cub,clb,lb,ub,options)
%
% Output:
%
% Outputs augmented system (Abar, bbar, and gbar):
%
% min   ggbar'x
%  x
% s.t.  Abar*x = bbar
%
% The system can then be solved with linprog(gbar,[],[],Abar,bbar,zeros(1,16)',[]);
%
% exitflag:
% 1 if optimal solution found
% -1 if problem is unbounded
%
% Created: 15.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------
    
    % Test if the print should be verbose (and get maximum iterations)
    verbose = options.verbose;
    maxiter = options.maxIterations;
    negativeBounds = information.doubleX;

    % Start simplex
    if verbose > 0
        fprintf("\n\n-------- Starting Simplex --------\n\n");
    end

    n = information.n;
    m = information.m;
    meq = information.meq;
    miq = information.miq;
    nOriginal = information.nOriginal;

    if verbose > 1
        fprintf("\nCombined system with slack and surplus: \n\n");
        disp([Abar repelem("=",m)' bbar]);
    end

    % REQUIRE:

    % Set initial basis
    basic = options.initialBasis;
    
    % Compute nonBasic
    I = 1:n; % -> includes surplus +/- and slack variables
    nonBasic = setdiff(I,basic);

    % Computing x's
    B = Abar(:,basic);
    N = Abar(:,nonBasic);
    x = zeros(n,1);
    x(basic) = B^(-1)*bbar;

    % Convergence
    converged = false;
    iter = 0;

    % Set exitflag as not converged
    exitflag = 0;

    if verbose > 0
        fprintf("\nDescription of iterations:\n");
    end

    while ~converged && iter < maxiter

        % Starting first iteration
        iter = iter + 1;

        % Pick basic and nonBasic variables
        xB = x(basic);
        xN = x(nonBasic);

        % Solve for mu (dual for equality constraints) and lambda (dual for
        % nonnegativity constraints)
        gB = gbar(basic);
        mu = B'\gB;
        gN = gbar(nonBasic);
        lambdaN = gN - N'*mu; %("Pricing")

        % Printing
        if verbose > 1
            fprintf("xB:\n");
            disp(xB);            
            fprintf("gB:\n");
            disp(gB);                    
            fprintf("mu:\n");
            disp(mu);
            fprintf("gN:\n");
            disp(gN);
            fprintf("lambdaN:\n");
            disp(lambdaN);
        end

        % Check if converged
        if all(lambdaN >= 0)
            converged = true;
            exitflag = 1;
        else
            % Find most negative entry (and corresponding index)
            [~,s] = min(lambdaN);
            iS = nonBasic(s);
            
            % Search for minimum ratio
            h = B\Abar(:,iS); % Called d in example 13.1 on page 371
            if verbose > 1
                fprintf("h:\n");
                disp(h);
            end
            minAux = Inf;
            J = NaN;
            for i = 1:m
                if h(i) > 0 && xB(i)/h(i) < minAux
                    minAux = xB(i)/h(i);
                    J = i;
                end
            end
            
            % Check if problem is unbounded
            if isnan(J)
                converged = true;
                exitflag = -3;
            else
                j = J(1);
                alpha = xB(j)/h(j); % Called d in example 13.1 on page 371
                x(basic) = x(basic) - alpha*h;
                iJ = basic(j);
                x(iJ) = 0;
                x(iS) = alpha;
                basic(j) = iS;
                nonBasic(s) = iJ;
                if verbose > 0
                    fprintf("Iteration %d: The variable x[%d] leaves the basis. x[%d] enters the basis.\n", iter, iJ, iS);
                end
                B = Abar(:,basic);
                N = Abar(:,nonBasic);
            end
        end
    end

    % Prepare primal variables to return
    x(basic) = xB;
    x(nonBasic) = xN;

    % Prepare fval
    fval = gbar'*x;

    if negativeBounds
        xPositive = x(1:nOriginal);
        xNegative = x(nOriginal+1:2*nOriginal);
        x = xPositive - xNegative;
    else
        x = x(1:nOriginal); 
    end

    % Prepare dual variables to return
    lambda = struct();

    % 1) Start with lower bound dual variables
    lowerLambda = zeros(n,1);
    lowerLambda(nonBasic) = lambdaN;
    lambda.lower = lowerLambda(1:nOriginal);

    % 2) Then upper bound dual variables (0 with our formulation)
    upperLambda = zeros(n,1);
    lambda.upper = upperLambda(1:nOriginal);

    % 3) Then dual variables on equality constraints (first meq components of mu)
    eqLambda = mu(1:meq);
    if isempty(eqLambda)
        eqLambda = [];
    end
    lambda.eqlin = -eqLambda(1:meq); % Lagrange is defined with opposite sign in linprog

    % 4) Then dual variables on inequality constraints (last components of mu)
    ineqLambda = mu(meq+1:meq+2*miq);
    if isempty(ineqLambda)
        ineqLambda = [];
    end
    lambda.ineqlin = ineqLambda; % Lagrange is defined with opposite sign in linprog

    % Prepare output struct
    output = struct();
    output.iterations = iter;
    output.algorithm = 'revised simplex';
    if exitflag == 1
        output.optimalBasis = basic;
        message = 'Optimal solution found.';
    elseif exitflag == 0
        message = 'Number of iterations exceeded.';
    elseif exitflag == -3
        message = 'Problem is unbounded.';
    end
    output.message = message;
    if verbose > 0
        fprintf("\nExit status from simplex:\n%s\n", message);
    end
end