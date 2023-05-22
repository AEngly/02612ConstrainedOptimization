function [x,fval,exitflag,output,lambdaOutput] = LP_InteriorPointPrimalDual(g,A,b,C,dub,dlb,lb,ub,options)

% ---------------- DESCRIPTION --------------
%
% Name: LP_InteriorPointPrimalDual.m
% Type: Implementation of Primal-Dual Interior Point method for LP
%
% Assumptions: 
% 
% 1) Equality constraint matrix A has full column rank.
%
% Problem structure:
%           min     g'x
%            x
%           s.t.    A'x + b = 0
%                   dl  <= C'x <= du
%                   l  <= x <= u
%
% Created: 12.05.2023
% INSPIRATION: Implementation by John Bagterp (from course 02612 Constrained Optimization)
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------

%%
    % Extract relevant constants
    [mE,n1]=size(A);
    [mI,n2] = size(C);
    n = max(n1,n2);
    mC = 2*mI+2*n;  

    % Flip coefficient matrices
    A = A';
    C = C';
    
    % Set algorithm specific constants
    maxit = 100;
    tolL = 1.0e-9;
    tolA = 1.0e-9;
    tolC = 1.0e-9;
    tols = 1.0e-9;
    
    eps = 0.99;
    
    x = zeros(n,1);
    lambda = ones(mE,1);
    mu = ones(2*mI + 2*n,1);
    s = ones(2*mI + 2*n,1);

    Cbar = [-C C -eye(n) eye(n)];
    d = [-dub;dlb;-ub;lb];

    % ====================================================================
    % Make sure none of the matrix operations will fail
    % ====================================================================
    
    while(any(s==0))
        x = x+0.000001;
        s = Cbar'*x - d;
    end


    % ====================================================================
    % Compute starting point
    % ====================================================================
    
    % Compute residuals
    rL = g - A*lambda - Cbar*mu;              % Lagrangian gradient
    rA = -A'*x - b;                           % Equality Constraint
    rC = -Cbar'*x + d + s;                    % Inequality Constraint
    rS = mu.*s;                               % Complementarity
    eta = sum(rS)/mC;                         % Duality gap

    % Form and Factorize Augmented System
    mudivs = mu./s;
    H = Cbar*diag(mudivs)*Cbar';
    augmentKKT = [H -A; -A' zeros(mE,mE)];
    [L,D,p] = ldl(augmentKKT,'lower','vector');
    
    % Affine Step
    temp = diag(mu./s)*(rC - rS./mu);
    rhs1 = rL - Cbar*temp;
    rhs2 = rA;
    rhs = -[rhs1;rhs2];

    dxdlAff = zeros(n+mE,1);
    dxdlAff(p) = L'\(D\(L\rhs(p))); % Should be solved with LDL
    dxAff = dxdlAff(1:n);
    dlAff = dxdlAff(n+1:end);
    dmuAff = -diag(mu./s)*Cbar'*dxAff + temp;
    dsAff = -rS./mu - (s./mu).*dmuAff;
    
    % Find starting point
    x = x;
    lambda = lambda;
    mu = max(1,abs(mu + dmuAff));
    s = max(1,abs(s + dsAff));

    % Recompute residuals
    rL = g - A*lambda - Cbar*mu;           % Lagrangian gradient
    rA = -A'*x - b;                        % Equality Constraint
    rC = -Cbar'*x + d + s;                 % Inequality Constraint
    rS = mu.*s;                            % Complementarity
    eta = sum(rS)/mC;                      % Duality gap

    % Converged
    Converged = (norm(rL,inf) <= tolL) && ...
                (norm(rA,inf) <= tolA) && ...
                (norm(rC,inf) <= tolC) && ...
                (abs(eta) <= tols);

    %%        
    iter = 0;
    while ~Converged && (iter<maxit)

        iter = iter+1;
        
        % ====================================================================
        % Form and Factorize Augmented System
        % ====================================================================
        
        mudivs = mu./s;
        H = Cbar*diag(mudivs)*Cbar';
        augmentKKT = [H -A; -A' zeros(mE,mE)];
        [L,D,p] = ldl(augmentKKT,'lower','vector');
        
        % ====================================================================
        % Affine Step
        % ====================================================================

        temp = diag(mu./s)*(rC - rS./mu);
        rhs1 = rL - Cbar*temp;
        rhs2 = rA;
        rhs = -[rhs1;rhs2];
        
        dxdlAff = zeros(n+mE,1);
        dxdlAff(p) = L'\(D\(L\rhs(p))); % Should be solved with LDL
        %dxdlAff = augmentKKT\rhs; % Should be solved with LDL
        dxAff = dxdlAff(1:n);
        dlAff = dxdlAff(n+1:end);
        dmuAff = -diag(mu./s)*Cbar'*dxAff + temp;
        dsAff = -rS./mu - (s./mu).*dmuAff;
        
        % Step length
        idx = find(dmuAff < 0.0);
        alpha = min([1.0; -mu(idx,1)./dmuAff(idx,1)]);
        
        idx = find(dsAff < 0.0);
        beta = min([1.0; -s(idx,1)./dsAff(idx,1)]);
        
        % ====================================================================
        % Center Parameter
        % ====================================================================
        
        muAff = mu + alpha*dmuAff;
        sAff = s + beta*dsAff;
        etaAff = sum(muAff.*sAff)/mC;
        sigma = (etaAff/eta)^3;
        tau = sigma*eta;    
    
        % ====================================================================
        % Center-Corrector Step
        % ====================================================================
        
        rSbar = rS + dsAff.*dmuAff - tau;
        temp = diag(mu./s)*(rC - rSbar./mu);
        rLbar = rL - Cbar*temp;
        rhs = -[rLbar;rA];

        % Solve augmented system
        dxdl = zeros(n+mE,1);
        dxdl(p) = L'\(D\(L\rhs(p))); % Should be solved with LDL
        %dxdl = augmentKKT\rhs; % Should be solved with LDL
        dx = dxdl(1:n);
        dl = dxdl(n+1:end);

        % Compute the rest of the results
        dmu = -diag(mu./s)*Cbar'*dx + temp;
        ds = -rSbar./mu - (s./mu).*dmu;
        
        % Step length
        idx = find(dmu < 0.0);
        alpha = min([1.0; -mu(idx,1)./dmu(idx,1)]);
        
        idx = find(ds < 0.0);
        beta = min([1.0; -s(idx,1)./ds(idx,1)]);

        stepFactor = min(alpha,beta);
    
        % ====================================================================
        % Take step 
        % ====================================================================

        x = x + (eps*stepFactor)*dx;
        lambda = lambda + (eps*stepFactor)*dl;
        mu = mu + (eps*stepFactor)*dmu;
        s = s + (eps*stepFactor)*ds;
        
        % ====================================================================
        % Residuals and Convergence
        % ====================================================================

        rL = g - A*lambda - Cbar*mu;              % Lagrangian gradient
        rA = -A'*x - b;                                 % Equality Constraint
        rC = -Cbar'*x + d + s;  % Inequality Constraint
        rS = mu.*s;                                   % Complementarity
        eta = sum(rS)/mC;                             % Duality gap
        
        % Converged
        Converged = (norm(rL,inf) <= tolL) && ...
                    (norm(rA,inf) <= tolA) && ...
                    (norm(rC,inf) <= tolC) && ...
                    (abs(eta) <= tols);

    end
    
    %%

    % Compute function value
    fval = g'*x;

    % Prepare dual variables for output
    lambdaOutput = struct();
    output = struct();

    % 1) Start with lower bound dual variables
    lambdaOutput.lower = mu(2*mI+n+1:end);
    lambdaOutput.upper = mu(2*mI+1:2*mI+n);
    lambdaOutput.eqlin = -lambda;
    lambdaOutput.ineqlin = mu(1:2*mI); 
   
    if Converged
        exitflag = 1;
        output.iterations = iter;
    else
        exitflag = -1;
        output.iterations = NaN;
        lambdaOutput.lower = [];
        lambdaOutput.upper = [];
        lambdaOutput.eqlin = []; 
        lambdaOutput.ineqlin = [];
    end

end
    