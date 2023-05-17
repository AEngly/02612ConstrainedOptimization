function [x,fval,exitflag,output,lambda,Abar,bbar,gbar] = simplex(g,Aeq,beq,A,cub,clb,lb,ub,options)

    % PHASE 1: Find a basic feasible solution

    % 1) Construct system as usually

    verbose = options.verbose;
    maxiter = options.maxIterations;
    negativeBounds = false;

    % Start simplex
    if verbose > 0
        fprintf("\n\n-------- Starting Simplex --------\n\n");
    end

    % Compute relevant lengths
    if ~isempty(Aeq)
        meq = size(Aeq,1);
        n = size(Aeq,2);
    else
        meq = 0;
    end
    if ~isempty(A)
        miq = size(A,1);
        n = size(A,2);
    else
        miq = 0;
    end

    slack1 = size(cub,1);
    slack2 = size(clb,1);
    slack3 = size(ub,1);
    slack4 = size(lb,1);
    slackAll = slack1 + slack2 + slack3 + slack4;

    if any(lb < 0) || any(ub < 0)
        negativeBounds = true;
        if meq > 0
            row1 = [Aeq -Aeq zeros(meq,slack1) zeros(meq,slack2) zeros(meq,slack3) zeros(meq,slack4)];
            if verbose > 0
                fprintf("Adding row for equality constraints: (%d,%d)\n", size(Aeq,1),2*n + slackAll);
            end
        else
            row1 = [];
        end
        if miq > 0 && ~isempty(cub)
            row2 = [-A A -eye(slack1) zeros(miq,slack2) zeros(miq,slack3) zeros(miq,slack4)];
            if verbose > 0
                fprintf("Adding row for upper bound on inequality constraints: (%d,%d)\n", size(A,1),2*n + slackAll);
            end
        else
            row2 = [];
        end
        if miq > 0 && ~isempty(clb)
            row3 = [A -A zeros(miq,slack1) -eye(slack2) zeros(miq,slack3) zeros(miq,slack4)];
            if verbose > 0
                fprintf("Adding row for lower bound on inequality constraints: (%d,%d)\n", size(A,1),2*n + slackAll);
            end
        else
            row3 = [];
        end
        if ~isempty(ub)
            row4 = [-eye(n) eye(n) zeros(n,slack1) zeros(n,slack2) -eye(n) zeros(n,slack4)];
            if verbose > 0
                fprintf("Adding row for upper bound on variables: (%d,%d)\n", size(ub,1),2*n + slackAll);
            end
        else
            row4 = [];
        end
        if ~isempty(lb)
            row5 = [eye(n) -eye(n) zeros(n,slack1) zeros(n,slack2) zeros(n,slack3) -eye(n)];
            if verbose > 0
                fprintf("Adding row for lower bound on variables: (%d,%d)\n", size(lb,1),2*n + slackAll);
            end
        else
            row5 = [];
        end
        gbar = [g' -g' zeros(1,slackAll)]';
    else
        if meq > 0
            row1 = [Aeq zeros(meq,slack1) zeros(meq,slack2) zeros(meq,slack3) zeros(meq,slack4)];
            if verbose > 0
                fprintf("Adding row for equality constraints: (%d,%d)\n", size(Aeq,1),2*n + slackAll);
            end
        else
            row1 = [];
        end
        if miq > 0 && ~isempty(cub)
            row2 = [-A -eye(slack1) zeros(miq,slack2) zeros(miq,slack3) zeros(miq,slack4)];
            if verbose > 0
                fprintf("Adding row for upper bound on inequality constraints: (%d,%d)\n", size(A,1),2*n + slackAll);
            end
        else
            row2 = [];
        end
        if miq > 0 && ~isempty(clb)
            row3 = [A zeros(miq,slack1) -eye(slack2) zeros(miq,slack3) zeros(miq,slack4)];
            if verbose > 0
                fprintf("Adding row for lower bound on inequality constraints: (%d,%d)\n", size(A,1),2*n + slackAll);
            end
        else
            row3 = [];
        end
        if ~isempty(ub)
            row4 = [-eye(n) zeros(n,slack1) zeros(n,slack2) -eye(n) zeros(n,slack4)];
            if verbose > 0
                fprintf("Adding row for upper bound on variables: (%d,%d)\n", size(ub,1),2*n + slackAll);
            end
        else
            row4 = [];
        end
        if ~isempty(lb)
            row5 = [eye(n) zeros(n,slack1) zeros(n,slack2) zeros(n,slack3) -eye(n)];
            if verbose > 0
                fprintf("Adding row for lower bound on variables: (%d,%d)\n", size(lb,1),2*n + slackAll);
            end
        else
            row5 = [];
        end
        gbar = [g' zeros(1,slackAll)]';
    end

    % Save original variables
    nOriginal = n;
    
    % Construct the program (with surplus and artificial variables)
    Abar = [row1; row2; row3; row4; row5];
    bbar = [-beq; -cub; clb; -ub; lb];

    % Find rows and columns
    m = size(Abar,1);
    n = size(Abar,2);

    % Now the system has been brough into Ax = b with x >= 0.
    % Then we can prepare phase 2
    diagE = ones(1,m);
    diagE(bbar < 0) = -1;
    E = diag(diagE);
    A_phaseOne = [Abar E];
    b_phaseOne = bbar;
    g_phaseOne = [zeros(1,n) ones(1,m)]';

    if verbose > 1
        fprintf("\nCombined system with slack and surplus: \n\n");
        disp([Abar repelem("=",m)' bbar]);
    end

    % 2) Find feasible solution
    % Run our implementation
    options = struct();
    options.maxIterations = 100;
    options.verbose = 1;
    options.initialBasis = (n+1):(n+m);
    [x,fval,exitflag,output,lambda,Abar,bbar,gbar] = simplexCore2(g_phaseOne,A_phaseOne,b_phaseOne,[],[],[],[],[],options);
    disp(x);
    % PHASE 2: Find a basic feasible solution

    % Run our implementation
    % The book states that we need to work on the same system as the basis
    % can contain some of those
    optionsPhaseTwo = struct();
    optionsPhaseTwo.maxIterations = 100;
    optionsPhaseTwo.verbose = 1;
    optionsPhaseTwo.initialBasis = output.optimalBasis;
    %[x,fval,exitflag,output,lambda,Abar,bbar,gbar] = simplexCore2(g_phaseOne,A_phaseOne,b_phaseOne,[],[],[],[],[],optionsPhaseTwo);

    
end