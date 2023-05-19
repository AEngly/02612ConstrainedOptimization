function [x,fval,exitflag,output,lambda] = simplex(g,Aeq,beq,A,cub,clb,lb,ub,options)

    % PHASE 1: Find a basic feasible solution

    % 1) Construct system as usually

    [AbarPhase1, bbarPhase1, gbar, information] = standardForm(g,Aeq,beq,A,cub,clb,lb,ub);
    maxiter = options.maxIterations;

    % 2) Run simplexStep on l_inf-regression to find initial basis

    m = information.m;
    n = information.n;
    nOriginal = information.nOriginal;
    diagE = ones(1,m);
    diagE(bbarPhase1 < 0) = -1;
    E = diag(diagE);
    A_phaseOne = [AbarPhase1 E];
    b_phaseOne = bbarPhase1;
    g_phaseOne = [zeros(1,n) ones(1,m)]';

    % 2.2) Find feasible solution
    options.initialBasis = (n+1):(n+m);
    [x,fval,exitflag,output,lambda,Abar,bbar,~] = simplexStep(g_phaseOne,A_phaseOne,b_phaseOne,information,options);
    
    % PHASE 2: Use basis to solve original problem

    if fval == 0

        % 2.3 If non of the artifical variables are in the basis, then we can
        % use the same basis

        information.n = size(AbarPhase1,2);
        information.m = size(AbarPhase1,1);
        g_phaseTwo = gbar;
        options.initialBasis = output.optimalBasis;
        [x,fval,exitflag,output,lambda,Abar,bbar,gbar] = simplexStep(g_phaseTwo,AbarPhase1,b_phaseOne,information,options);
    
    else

        output.exitflag = -2;
        output.message = 'Problem is infeasible.';

    end

    % TERMINATED

end