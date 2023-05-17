function [Abar, bbar, gbar, information] = standardForm(g,Aeq,beq,A,cub,clb,lb,ub)

    % Add information
    information = struct();
    information.doubleX = false;
    verbose = 0;

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
        information.doubleX = true;
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

    % Construct the program (with surplus and artificial variables)
    Abar = [row1; row2; row3; row4; row5];
    bbar = [-beq; -cub; clb; -ub; lb];
    gbar = gbar;

end