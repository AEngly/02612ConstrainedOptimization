% The output of this functions follows the same format
% as the nonlinear constraint function used in fmincon.

function [c,ceq,GC,GCeq] = conHimmelblau(x)

    % 1) Inequalities
    x1 = x(1,1);
    x2 = x(2,1);
    c = zeros(2,1);

    % Evaluate constraints
    c(1) = (x1 + 2)^2 + x2;
    c(2) = -4*x1 + 10*x2;

    
    
    % Compute gradients
    GC1 = [2*(x1 + 2); -1];
    GC2 = [-4; 10];
    GC = [GC1 GC2];

    % 2) Equalities
    ceq = [];
    GCeq = [];

end