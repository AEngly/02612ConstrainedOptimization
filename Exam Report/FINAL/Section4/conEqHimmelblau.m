function [c,dc,d2c] = con2Himmelblau(x)

    x1 = x(1,1);
    x2 = x(2,1);

    % Compute constraints
    c = -4*x1 + 10*x2;
    % Compute gradient
    dc = [-4; 10];
    % Compute Hessian
    d2c = [0 0; 0 0];
    
end