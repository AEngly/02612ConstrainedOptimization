function [c,dc,d2c] = conIqHimmelblau(x)

    x1 = x(1,1);
    x2 = x(2,1);
    c = zeros(2,1);

    % Compute constraints
    c(1) = (x1 + 2)^2 + x2;
    c(2) = -4*x1 + 10*x2;

    % Compute gradients
    dc1 = [2*(x1 + 2); -1];
    dc2 = [-4; 10];
    dc = [dc1 dc2];

    % Compute Hessian of the individual constraints
    d2c1 = [2 0; 0 0];
    d2c2 = [0 0; 0 0];
    d2c = [d2c1; d2c2];
    
end