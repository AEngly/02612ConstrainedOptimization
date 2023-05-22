function [c,ceq] = confunHimmelblau(x)

    c = zeros(2,1);
    ceq = zeros(0,1);

    % Bounds
    lower = 0;
    upper = 47;

    % Inequality constraints c(x) <= 0
    tmp = x(1)+2;
    c(1,1) = -(tmp*tmp - x(2)) + lower;
    c(2,1) = (tmp*tmp - x(2)) - upper;
    
end