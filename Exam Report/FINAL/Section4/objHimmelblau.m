function [f, df, d2f] = objHimmelblau(x)

    % Fetch elements from column vector
    x1 = x(1,1);
    x2 = x(2,1);

    % Compute objective
    f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2;

    %Compute gradient
    df = [4*x1*(x1^2 + x2 - 11) + 2*(x1 + x2^2 - 7); 2*(x1^2 + x2 - 11) + 4*x2*(x1 + x2^2 - 7)];

    %Compute Hessian
    d2f = zeros(2,2);
    d2f(1,1) = 4*(x1^2 + x2 - 11);
    d2f(1,2) = 4*(x1 + x2);
    d2f(2,1) = 4*(x1 + x2);
    d2f(2,2) = 4*(x1 + x2^2 - 7) + 8*x2^2 + 2;

end