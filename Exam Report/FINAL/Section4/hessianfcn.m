function Hout = hessianfcn(x, lambda, v)

    % Hessian of objective
    H = [12*x(1)^2 + 4*x(2) - 42, 4*x(1) + 4*x(2);
         4*x(1) + 4*x(2), 12*x(2)^2 + 4*x(1) - 26];
    % Hessian of nonlinear inequality constraint
    gradC1 = [-2 0; 0 0];
    gradC2 = [0 0; 0 0];

    if isempty(v)

        Hout = H + lambda.ineqnonlin(1)*gradC1 + lambda.ineqnonlin(2)*gradC2;

    else

        Hout = H + lambda.ineqnonlin(1)*gradC1 + lambda.ineqnonlin(2)*gradC2;
        Hout = Hout*v;

    end

end