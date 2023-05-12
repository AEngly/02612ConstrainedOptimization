% 'y' is lagrangian multipliers for the equality constraints
% 'z' are the lagrangian multipliers for the inequality constraints
% This implementation follows the slides from
% "SequentialQuadraticProgramming".

function [alpha_final] = lineSearch(fun,x,xl,xu,cl,cu,con,options,y,z)

    % Starting value
    alpha = 1;
    armijo_constant = 0.1;

    % Powell's update
    lambda = max(abs(y));
    mu = max(abs(z));

    % Compute merit values
    [f_0,df_0] = fun(x);
    [f_alpha,~] = fun(x + alpha*p);
    [c_0,ceq_0] = con(x);
    [c_alpha,ceq_alpha] = con(x + alpha*p);
    phi_0 = phi(f_0,ceq_0,c_0,mu,lambda);
    phi_alpha = phi(f_alpha,ceq_alpha,c_alpha,mu,lambda);
    directional_derivative_phi_0 = phi(transpose(df_0)*p,ceq_0,c_0,-lambda,-mu);

    % Check Armijo condition
    while phi_alpha > phi_0 + armijo_constant*alpha*directional_derivative_phi_0
        a = phi_alpha + (phi_0 + directional_derivative_phi_0*alpha);
        alpha_min = (-directional_derivative_phi_0/(2*a));
        alpha = min(0.9*alpha, max(alpha_min, 0.1*alpha));
    end

    alpha_final = alpha;

end

% This penality function is defined in (15.24)
function val = phi(f,ceq,c,lambda,mu)
    val = f + lambda'*abs(ceq) + mu'*abs(min(0,c));
end