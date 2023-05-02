%% EXERCISE 1.1 (FMINCON)

x0 = [0;0];
xl = [-5;-5];
xu = [5;5];

% First we pretend that there is no linear equality/inequality constraints
A = zeros(0,2);
b = zeros(0,1);
Aeq = zeros(0,2);
beq = zeros(0,1);

% Parameters (in this case we do not use parameters)
p=[];

% Call fmincon
options = optimoptions( 'fmincon',...
    'Display','none',...
    'Algorithm','interior-point');

[x,fval,exitflag,output]=fmincon(...
    @objfunHimmelblau, x0, ...
    A, b, Aeq, beq, ...
    xl, xu, ...
    @confunHimmelblau, ...
    options, ...
    p);

sprintf("x[1] = %f", x(1))
sprintf("x[2] = %f", x(2))
sprintf("f(x) = %f", fval)
disp(output)

%% EXERCISE 1.2

x0 = [0;0];
xl = [-5;-5];
xu = [5;5];

% First we pretend that there is no linear equality/inequality constraints
A = zeros(0,2);
b = zeros(0,1);
Aeq = zeros(0,2);
beq = zeros(0,1);

% Parameters (in this case we do not use parameters)
p=[];

% Call fmincon
options = optimoptions( 'fmincon',...
    'Display','none',...
    'Algorithm','interior-point', ...
    'SpecifyObjectiveGradient',true, ...
    'HessianFcn',@hessianfcn, ...
    'TolProjCG', 0.1, ...
    'TolProjCGAbs', 1.0000e-4);

[x,fval,exitflag,output]=fmincon(...
    @objfunHimmelblau, x0, ...
    A, b, Aeq, beq, ...
    xl, xu, ...
    @confunHimmelblau, ...
    options, ...
    p);

sprintf("x[1] = %f", x(1))
sprintf("x[2] = %f", x(2))
sprintf("f(x) = %f", fval)
disp(output)

%% EXERCISE 2.1 (QUADPROG) 

H = [1 -1; -1 2]; 
f = [-2; -6];
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];

[x,fval,exitflag,output,lambda] = ...
   quadprog(H,f,A,b);

output, x, fval, exitflag, lambda.ineqlin


