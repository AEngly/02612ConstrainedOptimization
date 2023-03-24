
%% --------- VERSION 1 -----------

  x0 = [0;0];
  xl = [-5;-5];
  xu = [5;5];

  % -4x1 + 10 x2>= 0 represented as A x<= b
  A = [4 -10];
  b = 0;
  Aeq = zeros(0,2);
  beq = zeros(0,1);
  % initial point % lower bounds % upper bounds

  % Parameters (in this case we do not use parameters)
  p = [];

  % Call fmincon
  options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, ...
                         'SpecifyConstraintGradient', true, 'Display', ...
                         'none', 'Algorithm','interior-point');
 
  [x,fval,exitflag,output]=fmincon(...
  @objfungradHimmelblau, x0, ...
  A, b, Aeq, beq, ...
  xl, xu, ...
  @confungradHimmelblau, ...
  options, ...
  p);

x,fval,output
