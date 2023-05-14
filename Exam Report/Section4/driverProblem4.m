%%
% ---------------- DESCRIPTION --------------
%
% Name: driverProblem4.m
% Type: Driver for all testing and plots from problem 4 (NLP)
%
% Assumptions: 
% 
% 1) Functions are sufficiently smooth and [nabla(h) nabla(g)] has full
% column rank.
%
% Problem structure:
%           min     f(x)
%            x
%           s.t.    h(x) = 0
%                   cl  <= g(x) <= cu
%                   l  <= x <= u
%
% Created: 12.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------

%% Global setup

% This contains e.g. options for quadprog or other solvers

%% 4.4) Plot of Himmelblau's Test Problem

% This section tests whether the test problems are generated correctly.
% The function for the objective and constraints are found in
% files named: conHimmelblau.m and objHimmelblau.m

% SETTINGS FOR LABELS, AXIS' AND FILL

upper_colorbar = 200;
lower_colorbar = 0;
granularity_colorbar = 10;

% BOUNDS FOR HIMMELBLAU

c1_l = 0;
c1_u = 47;

c2_l = 0;
c2_u = 70;

x1_l = -5;
x1_u = 5;

x2_l = -5;
x2_u = 5;

% OBJECTIVE VALUES ON GRID

x1 = x1_l:0.05:x1_u;
x2 = x2_l:0.05:x2_u;
[X1, X2] = meshgrid(x1,x2);
F = objfunHimmelblau(X1, X2);

v = lower_colorbar:granularity_colorbar:upper_colorbar;
contour(X1,X2,F,v,"linewidth",2);
colorbar;

% CONSTRAINT BOUNDARIES

yc11 = (x1 + 2).^2 - c1_l; % >= x2
yc12 = (x1 + 2).^2 - c1_u; % <= x2 - c1_u
yc21 = (4 .* x1 + c2_l)./10; % <= x2
yc22 = (4 .* x1 + c2_u)./10; % >= x2

% CONSTRAINT COLORS AND TRANSPARANCY

% ORANGE: [0.8500 0.3250 0.0980]
% BLUE: [0.6350 0.0780 0.1840]

yc1_color = [0 0 0];
yc1_density_l = 0.7; 
yc1_density_u = 0.7; 

yc2_color = [0 0 0];
yc2_density_l = 0.7;
yc2_density_u = 0.7;

% MAKE PLOT

hold on

    % Constraint 1
    h1 = fill([x1_l x1],[x2_u yc11], yc1_color, "facealpha",yc1_density_l);
    h2 = fill([x1_l x1 x1_u],[x2_l yc12 x2_l], yc1_color, "facealpha",yc1_density_u);

    % Constraint 2
    h3 = fill([x1_l x1 x1_u],[x2_l yc21 x2_l], yc2_color, "facealpha",yc2_density_l);
    h4 = fill([x1_l x1 x1_u],[x2_u yc22 x2_u], yc2_color, "facealpha",yc2_density_u);

    % Points
    h5 = plot(-3.5485, -1.4194,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h6 = plot(-0.2983,  2.8956,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h7 = plot(-3.6546,  2.7377,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h8 = plot(3.216440661, 1.286576264,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');
    h9 = plot(3,2,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h10 = plot(-1.4242,0.3315,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');
    h11 = plot(-3.0730,-0.0814,'black', 'MarkerSize', 16, 'Marker', 'diamond', 'MarkerFaceColor', '#D95319');
    h12 = plot(0.0867, 2.8843,'black', 'MarkerSize', 16, 'Marker', 'diamond', 'MarkerFaceColor', '#D95319');
    h13 = plot(-0.4870, -0.1948,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');

hold off

legend([h5, h11, h13],{'Local Minimum', 'Saddle Point', 'Local Maximum'})

xlim([x1_l x1_u])
ylim([x2_l x2_u])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')


saveas(gcf,'./Plots/441_ContourHB.png')

%% 4.5) Scalable test problems (TO DO)

% This section contains tests for adjusted general QP solvers.
% The adjusted implementation are found is separate files named:
% XXXXXXXXXXXXXXX

%% 4.6) Testing of Library Functions (fmincon vs CasADi)

% This section contains tests for primal-dual interior point algorithm for LP.
% The implementation is found in the file named: XXXX

% The test should be on some general LP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics of:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

% Starting point
x0 = [-3; 0];

% Specification of constraints
c1_l = 0;
c1_u = 47;

c2_l = 0;
c2_u = 70;

x1_l = -5;
x1_u = 5;

x2_l = -5;
x2_u = 5;

% SOLUTION: fmincon

% Structure of constraints
A = [-4 10; 4 -10];
b = [c2_u; -c2_l];
Aeq = [];
beq = [];
lb = [x1_l; x2_l];
ub = [x1_u; x2_u];

% Solve with fmincon
optimoptions('fmincon','Display', 'iter', 'SpecifyObjectiveGradient', true, 'TolFun', 1e-8, 'Algorithm', 'interior-point');
[x_fmincon,fval,exitflag,output_fmincon,lambda,grad,hessian] = fmincon(@objfminconHimmelblau, x0, A, b, Aeq, beq, lb, ub, @confunHimmelblau);

% SOLUTION: CasADi with IPOPT

% Uncomment on MacOS with M1 processor
addpath('./../OptimizationSoftware/casadi-3.6.2-osx64-matlab2018b');

% OBS: The solver might be blocked by your operating system. Go to security
% settings to fix this.

% Import Casadi
import casadi.*

x = MX.sym('x', 2);

% Define objective function
t1 = x(1)*x(1)+x(2)-11;
t2 = x(1)+x(2)*x(2)-7;
f = t1*t1 + t2*t2;

% Constraints
c1 = (x(1)+2)^2 - x(2);
c2 = -4*x(1) + 10*x(2);
g = [c1; c2];

nlp = struct;
nlp.x = x;
nlp.f = f;
nlp.g = g;

% Specify options
opts = struct("ipopt", struct('max_iter', 100, 'print_level', 5));

% Create solver instance
S = nlpsol('S', 'ipopt', nlp, opts);

sol_ipopt = S('x0',x0,'ubg',[c1_u; c2_u], ...
    'lbg',[c1_l; c2_l],'lbx',[x1_l;x2_l],'ubx',[x1_u;x2_u]);

disp('----------- RESULTS FOR HIMMELBLAU WITH FMINCON AND IPOPT THROUGH CASADI -----------')
fprintf('\n----- FMINCON ------\n')
fprintf('Iterations: %d\n\n', output_fmincon.iterations);
disp([x_fmincon(1) x_fmincon(2)]);
fprintf('\n----- IPOPT ------\n')
fprintf('See iterations in output from algorithm.\n');
disp(sol_ipopt.x);

%% 4.7) SQP Procedure with Line Search, Damped BFGS and Infeasibility Handling

% This section contains tests for a SQP with Line Search, BFGS
% approximation and infeasibility handling.
% The implementation is found in the file named: SQPLineSearchBFGS.m

% The plots should show the paths on Himmelblau's Test Problem.
%
% We should also provide statistics from scalable problems:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

% Test points (starting points)
%x1 = [-1.35;-0.35];
x1 = [-1.4242;0.3315];
x2 = [-0.85;0.2];
x3 = [-1.8;-0.25];
%x3 = [-2;3];
x4 = [-3.25;0.35];

% This code tests the core algorithm
xlb = [-5; -5];
xub = [5; 5];
cub = [70; 70];
clb = [0; 0];
nonlcon = @(x) conHimmelblau(x);
fun = @(x) objHimmelblau(x);

% Set options
options = struct();
options.maxit = 1000;
options.BFGS = true;
options.stepSolver = "quadprog";
options.l1Penalty = 1000;
options.lineSearch = "all";
options.convergenceRequirement = sqrt(1e-16);
options.infeasibilityPenalty = 100;

[primal_final_1, dual_final_1, solverInformation_1] = SQPLineSearchBFGS(fun,x1,xlb,xub,clb,cub,nonlcon,options);
[primal_final_2, dual_final_2, solverInformation_2] = SQPLineSearchBFGS(fun,x2,xlb,xub,clb,cub,nonlcon,options);
[primal_final_3, dual_final_3, solverInformation_3] = SQPLineSearchBFGS(fun,x3,xlb,xub,clb,cub,nonlcon,options);
[primal_final_4, dual_final_4, solverInformation_4] = SQPLineSearchBFGS(fun,x4,xlb,xub,clb,cub,nonlcon,options);

% Plotting sequence
% SETTINGS FOR LABELS, AXIS' AND FILL

upper_colorbar = 200;
lower_colorbar = 0;
granularity_colorbar = 10;

% BOUNDS FOR HIMMELBLAU

c1_l = 0;
c1_u = 47;

c2_l = 0;
c2_u = 70;

x1_l = -5;
x1_u = 5;

x2_l = -5;
x2_u = 5;

% OBJECTIVE VALUES ON GRID

x1 = x1_l:0.05:x1_u;
x2 = x2_l:0.05:x2_u;
[X1, X2] = meshgrid(x1,x2);
F = objfunHimmelblau(X1, X2);

v = lower_colorbar:granularity_colorbar:upper_colorbar;
contour(X1,X2,F,v,"linewidth",2);
colorbar;

% CONSTRAINT BOUNDARIES

yc11 = (x1 + 2).^2 - c1_l; % >= x2
yc12 = (x1 + 2).^2 - c1_u; % <= x2 - c1_u
yc21 = (4 .* x1 + c2_l)./10; % <= x2
yc22 = (4 .* x1 + c2_u)./10; % >= x2

% CONSTRAINT COLORS AND TRANSPARANCY

% ORANGE: [0.8500 0.3250 0.0980]
% BLUE: [0.6350 0.0780 0.1840]

yc1_color = [0 0 0];
yc1_density_l = 0.7; 
yc1_density_u = 0.7; 

yc2_color = [0 0 0];
yc2_density_l = 0.7;
yc2_density_u = 0.7;

% MAKE PLOT

hold on

    % Constraint 1
    h1 = fill([x1_l x1],[x2_u yc11], yc1_color, "facealpha",yc1_density_l);
    h2 = fill([x1_l x1 x1_u],[x2_l yc12 x2_l], yc1_color, "facealpha",yc1_density_u);

    % Constraint 2
    h3 = fill([x1_l x1 x1_u],[x2_l yc21 x2_l], yc2_color, "facealpha",yc2_density_l);
    h4 = fill([x1_l x1 x1_u],[x2_u yc22 x2_u], yc2_color, "facealpha",yc2_density_u);

    % Points
    h5 = plot(-3.5485, -1.4194,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h6 = plot(-0.2983,  2.8956,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h7 = plot(-3.6546,  2.7377,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h8 = plot(3.216440661, 1.286576264,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');
    h9 = plot(3,2,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h10 = plot(-1.4242,0.3315,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');
    h11 = plot(-3.0730,-0.0814,'black', 'MarkerSize', 16, 'Marker', 'diamond', 'MarkerFaceColor', '#D95319');
    h12 = plot(0.0867, 2.8843,'black', 'MarkerSize', 16, 'Marker', 'diamond', 'MarkerFaceColor', '#D95319');
    h13 = plot(-0.4870, -0.1948,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');

    % Plot solution point
    h14 = plot(solverInformation_1.primalSequence(1,end), solverInformation_1.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#FF0000");
    h15 = plot(solverInformation_1.primalSequence(1,1:end-1), solverInformation_1.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);

    % Plot solution point
    h16 = plot(solverInformation_2.primalSequence(1,end), solverInformation_2.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#0072BD");
    h17 = plot(solverInformation_2.primalSequence(1,1:end-1), solverInformation_2.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);
    
    % Plot solution point
    h18 = plot(solverInformation_3.primalSequence(1,end), solverInformation_3.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#D95319");
    h19 = plot(solverInformation_3.primalSequence(1,1:end-1), solverInformation_3.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);
    
    % Plot solution point
    h20 = plot(solverInformation_4.primalSequence(1,end), solverInformation_4.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#77AC30");
    h21 = plot(solverInformation_4.primalSequence(1,1:end-1), solverInformation_4.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);

hold off

legend([h5, h11, h13, h15, h14, h16, h18, h20],{'Local Minimum', 'Saddle Point', 'Local Maximum', 'Paths', sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x1(1), x1(2)), sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x2(1), x2(2)), sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x3(1),x3(2)), sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x4(1),x4(2))}, 'Location','southeast')

xlim([x1_l x1_u])
ylim([x2_l x2_u])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')

saveas(gcf,'./Plots/SQPLineSearchBFGSInfeasibilityHimmelblau.png')

%% 4.8) Testing of SQP with Trust Region

% This section contains tests for a SQP with Line Search and BFGS
% approximation.
% The implementation is found in the file named: SQPSimpleDamedBFGS.m

% The plots should show the paths on Himmelblau's Test Problem.
%
% We should also provide statistics from scalable problems:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

% Test points (starting points)
%x1 = [-1.35;-0.35];
x1 = [-1.4242;0.3315];
x2 = [0.0;0.90];
x3 = [-1.8;-0.25];
%x3 = [-2;3];
x4 = [-4.25;0.00];

% This code tests the core algorithm
xlb = [-5; -5];
xub = [5; 5];
cub = [70; 70];
clb = [0; 0];
nonlcon = @(x) conHimmelblau(x);
fun = @(x) objHimmelblau(x);

% Set options
options = struct();
options.maxit = 1000;
options.BFGS = true;
options.stepSolver = "quadprog";
options.l1Penalty = 0.1;
options.lineSearch = "all";
options.trustRegion = 0.1;
options.adaptiveTrustRegion = true;
options.infeasibilityPenalty = 50;
options.acceptanceMargin = 0.001;
options.convergenceRequirement = sqrt(1e-16);

[primal_final_1, dual_final_1, solverInformation_1] = SQPTrustRegion(fun,x1,xlb,xub,clb,cub,nonlcon,options);
[primal_final_2, dual_final_2, solverInformation_2] = SQPTrustRegion(fun,x2,xlb,xub,clb,cub,nonlcon,options);
[primal_final_3, dual_final_3, solverInformation_3] = SQPTrustRegion(fun,x3,xlb,xub,clb,cub,nonlcon,options);
[primal_final_4, dual_final_4, solverInformation_4] = SQPTrustRegion(fun,x4,xlb,xub,clb,cub,nonlcon,options);

% Plotting sequence
% SETTINGS FOR LABELS, AXIS' AND FILL

upper_colorbar = 200;
lower_colorbar = 0;
granularity_colorbar = 10;

% BOUNDS FOR HIMMELBLAU

c1_l = 0;
c1_u = 47;

c2_l = 0;
c2_u = 70;

x1_l = -5;
x1_u = 5;

x2_l = -5;
x2_u = 5;

% OBJECTIVE VALUES ON GRID

x1 = x1_l:0.05:x1_u;
x2 = x2_l:0.05:x2_u;
[X1, X2] = meshgrid(x1,x2);
F = objfunHimmelblau(X1, X2);

v = lower_colorbar:granularity_colorbar:upper_colorbar;
contour(X1,X2,F,v,"linewidth",2);
colorbar;

% CONSTRAINT BOUNDARIES

yc11 = (x1 + 2).^2 - c1_l; % >= x2
yc12 = (x1 + 2).^2 - c1_u; % <= x2 - c1_u
yc21 = (4 .* x1 + c2_l)./10; % <= x2
yc22 = (4 .* x1 + c2_u)./10; % >= x2

% CONSTRAINT COLORS AND TRANSPARANCY

% ORANGE: [0.8500 0.3250 0.0980]
% BLUE: [0.6350 0.0780 0.1840]

yc1_color = [0 0 0];
yc1_density_l = 0.7; 
yc1_density_u = 0.7; 

yc2_color = [0 0 0];
yc2_density_l = 0.7;
yc2_density_u = 0.7;

% MAKE PLOT

hold on

    % Constraint 1
    h1 = fill([x1_l x1],[x2_u yc11], yc1_color, "facealpha",yc1_density_l);
    h2 = fill([x1_l x1 x1_u],[x2_l yc12 x2_l], yc1_color, "facealpha",yc1_density_u);

    % Constraint 2
    h3 = fill([x1_l x1 x1_u],[x2_l yc21 x2_l], yc2_color, "facealpha",yc2_density_l);
    h4 = fill([x1_l x1 x1_u],[x2_u yc22 x2_u], yc2_color, "facealpha",yc2_density_u);

    % Points
    h5 = plot(-3.5485, -1.4194,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h6 = plot(-0.2983,  2.8956,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h7 = plot(-3.6546,  2.7377,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h8 = plot(3.216440661, 1.286576264,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');
    h9 = plot(3,2,'black', 'MarkerSize', 16, 'Marker', 'v', 'MarkerFaceColor', '#EDB120');
    h10 = plot(-1.4242,0.3315,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');
    h11 = plot(-3.0730,-0.0814,'black', 'MarkerSize', 16, 'Marker', 'diamond', 'MarkerFaceColor', '#D95319');
    h12 = plot(0.0867, 2.8843,'black', 'MarkerSize', 16, 'Marker', 'diamond', 'MarkerFaceColor', '#D95319');
    h13 = plot(-0.4870, -0.1948,'black', 'MarkerSize', 16, 'Marker', '^', 'MarkerFaceColor', '#A2142F');

    % Plot solution point
    h14 = plot(solverInformation_1.primalSequence(1,end), solverInformation_1.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#FF0000");
    h15 = plot(solverInformation_1.primalSequence(1,1:end-1), solverInformation_1.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);

    % Plot solution point
    h16 = plot(solverInformation_2.primalSequence(1,end), solverInformation_2.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#0072BD");
    h17 = plot(solverInformation_2.primalSequence(1,1:end-1), solverInformation_2.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);
    
    % Plot solution point
    h18 = plot(solverInformation_3.primalSequence(1,end), solverInformation_3.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#D95319");
    h19 = plot(solverInformation_3.primalSequence(1,1:end-1), solverInformation_3.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);
    
    % Plot solution point
    h20 = plot(solverInformation_4.primalSequence(1,end), solverInformation_4.primalSequence(2,end), 'black', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', "#77AC30");
    h21 = plot(solverInformation_4.primalSequence(1,1:end-1), solverInformation_4.primalSequence(2,1:end-1), 'black', 'LineStyle', "--", 'LineWidth', 2);

hold off

legend([h5, h11, h13, h15, h14, h16, h18, h20],{'Local Minimum', 'Saddle Point', 'Local Maximum', 'Paths', sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x1(1), x1(2)), sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x2(1), x2(2)), sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x3(1),x3(2)), sprintf('Solution with x_{0} = (%0.2f, %0.2f)',x4(1),x4(2))}, 'Location','southeast')

xlim([x1_l x1_u])
ylim([x2_l x2_u])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')

saveas(gcf,'./Plots/SQPTrustRegionHimmelblau.png')