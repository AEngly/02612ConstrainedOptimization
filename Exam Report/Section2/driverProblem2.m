% ---------------- DESCRIPTION --------------
%
% Name: driverProblem2.m
% Type: Driver for all testing and plots from problem 2 (quadratic programming)
%
% Assumptions: 
% 
% 1) Equality constraint matrix A has full column rank.
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%                   dl <= C'x <= du
%                   l <= x <= u
%
% Created: 18.04.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
%
% ---------------- IMPLEMENTATION --------------

%% Global setup

options_as = optimoptions('quadprog','Algorithm', 'active-set', 'MaxIterations', 1000);
options_ip = optimoptions('quadprog','MaxIterations', 1000);
addpath(genpath('../TestTools'));
% This contains e.g. options for quadprog or other solvers

%% 2.1) Construction of random size dependent test problems for QP

% This section tests whether the test problems are generated correctly.
% The function for generating the actual problems are found in separate
% files named:
% GeneratorQPRandom.m (from lecture 6 in course material)
% GeneratorQP
% GeneratorHuberFittingQP.m
% GeneratorOptimalControlQP.m
% GeneratorPortfolioOptimizationQP.m
% GeneratorSupportVectormachineQP.m

%% 2.5a) Primal active set QP solver (initial point linprog)

% This section contains tests for primal active set QP solver.
% The adjusted implementation are found is separate files named:
% QP_primalActiveSet.m
%
% The test should be on some general QP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics of:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

%% Test minimal QP and plot path

maxiter = 100;
[H,g,A,b,C,dl,du,l,u] = MinimalQP();
[xas, iter, x_k] = QP_PrimalActiveSet_linprog(H, g, A, b, C, dl, du, l, u, maxiter,'verbose', 'x_k');
xstar = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);
disp(norm(xas-xstar) <1e-9)

% Then we can plot the points on the contours
x = -2:0.05:3;
y = -2:0.05:3;
[X,Y] = meshgrid(x,y);

F = (X - 1).^2 + (Y - 2.5).^2;

v = -20:2:20;
figure,
[c,h]=contour(X,Y,F,v,"linewidth",2, 'DisplayName','Contour');
colorbar

yc1 = (x + 2)./2; % larger than equal ... x >= 2y - 2
yc2 = (-x + 2)./(2); % y less ...
yc3 = (x - 2)./2; % y greater ...
yc4 = (-x - 1)./(2); % y less ...

hold on
    plot(x,y,'-', 'LineWidth', 5, 'Color', 'magenta')
    fill([x, x(end), x(1)], [yc1, y(end), y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc2, y(end), y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc3, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc4, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    for i = 1:iter-1
        plot(x_k(1,i:(i+1)),x_k(2,i:(i+1)),':', 'LineWidth', (7-2*i))
    end
    plot(x_k(1,iter),x_k(2,iter),'.', 'Color', 'black', 'MarkerSize', 20)
hold off

%legend('','','','','','','Minimizer')

xlim([-1 1])
ylim([-1 1])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')
%legend()
saveas(gcf,'./Plots/2511_ActiveSetLinprogPath.png')


%% Solution quality vs quadprog plot numerical performance

% test setup
number_of_problems = 3; %Change this number as algorithms are added to the loop
n_max = 50;
step = 2;
Variable_names = ["QP","BoxQP","HuberFitting"];

%QP setup
alpha = 0.01;
beta = 4; % setting from OSQP paper
density = 0.15; % 15% must be non-zero

% solver setup
maxiter=1000;

% loop setup
problem_sizes = 10:step:n_max;
l = size(problem_sizes,2);
residual = zeros(l,number_of_problems);
iterations = zeros(l,number_of_problems);
cpu_t = zeros(l,number_of_problems);
j = 1;

for n = problem_sizes
    [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,1)] = QP_PrimalActiveSet_linprog(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,1) = toc;
    residual(j,1) = norm(x-xstar);

    [H,g,A,b,C,dl,du,l,u] = GeneratorBoxQP(n,alpha,beta, density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,2)] = QP_PrimalActiveSet_linprog(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,2) = toc;
    residual(j,2) = norm(x-xstar);
    
    [H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP_bn(n,beta,density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,3)] = QP_PrimalActiveSet_linprog(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,3) = toc;
    residual(j,3) = norm(x-xstar);

    j = j + 1;
end


figure,
hold on

scatter(problem_sizes,residual(:,1), '.');
scatter(problem_sizes,residual(:,2), 'o');
scatter(problem_sizes,residual(:,3), 'x');
%scatter(problem_sizes,residual(:,4), '+');
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("residual")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/2512_Numerical_performance_QP_PrimalActiveSet_linprog.png')

CPU_T = table(problem_sizes',cpu_t(:,1),cpu_t(:,2),cpu_t(:,3),'VariableNames',['n',Variable_names]);
table2latex(CPU_T, './2511cpuTime.tex')

Iter_T = table(problem_sizes',iterations(:,1),iterations(:,2),iterations(:,3),'VariableNames',['n',Variable_names]);
table2latex(Iter_T, './2512Iterations.tex')

%% 2.5b) Primal active set QP solver (initial point primal active set)

% This section contains tests for primal active set QP solver.
% The adjusted implementation are found is separate files named:
% QP_primalActiveSet.m
%
% The test should be on some general QP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics of:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

%% Test minimal QP and plot path

maxiter = 100;
[H,g,A,b,C,dl,du,l,u] = MinimalQP();
[xas, iter, x_k] = QP_PrimalActiveSet_PAS(H, g, A, b, C, dl, du, l, u, maxiter,'verbose', 'x_k');
xstar = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);
disp(norm(xas-xstar) <1e-9)

% Then we can plot the points on the contours
x = -2:0.05:3;
y = -2:0.05:3;
[X,Y] = meshgrid(x,y);

F = (X - 1).^2 + (Y - 2.5).^2;

v = -20:2:20;
figure,
[c,h]=contour(X,Y,F,v,"linewidth",2, 'DisplayName','Contour');
colorbar

yc1 = (x + 2)./2; % larger than equal ... x >= 2y - 2
yc2 = (-x + 2)./(2); % y less ...
yc3 = (x - 2)./2; % y greater ...
yc4 = (-x - 1)./(2); % y less ...

hold on
    plot(x,y,'-', 'LineWidth', 5, 'Color', 'magenta')
    fill([x, x(end), x(1)], [yc1, y(end), y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc2, y(end), y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc3, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc4, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    for i = 1:iter-1
        plot(x_k(1,i:(i+1)),x_k(2,i:(i+1)),':', 'LineWidth', (7-2*i))
    end
    plot(x_k(1,iter),x_k(2,iter),'.', 'Color', 'black', 'MarkerSize', 20)
hold off

%legend('','','','','','','Minimizer')

xlim([-1 1])
ylim([-1 1])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')
%legend()
saveas(gcf,'./Plots/2521_ActiveSetLinprogPath.png')

%% Solution quality vs quadprog plot numerical performance

% test setup
number_of_problems = 3; %Change this number as algorithms are added to the loop
n_max = 50;
step = 2;
Variable_names = ["QP","BoxQP","HuberFitting"];

%QP setup
alpha = 0.01;
beta = 4; % setting from OSQP paper
density = 0.15; % 15% must be non-zero

% solver setup
maxiter=1000;

% loop setup
problem_sizes = 10:step:n_max;
l = size(problem_sizes,2);
residual = zeros(l,number_of_problems);
iterations = zeros(l,number_of_problems);
cpu_t = zeros(l,number_of_problems);
j = 1;

for n = problem_sizes
    [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,1)] = QP_PrimalActiveSet_PAS(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,1) = toc;
    residual(j,1) = norm(x-xstar);

    [H,g,A,b,C,dl,du,l,u] = GeneratorBoxQP(n,alpha,beta, density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,2)] = QP_PrimalActiveSet_PAS(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,2) = toc;
    residual(j,2) = norm(x-xstar);
    
    [H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP_bn(n,beta,density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,3)] = QP_PrimalActiveSet_PAS(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,3) = toc;
    residual(j,3) = norm(x-xstar);

    j = j + 1;
end


figure,
hold on

scatter(problem_sizes,residual(:,1), '.');
scatter(problem_sizes,residual(:,2), 'o');
scatter(problem_sizes,residual(:,3), 'x');
%scatter(problem_sizes,residual(:,4), '+');
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("residual")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/2522_Numerical_performance_QP_PrimalActiveSet_PAS.png')

CPU_T = table(problem_sizes',cpu_t(:,1),cpu_t(:,2),cpu_t(:,3),'VariableNames',['n',Variable_names]);
table2latex(CPU_T, './2521cpuTime.tex')

Iter_T = table(problem_sizes',iterations(:,1),iterations(:,2),iterations(:,3),'VariableNames',['n',Variable_names]);
table2latex(Iter_T, './2522Iterations.tex')
%% 2.7) Testing of primal-dual interior point tailored for general LP

% This section contains tests for primal-dual interior point QP solver.
% The adjusted implementation are found is separate files named:
% QP_primalActiveSet.m
%
% The test should be on some general QP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics of:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

%% Test minimal QP

maxiter = 100;
[H,g,A,b,C,dl,du,l,u] = MinimalQP();
[xas, iter, x_k] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);
xstar = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);
norm(xas-xstar)

% Then we can plot the points on the contours
x = -2:0.05:3;
y = -2:0.05:3;
[X,Y] = meshgrid(x,y);

F = (X - 1).^2 + (Y - 2.5).^2;

v = -20:2:20;
figure,
[c,h]=contour(X,Y,F,v,"linewidth",2, 'DisplayName','Contour');
colorbar

yc1 = (x + 2)./2; % larger than equal ... x >= 2y - 2
yc2 = (-x + 2)./(2); % y less ...
yc3 = (x - 2)./2; % y greater ...
yc4 = (-x - 1)./(2); % y less ...

hold on
    plot(x,y,'-', 'LineWidth', 5, 'Color', 'magenta')
    fill([x, x(end), x(1)], [yc1, y(end), y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc2, y(end), y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc3, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc4, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    for i = 1:iter-1
        plot(x_k(1,i:(i+1)),x_k(2,i:(i+1)),':', 'LineWidth', 5)
    end
    plot(x_k(1,iter),x_k(2,iter),'.', 'Color', 'black', 'MarkerSize', 20)
hold off

%legend('','','','','','','Minimizer')

xlim([-1 1])
ylim([-1 1])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')
%legend()
saveas(gcf,'./Plots/271_InteriorPointPath.png')

%% Solution quality vs quadprog plot numerical performance

% test setup
number_of_problems = 3; %Change this number as algorithms are added to the loop
n_max = 50;
step = 2;
Variable_names = ["QP","BoxQP","HuberFitting"];

%QP setup
alpha = 0.01;
beta = 4; % setting from OSQP paper
density = 0.15; % 15% must be non-zero

% solver setup
maxiter=1000;

% loop setup
problem_sizes = 10:step:n_max;
l = size(problem_sizes,2);
residual = zeros(l,number_of_problems);
iterations = zeros(l,number_of_problems);
cpu_t = zeros(l,number_of_problems);
j = 1;

for n = problem_sizes
    [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,1)] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,1) = toc;
    residual(j,1) = norm(x-xstar);

    [H,g,A,b,C,dl,du,l,u] = GeneratorBoxQP(n,alpha,beta, density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,2)] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,2) = toc;
    residual(j,2) = norm(x-xstar);
    
    [H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP_bn(n,beta,density);
    x = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_as);
    tic;
    [xstar, iterations(j,3)] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,3) = toc;
    residual(j,3) = norm(x-xstar);

    j = j + 1;
end


figure,
hold on

scatter(problem_sizes,residual(:,1), '.');
scatter(problem_sizes,residual(:,2), 'o');
scatter(problem_sizes,residual(:,3), 'x');
%scatter(problem_sizes,residual(:,4), '+');
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("residual")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/271_Numerical_performance_QP_PrimalActiveSet_linprog.png')

CPU_T = table(problem_sizes',cpu_t(:,1),cpu_t(:,2),cpu_t(:,3),'VariableNames',['n',Variable_names]);
table2latex(CPU_T, './271cpuTime.tex')

Iter_T = table(problem_sizes',iterations(:,1),iterations(:,2),iterations(:,3),'VariableNames',['n',Variable_names]);
table2latex(Iter_T, './272Iterations.tex')
%% 2.8) Compare implementations with quadprog (MOSEK, Gurobi, qpOASES, OSQP, YALMIP and cvx)

% The test should be on some general QP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics from scalable problems:
% 1) Size of problem (n) on x-axis, and time-to-completion on y-axis.

%%

% test setup
number_of_problems = 5; %Change this number as algorithms are added to the loop
n_max =300;
step = 10;
Variable_names = ["quadprogAS","quadprogIP","PrimalActiveSet","InteriorPoint"];

%QP setup
alpha = 0.01;
beta = 4; % setting from OSQP paper
density = 0.15; % 15% must be non-zero

% solver setup
maxiter=1000;
num_tol = 1e-8;

% loop setup
problem_sizes = 10:step:n_max;
l = size(problem_sizes,2);
residual = zeros(l,number_of_problems);
iterations = zeros(l,number_of_problems);
cpu_t = zeros(l,number_of_problems);
j = 1;

for n = problem_sizes
    [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density);
    [x0, ma] = QP_PrimalActiveSet_PAS_initial_point(H, g, A, b, C, dl, du, l, u, maxiter);
    
    tic;
    [x,fval,exitflag,output] = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,x0,options_as);
    cpu_t(j,1) = toc;
    iterations(j,1) = output.iterations;

    tic;
    [x,fval,exitflag,output] = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_ip);
    cpu_t(j,2) = toc;
    iterations(j,2) = output.iterations;
    
    Abar = [A C -C eye(length(l)) -eye(length(u))];
    bbar = [-b; dl; -du; l; -u];
    tic;
    [xstar, iterations(j,3)] = QP_primalActiveSet_core(H, g, Abar, bbar, x0, ma, maxiter, num_tol);
    cpu_t(j,3) = toc;
    
    tic;
    [xstar, iterations(j,4)] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,4) = toc;

    j = j + 1;
end


figure,
hold on

plot(problem_sizes,cpu_t(:,1),  'LineWidth', 3);
plot(problem_sizes,cpu_t(:,2), 'LineWidth', 3);
plot(problem_sizes,cpu_t(:,3), 'LineWidth', 3);
plot(problem_sizes,cpu_t(:,4), 'LineWidth', 3);
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("cpu time [s]")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/281_Speed_QP.png')
figure,
hold on

plot(problem_sizes,iterations(:,1),  'LineWidth', 5);
plot(problem_sizes,iterations(:,2), 'LineWidth', 3);
plot(problem_sizes,iterations(:,3), 'LineWidth', 3);
plot(problem_sizes,iterations(:,4), 'LineWidth', 3);
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("iterations")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/282_No_iterations_QP.png')

%%

% test setup
number_of_problems = 4; %Change this number as algorithms are added to the loop
n_max =200;
step = 10;
Variable_names = ["quadprogAS","quadprogIP","PrimalActiveSet","InteriorPoint"];

%QP setup
alpha = 0.01;
beta = 4; % setting from OSQP paper
density = 0.15; % 15% must be non-zero

% solver setup
maxiter=1000;
num_tol = 1e-8;

% loop setup
problem_sizes = 10:step:n_max;
l = size(problem_sizes,2);
residual = zeros(l,number_of_problems);
iterations = zeros(l,number_of_problems);
cpu_t = zeros(l,number_of_problems);
j = 1;

for n = problem_sizes
    [H,g,A,b,C,dl,du,l,u] = GeneratorBoxQP(n,alpha,beta,density);
    [x0, ma] = QP_PrimalActiveSet_PAS_initial_point(H, g, A, b, C, dl, du, l, u, maxiter);
    
    tic;
    [x,fval,exitflag,output] = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,x0,options_as);
    cpu_t(j,1) = toc;
    iterations(j,1) = output.iterations;

    tic;
    [x,fval,exitflag,output] = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_ip);
    cpu_t(j,2) = toc;
    iterations(j,2) = output.iterations;
    
    Abar = [A C -C eye(length(l)) -eye(length(u))];
    bbar = [-b; dl; -du; l; -u];
    tic;
    [xstar, iterations(j,3)] = QP_primalActiveSet_core(H, g, Abar, bbar, x0, ma, maxiter, num_tol);
    cpu_t(j,3) = toc;
    
    tic;
    [xstar, iterations(j,4)] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,4) = toc;

    j = j + 1;
end


figure,
hold on

plot(problem_sizes,cpu_t(:,1),  'LineWidth', 3);
plot(problem_sizes,cpu_t(:,2), 'LineWidth', 3);
plot(problem_sizes,cpu_t(:,3), 'LineWidth', 3);
plot(problem_sizes,cpu_t(:,4), 'LineWidth', 3);
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("cpu time [s]")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/283_Speed_boxQP.png')
figure,
hold on

plot(problem_sizes,iterations(:,1),  'LineWidth', 3);
plot(problem_sizes,iterations(:,2), 'LineWidth', 3);
plot(problem_sizes,iterations(:,3), 'LineWidth', 3);
plot(problem_sizes,iterations(:,4), 'LineWidth', 3);
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("iterations")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/284_No_iterations_boxQP.png')
%%

% test setup
number_of_problems = 4; %Change this number as algorithms are added to the loop
n_max =50;
step = 5;
Variable_names = ["quadprogAS","quadprogIP","PrimalActiveSet","InteriorPoint"];

%QP setup
alpha = 0.01;
beta = 4; % setting from OSQP paper
density = 0.15; % 15% must be non-zero

% solver setup
maxiter=1000;
num_tol = 1e-8;

% loop setup
problem_sizes = 10:step:n_max;
l = size(problem_sizes,2);
residual = zeros(l,number_of_problems);
iterations = zeros(l,number_of_problems);
cpu_t = zeros(l,number_of_problems);
j = 1;

for n = problem_sizes
    [H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP_bn(n,beta,density);
    [x0, ma] = QP_PrimalActiveSet_PAS_initial_point(H, g, A, b, C, dl, du, l, u, maxiter);
    
    tic;
    [x,fval,exitflag,output] = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,x0,options_as);
    cpu_t(j,1) = toc;
    iterations(j,1) = output.iterations;

    tic;
    [x,fval,exitflag,output] = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u,[],options_ip);
    cpu_t(j,2) = toc;
    iterations(j,2) = output.iterations;
    
    Abar = [A C -C eye(length(l)) -eye(length(u))];
    bbar = [-b; dl; -du; l; -u];
    tic;
    [xstar, iterations(j,3)] = QP_primalActiveSet_core(H, g, Abar, bbar, x0, ma, maxiter, num_tol);
    cpu_t(j,3) = toc;
    
    tic;
    [xstar, iterations(j,4)] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);
    cpu_t(j,4) = toc;

    j = j + 1;
end


figure,
hold on

plot(problem_sizes,cpu_t(:,1),  'LineWidth', 3);
plot(problem_sizes,cpu_t(:,2), 'LineWidth', 3);
plot(problem_sizes,cpu_t(:,3), 'LineWidth', 3);
plot(problem_sizes,cpu_t(:,4), 'LineWidth', 3);
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("cpu time [s]")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/285_Speed_HuberFitting.png')
figure,
hold on

plot(problem_sizes,iterations(:,1),  'LineWidth', 3);
plot(problem_sizes,iterations(:,2), 'LineWidth', 3);
plot(problem_sizes,iterations(:,3), 'LineWidth', 3);
plot(problem_sizes,iterations(:,4), 'LineWidth', 3);
%scatter(problem_sizes,residual(:,5), '*');
%scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend(Variable_names)
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("iterations")
%title("\beta = 4")

hold off

saveas(gcf,'./Plots/286_No_iterations_HuberFitting.png')