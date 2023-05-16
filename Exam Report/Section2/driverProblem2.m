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

%% Test minimal QP

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
saveas(gcf,'./Plots/2_5A_ActiveSetTour.png')

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

%% Test minimal QP

maxiter = 100;
[H,g,A,b,C,dl,du,l,u] = MinimalQP();
[xas, iter, x_k] = QP_PrimalActiveSet_linprog(H, g, A, b, C, dl, du, l, u, maxiter);
xstar = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);
norm(xas-xstar)

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
saveas(gcf,'./Plots/2_5A_ActiveSetTour.png')
%% 2.8) Compare implementations with quadprog (MOSEK, Gurobi, qpOASES, OSQP, YALMIP and cvx)

% The test should be on some general QP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics from scalable problems:
% 1) Size of problem (n) on x-axis, and time-to-completion on y-axis.

