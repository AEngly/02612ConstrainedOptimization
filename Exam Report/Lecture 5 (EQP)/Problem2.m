%% 2.1)

%% Problem 2 - Linear Optimization

x = -10:0.05:10;
y = -10:0.05:10;
[X,Y] = meshgrid(x,y);

F = (X - 1).^2 + (Y - 2.5).^2;

v = -20:2:20;
[c,h]=contour(X,Y,F,v,"linewidth",2);
colorbar

yc1 = (x + 2)./2; % larger than equal ... x >= 2y - 2
yc2 = (x - 6)./(-2); % y less ...
yc3 = (x - 2)./2; % y greater ...
xc4 = 0; % x larger than 0
yc5 = 0; % y larger than 0

hold on
    fill([x(1),repelem(xc4, 401), x(1)],[y(1), y, y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x(1), x, x(end)],[y(1), repelem(yc5, 401), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc1, y(end), x(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc2, y(end), x(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x, x(end), x(1)], [yc3, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7)
    %plot(2.5,4.5,'r.', 'MarkerSize',20)
hold off

%legend('','','','','','','Minimizer')

xlim([-10 10])
ylim([-10 10])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')

saveas(gcf,'./Lecture 5 (EQP)/ContourProblem21.png')

%% 2.4)

H = [2 0; 0 2];
g = [-2; -5];

A_all = [1 -2;
     -1 -2;
     -1 2;
     1 0;
     0 1];

A_all = A_all';
[n_row, n_col] = size(A_all);

b = [-2; -6; -2; 0; 0];

% Test solution

x0 = [2;0];
active_constraints = find(~(A_all' * x0 - b));

% The iterates and the solution are not necessarily on the boundary of the
% feasible region (as opposed to Simplex).

% Three variants:
% 1) Primal
% 2) Dual
% 3) Primal-dual

% We change the working set to be only x2 >= 0.

% NEW ITERATION

Wk = [5];
A_all_t = A_all';
A_Wk = A_all_t(Wk,:);
b_Wk = b(Wk);

[LHS, RHS] = KKT_matrix(H, g, A_Wk', b_Wk);

% Then we can solve it

x = LHS \ RHS;
lambda = x(3:end);
x = x(1:2);


%% Testing active set algorithms

H = [2 0; 0 2];
g = [-2; -5];
A = [1 -2;
     -1 -2;
     -1 2;
     1 0;
     0 1];
A = A';
b = [-2; -6; -2; 0; 0];
x0 = [2;0];

x_final = ActiveSetMethodQP(H, g, A, b, x0);

