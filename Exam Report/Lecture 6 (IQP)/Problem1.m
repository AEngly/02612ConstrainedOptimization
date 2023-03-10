%% 1.1) Contour plot

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

saveas(gcf,'./ContourProblem11.png')

%% 1.2) Testing primal active-set QP

G = [2 0; 0 2];
g = [-2; -5];

A_t = [1 -2;
     -1 -2;
     -1 2;
     1 0;
     0 1];

A = A_t';
[n_row, n_col] = size(A);

b = [-2; -6; -2; 0; 0];

% Test solution

x0 = [2;0];

[xk, mu_star, active_constraints] = primalActiveSetQP(G, g, A, b, x0);

% Then we can plot the points on the contours
x = -10:0.05:10;
y = -10:0.05:10;
[X,Y] = meshgrid(x,y);

F = (X - 1).^2 + (Y - 2.5).^2;

v = -20:2:20;
[c,h]=contour(X,Y,F,v,"linewidth",2, 'DisplayName','Contour');
colorbar

yc1 = (x + 2)./2; % larger than equal ... x >= 2y - 2
yc2 = (x - 6)./(-2); % y less ...
yc3 = (x - 2)./2; % y greater ...
xc4 = 0; % x larger than 0
yc5 = 0; % y larger than 0

hold on
    fill([x(1),repelem(xc4, 401), x(1)],[y(1), y, y(end)],[0.7 0.7 0.7],"facealpha",0.7)
    fill([x(1), x, x(end)],[y(1), repelem(yc5, 401), y(1)],[0.7 0.7 0.7],"facealpha",0.7,'DisplayName','c_{2}')
    fill([x, x(end), x(1)], [yc1, y(end), x(end)],[0.7 0.7 0.7],"facealpha",0.7,'DisplayName','c_{3}')
    fill([x, x(end), x(1)], [yc2, y(end), x(end)],[0.7 0.7 0.7],"facealpha",0.7,'DisplayName','c_{4}')
    fill([x, x(end), x(1)], [yc3, y(1), y(1)],[0.7 0.7 0.7],"facealpha",0.7,'DisplayName','c_{5}')
    plot(xk(1,1),xk(2,1),'.','color','red', 'MarkerSize',16, 'DisplayName','x_{0}')
    plot(xk(1,2),xk(2,2),'*','color','red', 'MarkerSize',16, 'DisplayName','x_{1}')
    plot(xk(1,3),xk(2,3),'+','color','red', 'MarkerSize',16, 'DisplayName','x_{2}')
    plot(xk(1,4),xk(2,4),'o','color','red', 'MarkerSize',16, 'DisplayName','x_{3}')
    plot(xk(1,5),xk(2,5),'square','color','red', 'MarkerSize',16, 'DisplayName','x_{4}')
    plot(xk(1,6),xk(2,6),'^','color','red', 'MarkerSize',16, 'DisplayName','x_{5}')
    plot(xk(1,7),xk(2,7),'diamond','color','red', 'MarkerSize',16, 'DisplayName','x_{6}')
hold off

%legend('','','','','','','Minimizer')

xlim([-10 10])
ylim([-10 10])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('$x_{1}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('$x_{2}$','interpreter','latex', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')
legend()
saveas(gcf,'./ActiveSetTour.png')
