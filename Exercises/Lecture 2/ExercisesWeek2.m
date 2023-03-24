
%% Problem 2 - Linear Optimization

g = [1 -2]';
A = [1 0 1 1 -5; 0 1 -1 -5 1];
b = [0 0 -2 -20 -15];

x = -10:0.05:10;
y = -10:0.05:10;
[X,Y] = meshgrid(x,y);

F = g(1)*X + g(2)*Y;

v = [-20:2:20];
[c,h]=contour(X,Y,F,v,"linewidth",2);
colorbar

xc1 = 0; % x larger than 0
yc2 = 0; % y larger than 0
yc3 = x + 2; % y less ...
yc4 = (x + 20)/5; % y less ...
yc5 = (5*x - 15); % y greater ...



hold on
    fill([x(1),repelem(xc1,401), x(1)],[y(1), y, y(end)],[0.7 0.7 0.7],"facealpha",0.5)
    fill([x(1), x, x(end)],[y(1), repelem(yc2, 401), y(1)],[0.7 0.7 0.7],"facealpha",0.5)
    fill([x(1),x,x(end)],[y(end),yc3,y(end)],[0.7 0.7 0.7],"facealpha",0.5)
    fill([x(1),x,x(end)],[y(end),yc4,y(end)],[0.7 0.7 0.7],"facealpha",0.5)
    fill([x(end),x,x(end)],[y(1),yc5,y(end)],[0.7 0.7 0.7],"facealpha",0.5)
    plot(2.5,4.5,'r.', 'MarkerSize',20)
hold off

legend('','','','','','','Minimizer')

xlim([-10 10])
ylim([-10 10])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('x', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('y', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')

saveas(gcf,'./Lecture 2/ContourPlot1_minimizer.png')

%% Problem 3 - Nonlinear Optimization

g = [1 -2]';
A = [1 0 1 1 -5; 0 1 -1 -5 1];
b = [0 0 -2 -20 -15];

x = -10:0.05:10;
y = -10:0.05:10;
[X,Y] = meshgrid(x,y);

F = X.^2 + Y.^2 + 3*Y;

v = [-20:2:20];
[c,h]=contour(X,Y,F,v,"linewidth",2);
colorbar

c1 = @(x,y) x.^2 + (1+ y.^2) - 1;

hold on
    fill(y, xc1,[0.7 0.7 0.7],"facealpha",0.5)
    fimplicit(c1,[-10 10 -10 10])
    plot(2.5,4.5,'r.', 'MarkerSize',20)
hold off

legend('','','','','','','Minimizer')

xlim([-10 10])
ylim([-10 10])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('x', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('y', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')

saveas(gcf,'./Lecture 2/ContourPlotProblem3.png')

%% Problem 3 - Nonlinear Optimization (Part 2)

g = [1 -2]';
A = [1 0 1 1 -5; 0 1 -1 -5 1];
b = [0 0 -2 -20 -15];

x = -10:0.05:10;
y = -10:0.05:10;
[X,Y] = meshgrid(x,y);

F = X.^2 + Y.^2 + 3*Y;

v = [-20:2:20];
[c,h]=contour(X,Y,F,v,"linewidth",2);
colorbar

c1 = @(x,y) x.^2 + (1+ y.^2) - 1;

hold on
    fill(y, xc1,[0.7 0.7 0.7],"facealpha",0.5)
    fimplicit(c1,[-10 10 -10 10])
    plot(2.5,4.5,'r.', 'MarkerSize',20)
hold off

legend('','','','','','','Minimizer')

xlim([-10 10])
ylim([-10 10])
%title('Contour Plot of Linear Program', 'FontSize',20)
xlabel('x', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold') 
ylabel('y', 'FontSize',16,'Interpreter','LaTeX','Color','black','FontWeight','bold')

saveas(gcf,'./Lecture 2/ContourPlotProblem3.png')