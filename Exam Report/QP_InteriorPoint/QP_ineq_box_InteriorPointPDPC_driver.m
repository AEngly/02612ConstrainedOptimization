% Todo
% Add known test problem
% Test functionality
% Contour plot of path
% Performance plot

%% Add test problem path
addpath(genpath('../TestTools'))

%% Generate know test problem

%% Path on Contour plot

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

%% Generate test random test problem
n = 20;
alpha = 0.1;
density = 0.15;

[H,g,bl,A,bu,l,u] = RandomQP_ineq_box(n,alpha,density);

x0 = zeros(n,1);

%% find reference solution
Aq = [full(-A) full(A)]';
bq = [-bl; bu];
xr = quadprog(H,g,Aq,bq,[],[],l,u);

%%
z0 = ones(2*n*2+2*n,1);
s0 = ones(2*n*2+2*n,1);
[xi,z,s] = QP_ineq_box_InteriorPointPDPC(H,g,bl,A,bu,l,u,x0,z0,s0);
xi