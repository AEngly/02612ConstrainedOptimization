%% Plotting Himmelblau

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


saveas(gcf,'./Plots/ContourHB.png')

%% Solution to Himmelblau (fmincon)

% Starting point
x0 = [-3.2; 0];

% Specify starting condition
A = [-4 10; 4 -10];
b = [c2_u; -c2_l];
Aeq = [];
beq = [];
lb = [x1_l; x2_l];
ub = [x1_u; x2_u];

% Specify constraints
l = [-5 5];
u = [5 5];
cl = [0;0];
cu = [47;70];

% Solve with fmincon
optimoptions('fmincon','Display', 'iter', 'SpecifyObjectiveGradient', true, 'TolFun', 1e-3);
xfmin = fmincon(@objfminconHimmelblau, x0, A, b, Aeq, beq, lb, ub, @confunHimmelblau);
disp(xfmin)

%% 

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

saveas(gcf,'./Plots/ContourHB.png')