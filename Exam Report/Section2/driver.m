
% ---------------- DESCRIPTION --------------
%
% Name: Driver for Problem 2
% Type: Execution File
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
% ---------------- IMPLEMENTATION --------------


%% Test problems

% Initial test problem

maxiter = 100;
H = eye(2);
g = [-2;-5];
A = [1;-1];
b = zeros(1,1);
C = [ 1 -1; -2 -2];
dl = [-2; -6];
du = [2; -2];
l = zeros(2,1);
u = zeros(2,0);

[xas] = test(H, g, A, b, C, dl, du, l, u, maxiter);

xstar = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

norm(xas-xstar)

%% RandomQP settings from OSQP paper



%% 2.5 test Primal-active set
alpha = 0.01;
beta = 2; % setting from OSQP paper
density = 0.15; % 15% must be non-zero
n = 2;
maxiter=10000;
[H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density);

[all_xk] = test(H, g, A, b, C, dl, du, l, u, maxiter);

xstarx = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

x = all_xk(:,end);
norm(x-xstarx)

%% test huber fitting
n=5;
beta = 100;
maxiter=1000;
density = 0.15;

[H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP(n,beta,density);

[all_xk]  = test(H, g, A, b, C, dl, du, l, u,maxiter);

xstarx = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

x = all_xk(:,end);
norm(x-xstarx)
%% test portfolio

k = 10;
gamma = 1;
density = 0.50;
maxiter=100;

[H,g,A,b,C,dl,du,l,u] = GeneratorPortfolioOptimizationQP(k,gamma,density);

[x] = test(H,g,A,b,C,dl,du,l,u,maxiter);

xstarx = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

norm(x-xstarx)
%% 2.7 test for InteriorPoint

%% Test problems

% Initial test problem

maxiter = 100;
H = eye(2);
g = [-2;-5];
A = [1;-1];
b = zeros(1,1);
C = [ 1 -1; -2 -2];
dl = [-2; -6];
du = [2; -2];
l = zeros(2,1);
u = zeros(2,0);

[xas] = QP_InteriorPointPDPC(H, g, A, b, C, dl, du, l, u, maxiter);

xstar = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

norm(xas-xstar)

%% 
alpha = 0.01;
beta = 2; % setting from OSQP paper
density = 0.15; % 15% must be non-zero
n = 2;

maxiter=30;
[H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density);

[x,y,z,s] = QP_InteriorPointPDPC(H,g,A,b,C,dl,du,l,u,maxiter);

xstarx = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

norm(x-xstarx)

%% test for random qp

alpha = 0.01;
density = 0.15; % 15% must be non-zero
n = 100;

[H,g,A,b,C,dl,du,l,u] = RandomQP_(n,alpha,density);

[x,y,z,s] = QP_InteriorPointPDPC(H,g,A,b,C,dl,du,l,u,maxiter);

xstarx = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

norm(x-xstarx)

%% test huber fitting
n=2;
beta = 10;
maxiter=100;
density = 0.15;

[H,g,A,b,C,dl,du,l,u] = GeneratorHuberFittingQP(n,beta,density);
[x,y,z,s] = QP_InteriorPointPDPC(H,g,A,b,C,dl,du,l,u,maxiter);

xstarx = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

norm(x-xstarx)
%% test portfolio

k = 10;
gamma = 1;
density = 0.50;
maxiter=100;

[H,g,A,b,C,dl,du,l,u] = GeneratorPortfolioOptimizationQP(k,gamma,density);

[x,y,z,s] = QP_InteriorPointPDPC(H,g,A,b,C,dl,du,l,u,maxiter);

xstarx = quadprog(H,g,[-C'; C'],[-dl;du],A',-b,l,u);

norm(x-xstarx)

%%
% Options for quadprog
options =  optimoptions('quadprog','Display','off');

n_max = 310;
step = 20;

smoother = 3;
TTC = zeros(6,n,smoother);
TTC_avg = zeros(6,n);
j = 1;

problem_sizes = 10:step:n_max;

for n = problem_sizes

    % Display
    fprintf('Problem size: %d\n', n);

    [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density);

    for k = 1:smoother

        TTC(1,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
        TTC(2,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
        TTC(3,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
        TTC(4,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
        TTC(5,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
        TTC(6,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        TTC(7,j,k) = cpuTimer(@quadprog, H, g, A', b, [], [], [], [], [], options);

    end

    TTC_avg(1,j) = mean(TTC(1,j,:));
    TTC_avg(2,j) = mean(TTC(2,j,:));
    TTC_avg(3,j) = mean(TTC(3,j,:));
    TTC_avg(4,j) = mean(TTC(4,j,:));
    TTC_avg(5,j) = mean(TTC(5,j,:));
    TTC_avg(6,j) = mean(TTC(6,j,:));
    TTC_avg(7,j) = mean(TTC(7,j,:));

    j = j + 1;

end

%% Plotting the comparison (with quadprog)

figure,
hold on

    for i=1:7
        plot(problem_sizes, TTC_avg(i,:));
    end
    
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'quadprog')
    xlabel("n")
    ylabel("CPU time")

hold off

saveas(gcf,'./Plots/ComparisonSolvers1.png')

%% Plotting the comparison (without quadprog)

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg(i,:));
    end
    
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace')
    xlabel("n")
    ylabel("CPU time")

hold off

saveas(gcf,'./Plots/ComparisonSolvers2.png')

%% Plotting the comparison (with quadprog)

plot(problem_sizes, TTC_avg(6,:));
legend('RangeSpace')
xlabel("n")
ylabel("CPU time")

saveas(gcf,'./Plots/ComparisonSolversRangeSpace.png')
