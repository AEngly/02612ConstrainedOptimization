% ---------------- DESCRIPTION --------------
%
% Name: driverProblem1.m
% Type: Driver for all testing and plots from problem 1 (ECQP)
%
%
% Problem structure:
%           min     1/2x'Hx + g'x
%            x
%           s.t.    A'x + b = 0
%
% Created: 12.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark

% ---------------- IMPLEMENTATION --------------

%% 1.5.1) Testing of random generator

alpha = 0.15;
density = 0.15; % 15% must be non-zero
n = 100;
beta = 0.5;

[H,g,A,b,x,l1] = GeneratorECQP(n,alpha,beta,density);
[x1, lambda1] = EqualityQPSolver(H, g, A, b, "LUdense");
[x2, lambda2] = EqualityQPSolver(H, g, A, b, "LUsparse");
[x3, lambda3] = EqualityQPSolver(H, g, A, b, "LDLdense");
[x4, lambda4] = EqualityQPSolver(H, g, A, b, "LDLsparse");
[x5, lambda5] = EqualityQPSolver(H, g, A, b, "RangeSpace");
[x6, lambda6] = EqualityQPSolver(H, g, A, b, "NullSpace");
options = optimoptions('quadprog','Display','off');
[x7,fval,exitflag,output,lambda7] = quadprog(H,g,[],[],A',-b, [], [], [], options);

if norm(x1-x,2) < 1e-09
    fprintf("LUdense == quadprog\n");
end
if norm(x2-x,2) < 1e-09
    fprintf("LUsparse == quadprog\n");
end
if norm(x3-x,2) < 1e-09
    fprintf("LDLdense == quadprog\n");
end
if norm(x4-x,2) < 1e-09
    fprintf("LDLsparse == quadprog\n");
end
if norm(x5-x,2) < 1e-09
    fprintf("RangeSpace == quadprog\n");
end
if norm(x6-x,2) < 1e-09
    fprintf("NullSpace == quadprog\n");
end

% Problems with LU dense and sparse -> problem specific


%% 1.5.1) Test of Correctness and Stability

alpha = 0.15;
density = 0.15; % 15% must be non-zero
n = 300;
step = 10;
problem_sizes = 50:step:n;
l = size(problem_sizes,2);

% PROBLEM SPECIFICS
beta = 0.5;

% PREPARE PLACEHOLDER
res_x = zeros(l,6);
res_l = zeros(l,6);
j = 1;

for i = problem_sizes
    [H,g,A,b,x,l1] = GeneratorECQP(i,alpha,beta,density);
    [xstar,l] = EqualityQPSolver(H, g, A, b, "LUdense");
    res_x(j,1) = norm(x-xstar);
    res_l(j,1) = norm(l-l1);
    [xstar,l] = EqualityQPSolver(H, g, A, b, "LUsparse");
    res_x(j,2) = norm(x-xstar);
    res_l(j,2) = norm(l-l1);
    [xstar,l] = EqualityQPSolver(H, g, A, b, "LDLdense");
    res_x(j,3) = norm(x-xstar);
    res_l(j,3) = norm(l-l1);
    [xstar,l] = EqualityQPSolver(H, g, A, b, "LDLsparse");
    res_x(j,4) = norm(x-xstar);
    res_l(j,4) = norm(l-l1);
    [xstar,l] = EqualityQPSolver(H, g, A, b, "RangeSpace");
    res_x(j,5) = norm(x-xstar);
    res_l(j,5) = norm(l-l1);
    [xstar,l] = EqualityQPSolver(H, g, A, b, "NullSpace");
    res_x(j,6) = norm(x-xstar);
    res_l(j,6) = norm(l-l1);
    j = j + 1;
end

res_x = log(res_x);
res_l = log(res_l);

hold on
scatter(problem_sizes,res_x(:,1),200, '.','LineWidth',2);
scatter(problem_sizes,res_x(:,2), 'o','LineWidth',2);
scatter(problem_sizes,res_x(:,3), 'x','LineWidth',2);
scatter(problem_sizes,res_x(:,4), '+','LineWidth',2);
scatter(problem_sizes,res_x(:,5), '*','LineWidth',2);
scatter(problem_sizes,res_x(:,6), 's','LineWidth',2);
legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','southeast')
xlabel("n");
ylabel("log(residuals)");
hold off
saveas(gcf,'./Plots/1511_Residuals_x.png')
close;

hold on
scatter(problem_sizes,res_l(:,1),200, '.','LineWidth',2);
scatter(problem_sizes,res_l(:,2), 'o','LineWidth',2);
scatter(problem_sizes,res_l(:,3), 'x','LineWidth',2);
scatter(problem_sizes,res_l(:,4), '+','LineWidth',2);
scatter(problem_sizes,res_l(:,5), '*','LineWidth',2);
scatter(problem_sizes,res_l(:,6), 's','LineWidth',2);
legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','southeast')
xlabel("n");
ylabel("log(residuals)");
hold off
saveas(gcf,'./Plots/1511_Residuals_lambda.png')
close;

%% 1.5.2) Speed Testing with Random Generator (beta = 0.15)

options = optimoptions('quadprog','Display','off');

beta = 0.15;
alpha = 0.15;
density = 0.15;
n = 30;
step = 25;
smoother = 2;
TTC = zeros(7,smoother,n);
problem_sizes = step:step:n*step;

for k = 1:smoother
    j = 0;
    for i = problem_sizes

        % Generate problem
        [H,g,A,b,x,lambda] = GeneratorECQP(i,alpha,beta,density);
    
        % Display
        fprintf('Problem size: %d\n', i);
        j = j + 1;
    
        TTC(1,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
        TTC(2,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
        TTC(3,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
        TTC(4,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
        TTC(5,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
        TTC(6,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        TTC(7,k,j) = cpuTimer(@quadprog, H,g,[],[],A',-b, [], [], [], options);
    
    end
end

hold on

    for i=1:7
        plot(problem_sizes, log10(reshape(mean(TTC(i,:,:),2), [1,j])), 'LineWidth', 2);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("n")
    ylabel("CPU time [log(s)]")

hold off

saveas(gcf,'./Plots/152_PerformanceBeta015.png');
close;

%% 1.5.2) Speed Testing with Random Generator (beta = 0.30)

options = optimoptions('quadprog','Display','off');

beta = 0.30;
alpha = 0.15;
density = 0.15;
n = 30;
step = 25;
smoother = 2;
TTC = zeros(7,smoother,n);
problem_sizes = step:step:n*step;

for k = 1:smoother
    j = 0;
    for i = problem_sizes

        % Generate problem
        [H,g,A,b,x,lambda] = GeneratorECQP(i,alpha,beta,density);
    
        % Display
        fprintf('Problem size: %d\n', i);
        j = j + 1;
    
        TTC(1,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
        TTC(2,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
        TTC(3,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
        TTC(4,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
        TTC(5,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
        TTC(6,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        TTC(7,k,j) = cpuTimer(@quadprog, H,g,[],[],A',-b, [], [], [], options);
    
    end
end

hold on

    for i=1:7
        plot(problem_sizes, log10(reshape(mean(TTC(i,:,:),2), [1,j])), 'LineWidth', 2);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("n")
    ylabel("CPU time [log(s)]")

hold off

saveas(gcf,'./Plots/152_PerformanceBeta030.png');
close;

%% 1.5.2) Speed Testing with Random Generator (beta = 0.50)

options = optimoptions('quadprog','Display','off');

beta = 0.50;
alpha = 0.15;
density = 0.15;
n = 30;
step = 25;
smoother = 2;
TTC = zeros(7,smoother,n);
problem_sizes = step:step:n*step;

for k = 1:smoother
    j = 0;
    for i = problem_sizes

        % Generate problem
        [H,g,A,b,x,lambda] = GeneratorECQP(i,alpha,beta,density);
    
        % Display
        fprintf('Problem size: %d\n', i);
        j = j + 1;
    
        TTC(1,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
        TTC(2,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
        TTC(3,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
        TTC(4,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
        TTC(5,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
        TTC(6,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        TTC(7,k,j) = cpuTimer(@quadprog, H,g,[],[],A',-b, [], [], [], options);
    
    end
end

hold on

    for i=1:7
        plot(problem_sizes, log10(reshape(mean(TTC(i,:,:),2), [1,j])), 'LineWidth', 2);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("n")
    ylabel("CPU time [log(s)]")

hold off

saveas(gcf,'./Plots/152_PerformanceBeta050.png');
close;

%% 1.5.2) Speed Testing with Random Generator (beta = 0.75)

options = optimoptions('quadprog','Display','off');

beta = 0.75;
alpha = 0.15;
density = 0.15;
n = 30;
step = 25;
smoother = 2;
TTC = zeros(7,smoother,n);
problem_sizes = step:step:n*step;

for k = 1:smoother
    j = 0;
    for i = problem_sizes

        % Generate problem
        [H,g,A,b,x,lambda] = GeneratorECQP(i,alpha,beta,density);
    
        % Display
        fprintf('Problem size: %d\n', i);
        j = j + 1;
    
        TTC(1,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
        TTC(2,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
        TTC(3,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
        TTC(4,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
        TTC(5,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
        TTC(6,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        TTC(7,k,j) = cpuTimer(@quadprog, H,g,[],[],A',-b, [], [], [], options);
    
    end
end

hold on

    for i=1:7
        plot(problem_sizes, log10(reshape(mean(TTC(i,:,:),2), [1,j])), 'LineWidth', 2);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("n")
    ylabel("CPU time [log(s)]")

hold off

saveas(gcf,'./Plots/152_PerformanceBeta075.png');
close;

%% 1.5.2) Speed Testing with Random Generator (beta = 0.90)

options = optimoptions('quadprog','Display','off');

beta = 0.90;
alpha = 0.15;
density = 0.15;
n = 30;
step = 25;
smoother = 2;
TTC = zeros(7,smoother,n);
problem_sizes = step:step:n*step;

for k = 1:smoother
    j = 0;
    for i = problem_sizes

        % Generate problem
        [H,g,A,b,x,lambda] = GeneratorECQP(i,alpha,beta,density);
    
        % Display
        fprintf('Problem size: %d\n', i);
        j = j + 1;
    
        TTC(1,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
        TTC(2,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
        TTC(3,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
        TTC(4,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
        TTC(5,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
        TTC(6,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        TTC(7,k,j) = cpuTimer(@quadprog, H,g,[],[],A',-b, [], [], [], options);
    
    end
end

hold on

    for i=1:7
        plot(problem_sizes, log10(reshape(mean(TTC(i,:,:),2), [1,j])), 'LineWidth', 2);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("n")
    ylabel("CPU time [log(s)]")

hold off

saveas(gcf,'./Plots/152_PerformanceBeta090.png');
close;

%% 1.6.1) Testing of Recycle Problem

[H, g, A, b] = RecycleSystem(20);
[x1, lambda1] = EqualityQPSolver(H, g, A, b, "LUdense");
[x2, lambda2] = EqualityQPSolver(H, g, A, b, "LUsparse");
[x3, lambda3] = EqualityQPSolver(H, g, A, b, "LDLdense");
[x4, lambda4] = EqualityQPSolver(H, g, A, b, "LDLsparse");
[x5, lambda5] = EqualityQPSolver(H, g, A, b, "RangeSpace");
[x6, lambda6] = EqualityQPSolver(H, g, A, b, "NullSpace");
options = optimoptions('quadprog','Display','off');
[x7,fval,exitflag,output,lambda7] = quadprog(H,g,[],[],A',-b, [], [], [], options);

if norm(x1-x7,2) < 1e-09
    fprintf("LUdense == quadprog\n");
end
if norm(x2-x7,2) < 1e-09
    fprintf("LUsparse == quadprog\n");
end
if norm(x3-x7,2) < 1e-09
    fprintf("LDLdense == quadprog\n");
end
if norm(x4-x7,2) < 1e-09
    fprintf("LDLsparse == quadprog\n");
end
if norm(x5-x7,2) < 1e-09
    fprintf("RangeSpace == quadprog\n");
end
if norm(x6-x7,2) < 1e-09
    fprintf("NullSpace == quadprog\n");
end


%% 1.6.2) Testing of Recycle Problem

options = optimoptions('quadprog','Display','off');

n = 40;
step = 25;
smoother = 3;
TTC = zeros(7,3,n);
problem_sizes = step:step:n*step;

for k = 1:smoother
    j = 0;
    for i = problem_sizes
    
        % Display
        fprintf('Problem size: %d\n', i);
        j = j + 1;
    
        [H, g, A, b] = RecycleSystem(i);
        TTC(1,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
        TTC(2,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
        TTC(3,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
        TTC(4,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
        TTC(5,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
        TTC(6,k,j) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        TTC(7,k,j) = cpuTimer(@quadprog, H,g,[],[],A',-b, [], [], [], options);
    
    end
end

hold on

    for i=1:7
        plot(problem_sizes, log10(reshape(mean(TTC(i,:,:),2), [1,j])), 'LineWidth', 2);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("n")
    ylabel("CPU time [log(s)]")

hold off

saveas(gcf,'./Plots/161_RecycleSystem.png');
close;


