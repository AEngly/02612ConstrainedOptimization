
% ---------------- DESCRIPTION --------------
%
% Name: Driver for Problem 1
% Type: Execution File
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

%% Test numerical stability
alpha = 0.15;
density = 0.15; % 15% must be non-zero

n = 500;
step = 10;

problem_sizes = 10:step:n;

l = size(problem_sizes,2);

%% gather data beta 0.5

beta = 0.5;

residual = zeros(l,6);
j = 1;

for i = problem_sizes
    [H,g,A,b,x,lambda] = updated_GeneratorECQP(i,alpha,beta,density);
    [xstar] = EqualityQPSolverLUdense(H, g, A, b);
    residual(j,1) = norm(x-xstar);
    [xstar] = EqualityQPSolverLUsparse(H, g, A, b);
    residual(j,2) = norm(x-xstar);
    [xstar] = EqualityQPSolverLDLdense(H, g, A, b);
    residual(j,3) = norm(x-xstar);
    [xstar] = EqualityQPSolverLDLsparse(H, g, A, b);
    residual(j,4) = norm(x-xstar);
    [xstar] = EqualityQPSolverRS(H, g, A, b);
    residual(j,5) = norm(x-xstar);
    [xstar] = EqualityQPSolverNS(H, g, A, b);
    residual(j,6) = norm(x-xstar);
    j = j + 1;
end



%% plot numerical performance


figure,
hold on

scatter(problem_sizes,residual(:,1),200, '.','LineWidth',2);
scatter(problem_sizes,residual(:,2), 'o','LineWidth',2);
scatter(problem_sizes,residual(:,3), 'x','LineWidth',2);
scatter(problem_sizes,residual(:,4), '+','LineWidth',2);
scatter(problem_sizes,residual(:,5), '*','LineWidth',2);
scatter(problem_sizes,residual(:,6), 's','LineWidth',2);

set(gca,'yscale','log')
legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','southeast')
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("residual")
title("\beta = 0.5")

hold off

saveas(gcf,'./Plots/1511_Numerical_performance_Solvers_beta_0_5.png')

%%

alpha = 0.15;
density = 0.15; % 15% must be non-zero

n = 50;
step = 1;

problem_sizes = 10:step:n;
betas = 0.05:0.05:0.95;

l = size(problem_sizes,2);
m = size(betas, 2);

residuals = zeros(m,l,6);

k = 1;

for beta = betas
    j = 1;
    for i = problem_sizes
        [H,g,A,b,x,lambda] = updated_GeneratorECQP(i,alpha,beta,density);
    [xstar] = EqualityQPSolverLUdense(H, g, A, b);
    residuals(k,j,1) = norm(x-xstar);
    [xstar] = EqualityQPSolverLUsparse(H, g, A, b);
    residuals(k,j,2) = norm(x-xstar);
    [xstar] = EqualityQPSolverLDLdense(H, g, A, b);
    residuals(k,j,3) = norm(x-xstar);
    [xstar] = EqualityQPSolverLDLsparse(H, g, A, b);
    residuals(k,j,4) = norm(x-xstar);
    [xstar] = EqualityQPSolverRS(H, g, A, b);
    residuals(k,j,5) = norm(x-xstar);
    [xstar] = EqualityQPSolverNS(H, g, A, b);
    residuals(k,j,6) = norm(x-xstar);
        j = j + 1;
    end
    disp(k)
    k = k + 1;
end

%% Plot results
titles = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','southeast'};

limits = [1e-15,1e-12];
figure,
for i=1:6
    subplot(6,1,i)
    pcolor(problem_sizes,betas,residuals(:,:,i))
    axis ij
    xlabel("n")
    ylabel("betas")
    colorbar
    set(gca,'clim',limits([1,end]))
    title(titles(i))
end

saveas(gcf,'./Plots/1512_Numerical_issues.png')

%% Benchmarking solvers
alpha = 0.15;
density = 0.15; % 15% must be non-zero

n = 500;
step = 10;

smoother = 10;
problem_sizes = 50:step:n;

l = size(problem_sizes,2);

%% beta 0.05

beta = 0.05;

[TTC_avg1] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg1(i,:), 'LineWidth', 3);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("Problemsize - n")
    ylabel("CPU time [s]")
    title("\beta 0.05")


hold off

saveas(gcf,'./Plots/1521_ComparisonSolvers_beta_0_05.png')


% beta 0.2

beta = 0.2;

[TTC_avg2] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg2(i,:), 'LineWidth', 3);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("Problemsize - n")
    ylabel("CPU time [s]")
    title("\beta 0.20")


hold off

saveas(gcf,'./Plots/1522_ComparisonSolvers_beta_0_2.png')
% beta 0.5

beta = 0.5;

[TTC_avg3] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg3(i,:), 'LineWidth', 3);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("Problemsize - n")
    ylabel("CPU time [s]")
    title("\beta 0.50")


hold off

saveas(gcf,'./Plots/1523_ComparisonSolvers_beta_0_5.png')
% beta 0.8

beta = 0.8;

[TTC_avg4] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg4(i,:), 'LineWidth', 3);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("Problemsize - n")
    ylabel("CPU time [s]")
    title("\beta 0.80")


hold off

saveas(gcf,'./Plots/1524_ComparisonSolvers_beta_0_8.png')
% beta 0.95
beta = 0.95;

[TTC_avg5] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg5(i,:), 'LineWidth', 3);
    end
    
    %set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', 'Location','northwest')
    xlabel("Problemsize - n")
    ylabel("CPU time [s]")
    title("\beta 0.95")

hold off

saveas(gcf,'./Plots/1525_ComparisonSolvers_beta_0_95.png')

%% data gathering loop

function [TTC_avg] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density) 
    TTC = zeros(6,l,smoother);
    TTC_avg = zeros(6,l);
    j = 1;

    for i = problem_sizes
    
        % Display
        %fprintf('Problem size: %d\n', i);
    
       
    
        for k = 1:smoother

            [H,g,A,b] = updated_GeneratorECQP(i,alpha,beta,density);
    
            TTC(1,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
            TTC(2,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
            TTC(3,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
            TTC(4,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
            TTC(5,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
            TTC(6,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        end
    
        TTC_avg(1,j) = mean(TTC(1,j,:));
        TTC_avg(2,j) = mean(TTC(2,j,:));
        TTC_avg(3,j) = mean(TTC(3,j,:));
        TTC_avg(4,j) = mean(TTC(4,j,:));
        TTC_avg(5,j) = mean(TTC(5,j,:));
        TTC_avg(6,j) = mean(TTC(6,j,:));
    
        j = j + 1;
    
    end
end