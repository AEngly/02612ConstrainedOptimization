
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

n = 300;
step = 10;

problem_sizes = 10:step:n;

l = size(problem_sizes,2);

%% gather data beta 0.5

beta = 0.5;

residual = zeros(l,2);
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

scatter(problem_sizes,residual(:,1), '.');
scatter(problem_sizes,residual(:,2), 'o');
scatter(problem_sizes,residual(:,3), 'x');
scatter(problem_sizes,residual(:,4), '+');
scatter(problem_sizes,residual(:,5), '*');
scatter(problem_sizes,residual(:,6), 's');

set(gca,'yscale','log')
legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace')
%legend('NullSpace', 'RangeSpace')
xlabel("n")
ylabel("residual")
title("beta 0.5")

hold off

saveas(gcf,'./Plots/Numerical_performance_Solvers_beta_0_5.png')

%% Benchmarking solvers
alpha = 0.15;
density = 0.15; % 15% must be non-zero

n = 300;
step = 10;

smoother = 10;
problem_sizes = 10:step:n;

l = size(problem_sizes,2);

%% beta 0.05

beta = 0.05;

[TTC_avg1] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg1(i,:));
    end
    
    set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace')
    xlabel("n")
    ylabel("CPU time")
    title("beta 0.05")

hold off

saveas(gcf,'./Plots/ComparisonSolvers_beta_0_05.png')


%% beta 0.2

beta = 0.2;

[TTC_avg2] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg2(i,:));
    end
    
    set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace')
    xlabel("n")
    ylabel("CPU time")
    title("beta 0.2")

hold off

saveas(gcf,'./Plots/ComparisonSolvers_beta_0_2.png')
%% beta 0.5

beta = 0.5;

[TTC_avg3] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg4(i,:));
    end
    
    set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace')
    xlabel("n")
    ylabel("CPU time")
    title("beta 0.5")

hold off

saveas(gcf,'./Plots/ComparisonSolvers_beta_0_5.png')
%% beta 0.8

beta = 0.8;

[TTC_avg5] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg5(i,:));
    end
    
    set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace')
    xlabel("n")
    ylabel("CPU time")
    title("beta 0.8")

hold off

saveas(gcf,'./Plots/ComparisonSolvers_beta_0_8.png')
%% beta 0.95
beta = 0.95;

[TTC_avg] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density);

figure,
hold on

    for i=1:6
        plot(problem_sizes, TTC_avg(i,:));
    end
    
    set(gca,'yscale','log')
    legend('LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace')
    xlabel("n")
    ylabel("CPU time")
    title("beta 0.95")

hold off

saveas(gcf,'./Plots/ComparisonSolvers_beta_0_95.png')

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