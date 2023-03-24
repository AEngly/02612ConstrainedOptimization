
% Description: Driver for problem 1


%% 1.3)

[H, g, A, b] = RecycleSystem(3, 0.2, 1);

%% 1.4)

[LHS, RHS] = KKT_matrix(H, g, A, b);

%% 1.5)

u1_dense = LUSolver(LHS, RHS, "dense");

%% 1.6)

u2_dense = LDLSolver(LHS, RHS, "dense");

%% Null-Space

[x1, lambda1] = NullSpace(H, g, A, b);

%% Range-Space

[x2, lambda2] = RangeSpace(H, g, A, b);

%% LU sparse

u1_sparse = LUSolver(LHS, RHS, "sparse");

%% LDL sparse

start = cputime;
u3_sparse = LDLSolver(LHS, RHS, "sparse");
final = start - cputime;
display(final);

%% Comparison of solvers

n = 175;
step = 25;

TTC = zeros(7,n);
j = 1;

problem_sizes = 10:step:n*step;

for i = problem_sizes

    % Display
    fprintf('Problem size: %d\n', i);

    [H, g, A, b] = RecycleSystem(i, 0.2, 1);
    [LHS, RHS] = KKT_matrix(H, g, A, b);
    
    %TTC(1,j) = cpuTimer(@LUSolver, LHS, RHS, "dense");
    %TTC(2,j) = cpuTimer(@LDLSolver, LHS, RHS, "dense");
    TTC(3,j) = cpuTimer(@LUSolver, LHS, RHS, "sparse");
    TTC(4,j) = cpuTimer(@LDLSolver, LHS, RHS, "sparse");
    %TTC(5,j) = cpuTimer(@NullSpace, H, g, A, b);
    %TTC(6,j) = cpuTimer(@RangeSpace, H, g, A, b);
    %TTC(7,j) = cpuTimer(@quadprog, H,g,[],[],A',b);

    j = j + 1;

end

%% Plotting the comparison

hold on

for i=3:4
    plot(problem_sizes, TTC(i,:));
end

%legend('LUSolver', 'LDLSolver', 'LUSolver (sparse)', 'LDLSolver (sparse)', 'NullSpace', 'RangeSpace', 'Quadprog')
%legend('LUSolver', 'LDLSolver', 'NullSpace', 'RangeSpace')
legend('LUSolver (sparse)', 'LDLSolver (sparse)')
xlabel("n")
ylabel("CPU time")

hold off

saveas(gcf,'ComparisonSparseSolvers.png')

%% 1.10) Sparsity of KKT-system

[H, g, A, b] = RecycleSystem(3, 0.2, 1);
[LHS, RHS] = KKT_matrix(H, g, A, b);
spy(LHS);

%% 1.11) SPY

[x,fval,exitflag,output,lambda] = quadprog(H,g,[],[],A',b);

