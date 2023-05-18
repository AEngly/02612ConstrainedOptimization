% ---------------- DESCRIPTION --------------
%
% Name: driverProblem3.m
% Type: Driver for all testing and plots from problem 3 (linear programming)
%
% Assumptions: 
% 
% 1) Equality constraint matrix A has full column rank.
%
% Problem structure:
%           min     g'x
%            x
%           s.t.    A'x + b = 0
%                   dl  <= C'x <= du
%                   l  <= x <= u
%
% Created: 12.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------

%% Global setup

% This contains e.g. options for quadprog or other solvers

%% 3.3.1) Test of Simplex (size dependent)

% This section tests whether the test problems are generated correctly.
% The function for generating the actual problems are found in separate
% files named: problemGenerator.m

parameters = struct();
parameters.n = 2;
parameters.beta = 7;
parameters.density = 0.15;
parameters.sparse = false;
[f,Aeq,beq,A,clb,cub,lb,ub,solution] = problemGenerator("RandomLP", parameters);

% Run our implementation
options = struct();
options.maxIterations = 100;
options.verbose = 1;
options.initialBasis = [];
[x1,fval1,exitflag1,output1,lambda1,Abar,bbar,gbar] = simplexCore2(f,Aeq,beq,A,cub,clb,lb,ub,options);

% Run linprog
optionsLinprog = struct();
optionsLinprog.Display = 'off';
Aineq = [A; -A];
bineq = [cub; -clb];
[x2,fval2,exitflag2,output2,lambda2] = linprog(f,Aineq,bineq,[],[],lb,ub,optionsLinprog);

% Prepare software test
tests = 0;
totalTests = 5;

fprintf("\nComparison of solutions:\n");
if norm(x1-x2,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same primal variables.\n");
    tests = tests + 1;
end
if norm(fval1-fval2,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same objective value.\n");
    tests = tests + 1;
end
if norm(lambda1.lower-lambda2.lower,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for lower bounds.\n");
    tests = tests + 1;
end
if norm(lambda1.upper-lambda2.upper,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for upper bounds.\n");
    tests = tests + 1;
end
if norm(exitflag1-exitflag2,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same exit status.\n");
    tests = tests + 1;
end

fprintf("\nStatus on tests:\n");
fprintf("Our solver passes %d/%d tests.\n", tests, totalTests);

%% 3.3.2) Test of Simplex (page 9 in Linear Optimization Simplex from lecture 7)

% This section tests whether the test problems are generated correctly.
% The function for generating the actual problems are found in separate
% files named: problemGenerator.m

parameters = struct();
[f,Aeq,beq,A,clb,cub,lb,ub,solution] = problemGenerator("Slides (page 9/27) Linear Optimization Simplex", parameters);

% Run our implementation
options = struct();
options.maxIterations = 100;
options.verbose = 0;
options.initialBasis = [];
[x1,fval1,exitflag1,output1,lambda1,Abar,bbar,gbar] = simplexCore2(f,Aeq,beq,A,cub,clb,lb,ub,options);

% We need some restructuring for linprog
optionsLinprog = struct();
optionsLinprog.Display = 'off';
Aineq = [A; -A];
bineq = [cub; -clb];
%[x2,fval2,exitflag2,output2,lambda2] = linprog(gbar,[],[],Abar,bbar,zeros(1,size(Abar,2))',[], optionsLinprog);
[x2,fval2,exitflag2,output2,lambda2] = linprog(f,Aineq,bineq,Aeq,-beq,lb,ub,optionsLinprog);

% Prepare software test
tests = 0;
totalTests = 4;

fprintf("\nComparison of solutions:\n");
if norm(x1-x2,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same primal variables.\n");
    tests = tests + 1;
end
if norm(fval1-fval2,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same objective value.\n");
    tests = tests + 1;
end
if norm(lambda1.lower-lambda2.lower,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for lower bounds.\n");
    tests = tests + 1;
end
if norm(lambda1.upper-lambda2.upper,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for upper bounds.\n");
    tests = tests + 1;
end

fprintf("\nStatus on tests:\n");
fprintf("Our solver passes %d/%d tests.\n", tests, totalTests);

%% 3.3.3) Test of Simplex (Example 13.1 from page 371 in NW)

% This section tests whether the test problems are generated correctly.
% The function for generating the actual problems are found in separate
% files named: problemGenerator.m

% Construct program
parameters = struct();
[f,Aeq,beq,A,clb,cub,lb,ub,solution] = problemGenerator("Example 13.1", parameters);

% Solve with our implementation of simplex
options = struct();
options.initialBasis = [];
options.maxIterations = 100;
options.verbose = 0;
[x1,fval1,exitflag1,output1,lambda1,Abar,bbar,gbar] = simplexCore2(f,Aeq,-beq,A,cub,clb,lb,ub,options);

% Solve with linprog
optionsLinprog = struct();
optionsLinprog.Display = 'off';
Aineq = [A; -A];
bineq = [cub; -clb];
[x2,fval2,exitflag2,output2,lambda2] = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,optionsLinprog);

% Prepare software test
tests = 0;
totalTests = 6;

fprintf("\nComparison of solutions:\n");
if norm(x1-x2,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same primal variables.\n");
    tests = tests + 1;
end
if norm(fval1-fval2,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same objective value.\n");
    tests = tests + 1;
end
if norm(lambda1.lower-lambda2.lower,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for lower bounds.\n");
    tests = tests + 1;
end
if norm(lambda1.upper-lambda2.upper,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for upper bounds.\n");
    tests = tests + 1;
end
if norm(lambda1.eqlin-lambda2.eqlin,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for equality constraints.\n");
    tests = tests + 1;
end
if norm(lambda1.ineqlin-lambda2.ineqlin,2) < (1e-8) 
    fprintf("Our implementation 'simplexCore.m' and 'linprog' reaches the same duals for inequality constraints.\n");
    tests = tests + 1;
end
fprintf("\nStatus on tests:\n");
fprintf("Our solver passes %d/%d tests.\n", tests, totalTests);

%% Test OF TOTAL SIMPLEX


% This section tests whether the test problems are generated correctly.
% The function for generating the actual problems are found in separate
% files named: problemGenerator.m

parameters = struct();
[f,Aeq,beq,A,clb,cub,lb,ub,solution] = problemGenerator("Slides (page 9/27) Linear Optimization Simplex", parameters);
beq = -beq;

% Run our implementation
options = struct();
options.maxIterations = 100;
options.verbose = 0;
options.initialBasis = [];
[x1,fval1,exitflag1,output1,lambda1,Abar,bbar,gbar] = simplex(f,Aeq,beq,A,cub,clb,lb,ub,options);


%% 3.4) Adjusting QP algorithms to solve LP

% This section contains tests for adjusted general QP solvers.
% The adjusted implementation are found is separate files named:
% QP_InteriorPointPDPC_adjusted.m



%% 3.6) Testing of primal-dual interior point tailored for general LP

% This section contains tests for primal-dual interior point algorithm for LP.
% The implementation is found in the file named: XXXX

% The test should be on some general LP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics of:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

n = 5;
meq = 10;
miq = 10;

totalN = 2*meq + 2*miq + 2*n;

state = 1000;
%rand('state',state);

g = rand(n,1);
Aeq = randn(meq,n)*100;
beq = rand(meq,1);
A = randn(miq,n)*100;
clb = -10*rand(miq,1);
cub = 10*rand(miq,1);
lb = -5*rand(n,1);
ub = 5*rand(n,1);

% Then we construct the system matrix
[Abar, bbar, gbar, information] = standardForm(g,Aeq,beq,A,cub,clb,lb,ub);
n_new = size(Abar,2);
m_new = size(Abar,1);

% Make duals for inequality
lambda = rand(m_new,1);

% Make x (primal) and mu (dual for nonnegativity) such that complementarity
% holds

x = zeros(n_new,1);
x(1:m_new,1) = abs(rand(m_new,1));
mu = zeros(n_new,1);
mu(m_new+1:n_new,1) = abs(rand(n_new-m_new,1));

% Then reverse engineer bbar (bbar = [-beq; -cub; clb; -ub; lb];)
bbar_new = Abar*x;
beq = -bbar_new(1:meq);
cub = -bbar_new(1+meq:meq+miq);
clb = bbar_new(1+meq+miq:meq+2*miq);
ub = -bbar_new(1+meq+2*miq:meq+2*miq+n);
lb = bbar_new(1+meq+2*miq+n:meq+2*miq+2*n);
g = Abar'*lambda + mu;

% Then try to solve the system
[xlp,info,mulp,lambdalp,iter] = LP_InteriorPointPrimalDual(g,sparse(Abar),bbar_new,ones(n_new,1));

if info
    X = max(abs(xlp-x))
    M = max(abs(mulp-mu))
    L = max(abs(lambdalp-lambda))
end

%% 3.8) Testing of primal active-set for general LP (Simplex)

% This section contains tests for primal active-set for general LP (Simplex)
% The implementation is found in the file named: XXXX

% The test should be on some general LP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics from scalable problems:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

%% 3.9) Compare implementation with linprog, MOSEK, Gurobi, and cvx.

% The test should be on some general LP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics from scalable problems:
% 1) Size of problem (n) on x-axis, and time-to-completion on y-axis.

