% DriverLP   Computes solutions and plots for exercise 3 (Linear
% Programming)

% Created: 04.12.2007
% Author : Andreas Heidelbach Engly
%          IMM, Technical University of Denmark

%
% Test LP solver
%

state = 1000;
rand('seed', state);

% STEP 1: Construct problem
test_range = 10:5:200; % number of variables
m = 10; % number of constraints
upper_bound = 1;
lower_bound = 0;
counter = 1;

% Information to save
methods = 3;
iterations = zeros(methods, length(test_range));

% Options for linprog
options1 = optimoptions('linprog','Algorithm','interior-point', 'Display', 'off');
options2 = optimoptions('linprog','Algorithm','dual-simplex', 'Display', 'off');

for n = test_range
    
    A = randn(m,n);
    x = zeros(n,1);
    x(1:m,1) = abs(rand(m,1));
    lambda = zeros(n,1);
    lambda(m+1:n,1) = abs(rand(n-m,1));
    mu = rand(m,1);
    lb = ones(n,1)*lower_bound;
    ub = ones(n,1)*upper_bound;
    
    g = A'*mu + lambda;
    b = A*x;
    
    [xlp,info,mulp,lambdalp,iter] = LPippd(g,A,b,ones(n,1));
    [x1,fval1,exitflag1,output1,lambda1] = linprog(g,[],[],A,b,lb,ub,options1);
    [x2,fval2,exitflag2,output2,lambda2] = linprog(g,[],[],A,b,lb,ub,options2);

    iterations(1, counter) = iter;
    iterations(2, counter) = output1.iterations;
    iterations(3, counter) = output2.iterations;

    counter = counter + 1;

end

hold on

for i=1:methods
    plot(test_range, iterations(i,:));
end

legend('Own solver', 'Linprog (interior-point)', 'Linprog (dual-simplex)');
xlabel("n")
ylabel("Iterations")

hold off

saveas(gcf,'ComparisonSolversLP.png')
