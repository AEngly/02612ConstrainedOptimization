function [c,ceq,GC,GCeq] = conNLPTest(x,n)

    d = 1; % Must be >= 1
    k = floor(n/2)-d;

    % 1) Inequalities
    c = zeros(2+k,1);
    c(1,1) = sum(x) - n;
    c(2,1) = prod(x) - 1;
    for j = 3:k+3
        c(j,1) = sum((x(1:end-j) - x(1+j:end)).^2) - 1;
    end

    % 2) Compute gradients of inequalities
    GC = zeros(n,2+k);
    GC(:,1) = ones(n,1);
    for i=1:n
        GC(i,2) = prod(x(1:i-1))*1*prod(x(i+1:end));
    end
    for j=3:k+3
        GC(:,j) = [2*x(1:(j-2)); zeros(n-2*(j-2),1); -2*x((n-j+3):n)]';
    end

    % 3) Set empty equalities
    ceq = [];
    GCeq = [];

end