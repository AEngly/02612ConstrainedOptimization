
% Description:
% Authors: 
% Last Edit: 03/03-2023

function [H, g, A, b] = RecycleSystem(n, u, d0)

    H = eye(n+1);
    A = zeros(n,n+1);

    A(1,1) = -1;
    A(1,n) = 1;

    for j = 1:(n-2)
        i = j+1;
        A(i,i-1) = 1;
        A(i,i) = -1;
    end

    A(n, n-1) = 1;
    A(n, n) = -1;
    A(n, n+1) = -1;
    A = A';

    g = -u * ones(n+1,1);

    b = zeros(n,1);
    b(1) = -d0;

end