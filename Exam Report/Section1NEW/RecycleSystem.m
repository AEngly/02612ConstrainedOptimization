
function [H, g, A, b] = RecycleSystem(n)

% ---------------- DESCRIPTION --------------
%
% Name: RecycleSystem.m
% Type: Generates recycle problem for exercises week 5
%
% Syntax: [H, g, A, b] = RecycleSystem(n)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    ubar = 0.2;
    d0 = 1;

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

    g = -ubar * ones(n+1,1);

    b = zeros(n,1);
    b(1) = -d0;

end