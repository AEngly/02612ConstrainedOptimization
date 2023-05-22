function [H,g,A,b,C,dl,du,l,u] = MinimalQP()

H = eye(2);
g = [-2;-5];
A = [1;-1];
b = zeros(1,1);
C = [ 1 1; -2 2];
dl = [-2; -1];
du = [2; 2];
l = zeros(2,0);
u = zeros(2,0);

end