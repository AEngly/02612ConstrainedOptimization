
function [f, fGradient, fHessian] = FunEx1(x)
    f = x(1).^2 - 2*x(1) + 3*x(1).*x(2) + 4*x(2).^2;
    fGradient = [2*x(1) - 2 + 3*x(2); 3*x(1) + 8*x(2)];
    fHessian = [2 3; 3 8];
end