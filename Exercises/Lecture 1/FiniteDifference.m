
function [fVal, fGradient, fHessian] = FiniteDifference(f, x, h)

    fVal = f(x);

    fGradient = [(f([x(1) + h; x(2)]) - f(x))/h; (f([x(1); x(2) + h]) - f(x))/h];

    xx = (f([x(1) + h; x(2)]) - 2*f(x) + f([x(1) - h; x(2)]))/h.^2;
    yy = (f([x(1); x(2) + h]) - 2*f(x) + f([x(1); x(2) - h]))/h.^2;
    xy = (f([x(1) + h; x(2) + h]) - f([x(1) - h; x(2) + h]) - f([x(1) + h; x(2) - h]) + f([x(1) - h; x(2) - h]))/(4.*h.^2);

    fHessian = [xx xy; xy yy];

end