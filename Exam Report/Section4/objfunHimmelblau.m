function [f, Grad] = objfunHimmelblau(x, p)

    tmp1 = x(1)*x(1)+x(2)-11;
    tmp2 = x(1)+x(2)*x(2)-7;

    f = tmp1*tmp1 + tmp2*tmp2;

    if nargout > 1
        Grad = [2*x(1) + 4*x(1)*(x(1)^2 + x(2) - 11) + 2*x(2)^2 - 14;
                2*x(2) + 4*x(2)*(x(2)^2 + x(1) - 7) + 2*x(1)^2 - 22];
    end

end