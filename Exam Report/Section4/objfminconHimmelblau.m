function [f, Grad] = objfminconHimmelblau(x)

    x1 = x(1);
    x2 = x(2);

    tmp1 = x1.*x1+x2-11;
    tmp2 = x1+x2.*x2-7;

    f = -tmp1.*tmp1 + tmp2.*tmp2;

    if nargout > 1
        Grad = -[2.*x1 + 4.*x1.*(x1.^2 + x2 - 11) + 2.*x2.^2 - 14;
                2.*x2 + 4.*x2.*(x2.^2 + x1 - 7) + 2.*x1.^2 - 22];
    end

end