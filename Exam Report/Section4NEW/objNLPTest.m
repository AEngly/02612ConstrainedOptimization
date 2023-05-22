function [f, Grad] = objNLPTest(x,n)

    % Supply a vector x and the size of the problem n
    sizeVector = 1:n;
    f = sum((x.^2 - sizeVector').^2 + sin(sizeVector'));

    if nargout > 1
        Grad = 6*x + cos(x);
    end

end