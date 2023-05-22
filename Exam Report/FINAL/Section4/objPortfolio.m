function [f, Grad] = objPortfolio(x,mu,D,gamma,n,k)

    y = x(n+1:n+k);
    x = x(1:n);
    
    f = transpose(x)*D*x + transpose(y)*y + (gamma^(-1))*transpose(mu)*x;

    if nargout > 1
        Grad = 0;
    end

end