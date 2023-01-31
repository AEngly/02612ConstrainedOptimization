
function [c, dc, d2c] = FunEx2(f, x)
    c = f(x);
    dc = gradient(f, x);
    d2c = hessian(f, x);
end