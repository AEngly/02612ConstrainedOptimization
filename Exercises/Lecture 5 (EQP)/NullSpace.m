function [x, lambda] = NullSpace(H, g, A, b)

    [Q, R] = qr(A);
    [n, m] = size(A);
    Q1 = Q(:,1:m);
    Q2 = Q(:,m+1:n);
    R = R(1:m,1:m);
    xy = R' \ b;
    xz = Q2'*H*Q2 \ (-Q2'*(H * Q1 * xy + g));
    x = Q1 * xy + Q2 * xz;
    lambda = R \ Q1'*(H*x + g);

end