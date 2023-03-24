function [x, lambda] = RangeSpace(H, g, A, b)

    L = chol(H);
    HG = L \ (L' \ g);
    HA = L \ (L' \ A);
    lambda = (A' * HA) \ (b + A'*HG);
    x = HA*lambda - HG;
    
end