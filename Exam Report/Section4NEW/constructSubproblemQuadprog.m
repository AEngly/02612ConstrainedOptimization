function [A,b] = constructSubproblemQuadprog(x, con, cl, cu)
        [c,~,dc,~,~,~] = con(x);
        A = -transpose([dc dc]);
        b = [cu - c; c - cl];
end