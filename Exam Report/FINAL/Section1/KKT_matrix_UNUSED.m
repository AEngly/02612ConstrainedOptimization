function [LHS, RHS] = blabla(x, H, g, A, b)

    [n_row, ~] = size(A);

    row1 = [H A'];
    row2 = [A zeros(n_row, n_row)];

    LHS = [row1; row2];
    RHS = [H*x + g; A*x - b];

end