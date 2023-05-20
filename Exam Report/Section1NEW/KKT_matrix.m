function [LHS, RHS, n_row, n_col] = KKT_matrix(H, g, A, b)

    [n_row, n_col] = size(A);
    row1 = [H -A];
    row2 = [-A' zeros(n_col, n_col)];
    LHS = [row1; row2];
    RHS = [-g; b];

end