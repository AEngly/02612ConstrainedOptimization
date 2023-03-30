function [LHS, RHS, size_x] = KKT_matrix(H, g, A, b)

    [n_row, n_col] = size(A);

    row1 = [H -A];
    row2 = [-A' zeros(n_col, n_col)];

    LHS = [row1; row2];
    RHS = [-g; b];

    if nargout > 2
        size_x = n_row;
    end

end