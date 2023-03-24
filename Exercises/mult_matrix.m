
function [X] = mult_matrix(A,B)

    dim_A = size(A);
    rows_A = dim_A(1);
    columns_A = dim_A(2);

    dim_B = size(B);
    rows_B = dim_B(1);
    columns_B = dim_B(2);

    if not(columns_A == rows_B)
        fprintf('A has dimensions (%f, %f) and B (%f %f). A must have the same number of columns as B have rows.', ...
            rows_A, columns_A, rows_B, columns_B)
        X = zeros(rows_A, columns_B);
        return
    end

    X = zeros(rows_A, columns_B);

    for i = 1:rows_A
        for j = 1:columns_B
            for k = 1:columns_A
                    X(i,j) = X(i,j) + A(i,k) * B(k, j);
            end
        end
    end

end