
function x_final = ActiveSetMethodQP(H, g, A, b, x0)

    % Get number of constraints
    [x_n, n_constraints] = size(A);

    % Set max iterations
    max_iter = 20;

    % We start by checking the active constraints
    active_constraints = find(~(A' * x0 - b));
    n_active = length(active_constraints);

    % How many constraints to choose
    W0_size = randi([0,n_active-1],1,1) + 1;

    % Random constraints
    W0 = active_constraints(randsample(n_active, W0_size));

    % Set working set equal to initial
    W_k = W0;
    x_k = x0;
    p_k = x0;
    x_final = x0;

    for i=1:max_iter

        fprintf("Iteration: %0.0f\n", i);
        fprintf("W[%0.0f]\n", i);
        disp(W_k);
        fprintf("x[%0.0f]\n", i);
        disp(x_k);
        fprintf("p[%0.0f]\n", i);
        disp(p_k);

        % Then we subset the matrices
        A_t = A';
        A_k = A_t(W_k, :);
        b_k = b(W_k);
        g_k = H * x_k + g;
        
        % Then we create the corresponding KKT-matrix and solve the system
        % (16.39)
        [LHS, RHS] = KKT_matrix_modified(x0, H, g_k, A_k, b_k);
        
        p_k = LHS \ RHS;
        p_k = p_k(1:x_n);

        if all(p_k == 0)

            % Then we find the multipliers (16.42)
            % (16.39)
            lambda_k = A_k \ g_k;

            if all(lambda_k >= 0)
                x_final = x_k;
                break;
            else
                W_k = W_k(lambda_k ~= min(lambda_k));
                x_k = x_k;
            end

        else

            % Start by selecting those that are not in the set
            anti_W_k_indices = setdiff(1:n_constraints, W_k);
            anti_A_k = A(:,anti_W_k_indices);
            anti_b_k = b(anti_W_k_indices,:);

            % Then we need to check that a_i^T * p_k < 0
            cond2 = anti_A_k' * p_k < 0;
            disp(anti_A_k);
            disp(cond2);
            anti_A_k = A(:,cond2');
            disp(anti_b_k);
            anti_b_k = b(cond2,:);
            NOMINATOR = anti_b_k - anti_A_k' * x_k;
            %disp(NOMINATOR);
            DENOMINATOR = anti_A_k' * p_k;
            %disp(DENOMINATOR);
            alpha_k = min(1, min(NOMINATOR ./ DENOMINATOR));
            x_k = x_k + alpha_k*p_k;
            disp(alpha_k)

            active_constraints_not_in_W_k = find(~(anti_A_k' * x_k - anti_b_k));
            
            if any(active_constraints_not_in_W_k)

                W_k = union(W_k, active_constraints_not_in_W_k(randsample(1, length(active_constraints_not_in_W_k))));

            else

                W_k = W_k;

            end

        end

    end





end