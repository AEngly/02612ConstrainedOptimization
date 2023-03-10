
function x = LUSolver(LHS, RHS, type)

    if type == "dense"
        [L, U] = lu(LHS);
        y = L \ RHS; % Forward substitution
        x = U \ y;   % Backward substitution

    elseif type == "sparse"
        LHS = sparse(LHS);
        RHS = sparse(RHS);
        [L, U] = lu(LHS);
        y = L \ RHS; % Forward substitution
        x = U \ y;   % Backward substitution
    end

end