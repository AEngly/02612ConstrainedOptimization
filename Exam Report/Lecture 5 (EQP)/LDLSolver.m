
function x = LDLSolver(LHS, RHS, type)

    if type == "dense"
        [L,D,~] = ldl(LHS);
        z = L \ RHS; % Forward substitution
        y = D \ z; % Forward substitution
        x = L' \ y;   % Backward substitution

    elseif type == "sparse"
        LHS = sparse(LHS);
        [L,D,p] = ldl(LHS, 'lower','vector');
        z = L \ RHS(p); % Forward substitution
        y = D \ z; % Forward substitution
        x(p) = L' \ y;   % Backward substitution
        x = x';
    end

end