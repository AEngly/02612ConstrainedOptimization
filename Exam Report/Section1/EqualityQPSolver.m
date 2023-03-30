function [x, lambda] = EqualityQPSolver(H, g, A, b, solver)

% ---------------- DESCRIPTION --------------
%
% Name: Range Space Solver for EQP
% Type: Solver for Equality Constrained Quadratic Programs
%
% Problem structure:
%           min     0.5 x' H x + g' x
%            x
%           s.t.    A'x + b = 0
%
% Syntax: [x, lambda] = EqualityQPSolverRS(H, g, A, b)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    if solver == "LUdense"
        [x, lambda] = EqualityQPSolverLUdense(H, g, A, b);
    elseif solver == "LUsparse"
        [x, lambda] = EqualityQPSolverLUsparse(H, g, A, b);
    elseif solver == "LDLdense"
        [x, lambda] = EqualityQPSolverLDLdense(H, g, A, b);
    elseif solver == "LDLsparse"
        [x, lambda] = EqualityQPSolverLDLsparse(H, g, A, b);
    elseif solver == "RangeSpace"
        [x, lambda] = EqualityQPSolverRS(H, g, A, b);
    elseif solver == "NullSpace"
        [x, lambda] = EqualityQPSolverNS(H, g, A, b);
    else
        [x, lambda] = nan;
    end
    
end