% ---------------- DESCRIPTION --------------
%
% Name: driverProblem3.m
% Type: Driver for all testing and plots from problem 3 (linear programming)
%
% Assumptions: 
% 
% 1) Equality constraint matrix A has full column rank.
%
% Problem structure:
%           min     g'x
%            x
%           s.t.    A'x + b = 0
%                   dl  <= C'x <= du
%                   l  <= x <= u
%
% Created: 12.05.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%

% ---------------- IMPLEMENTATION --------------

%% Global setup

% This contains e.g. options for quadprog or other solvers

%% 3.3) Construction of size dependent test problems for LP

% This section tests whether the test problems are generated correctly.
% The function for generating the actual problems are found in separate
% files named: XXXXXX

%% 3.4) Adjusting QP algorithms to solve LP

% This section contains tests for adjusted general QP solvers.
% The adjusted implementation are found is separate files named:
% XXXXXXXXXXXXXXX

%% 3.6) Testing of primal-dual interior point tailored for general LP

% This section contains tests for primal-dual interior point algorithm for LP.
% The implementation is found in the file named: XXXX

% The test should be on some general LP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics of:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

%% 3.8) Testing of primal active-set for general LP (Simplex)

% This section contains tests for primal active-set for general LP (Simplex)
% The implementation is found in the file named: XXXX

% The test should be on some general LP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics from scalable problems:
% 1) Number of iterations
% 2) Average time per iteration
% 3) Compare difference between solutions of our solver and a library
%    function.

%% 3.9) Compare implementation with linprog, MOSEK, Gurobi, and cvx.

% The test should be on some general LP for different starting points.
% The plots should show the paths.
%
% We should also provide statistics from scalable problems:
% 1) Size of problem (n) on x-axis, and time-to-completion on y-axis.

