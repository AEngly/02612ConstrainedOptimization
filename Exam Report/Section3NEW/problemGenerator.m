function [f,A,b,C,dl,du,l,u,solution] = problemGenerator(problemName,parameters)

% ---------------- DESCRIPTION --------------
%
% Name: Generator for Random QP
% Type: Generates random programs
%
% Problem structure:
%           min     g' x
%            x
%           s.t.    A'x + b = 0
%                   dl <= C'x <= du
%                   l <= x <= u
%
% Syntax: [H,g,A,b,C,dl,du,l,u] = GeneratorQP(n,alpha,beta,density)
%
%
% Created: 30.03.2023
% Authors: Andreas Engly (s170303) and Karl Takeuchi-Storm (s130377)
%          Compute, Technical University of Denmark
%
% ---------------- IMPLEMENTATION --------------

    % Set large number for missing bounds
    M = 100;

    if problemName == "RandomLP"
        
        n = parameters.n;
        beta = parameters.beta;
        density = parameters.density;
        sparse = parameters.sparse;

        % Constraints as percentage of variables
        m = round(beta*n);

        % Create empty A and b, as no equality constraints are considered
        b = zeros(0,1);
        A = zeros(0,n);
        C = sprandn(m,n,density);
        if ~sparse
            C = full(C);
        end
        dl = -3*rand(m,1);
        du = 3*rand(m,1);

        % Create empty upper and lower bounds for x, as x is unconstrained
        l = -3*rand(n,1);
        u = 3*rand(n,1);
        f = randn(n,1);
        solution = struct();
        solution.primal = [];

    elseif problemName == "Example 13.1"

        % Example from page 43 in 'Introduction to OR'
        M = 1000;
        f = [-3; -2; 0; 0];
        A = [1 1 1 0;
             2 0.5 0 1];
        b = [5;8];
        C = [];
        dl = [];
        du = [];
        l = [0;0;0;0];
        u = [];
        solution = struct();
        solution.primal = [0.5000; 0.2000; 0.3000];

    elseif problemName == "Slides (page 9/27) Linear Optimization Simplex"

        % Example from page 43 in 'Introduction to OR'
        f = [-3; 2; -4];
        A = [1 1 1];
        b = [-1];
        C = [1 1 0;
             1 0 1];
        dl = [0.7;-M];
        du = [M;0.8];
        l = [0; 0; 0];
        u = [1; 1; 1];
        solution = struct();
        solution.primal = [0.5000; 0.2000; 0.3000];

    end

end