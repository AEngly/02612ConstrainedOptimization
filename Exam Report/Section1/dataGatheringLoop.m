function [TTC_avg] = dataGatheringLoop(beta, problem_sizes, smoother, l, alpha, density) 
    TTC = zeros(6,l,smoother);
    TTC_avg = zeros(6,l);
    j = 1;

    for i = problem_sizes
    
        % Display
        %fprintf('Problem size: %d\n', i);
    
       
    
        for k = 1:smoother

            [H,g,A,b] = GeneratorECQP(i,alpha,beta,density);
    
            TTC(1,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUdense");
            TTC(2,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LUsparse");
            TTC(3,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLdense");
            TTC(4,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "LDLsparse");
            TTC(5,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "RangeSpace");
            TTC(6,j,k) = cpuTimer(@EqualityQPSolver, H, g, A, b, "NullSpace");
        end
    
        TTC_avg(1,j) = mean(TTC(1,j,:));
        TTC_avg(2,j) = mean(TTC(2,j,:));
        TTC_avg(3,j) = mean(TTC(3,j,:));
        TTC_avg(4,j) = mean(TTC(4,j,:));
        TTC_avg(5,j) = mean(TTC(5,j,:));
        TTC_avg(6,j) = mean(TTC(6,j,:));
    
        j = j + 1;
    
    end
end