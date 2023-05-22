function [time, output] = cpuTimer(f, varargin)

    start = cputime;
    results = f(varargin{:});
    time = cputime - start;

    if argout > 1
        output = results;
    end

end