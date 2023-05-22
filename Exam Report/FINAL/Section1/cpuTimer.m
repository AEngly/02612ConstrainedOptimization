function time = cpuTimer(f, varargin)

    start = cputime;
    output = f(varargin{:});
    time = cputime - start;

end