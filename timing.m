trials = 100;
multLoopTimes = zeros(trials,1);
multTimes = zeros(trials,1);

% Timing multLoop and *
for trial = 1:trials
    n = 1000;
    A = randn(n,n);
    x = randn(n);

    startTime = cputime;
    multLoop(A,x);
    multLoopTimes(trial) = cputime - startTime;

    startTime = cputime;
    A * x;
    multTimes(trial) = cputime - startTime;
end

multAverage = mean(multTimes)
multLoopAverage = mean(multLoopTimes)

function b = multLoop(A,x)
    [m,n] = size(A);
    % Ensure correct dimensions
    assert(size(x,1) == n);

    b = zeros(1,m);
    for j = 1:n
        for i = 1:m
            b(i) = b (i) + A(i,j)*x(j);
        end
    end
end
