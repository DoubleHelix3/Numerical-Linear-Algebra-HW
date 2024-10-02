t = [-1 -0.75 -0.5 0 0.25 0.5 0.75]';
y = [1 0.8125 0.75 1 1.3125 1.75 2.3125]';

poly_regression(t,y,2)

t = [1000 1050 1060 1080 1110 1130];
y1 = [6010 6153 6421 6399 6726 6701];
y2 = [9422 9300 9220 9150 9042 8800];

poly_regression(t,y1,1)
poly_regression(t,y2,1)

% polynomial regression of degree d
function c = poly_regression(t,y,d)
    assert(size(t,1) == size(y,1));

    % Ac = y
    m = size(t,1);
    A = zeros(m,d+1);
    for i=1:m
        for j=1:d+1
            A(i,j) = t(i)^(j-1);
        end
    end
    
    c = least_squares(A,y);
end

function x = least_squares(A,b)
    % A' * Ax = A' * b
    % Least squares solution is R'Rx = R' Q' b
    [Q,R] = qr(A);

    % I wrote these functions on a previous homework assignment
    Rx = forwardsub(R', R' * Q' * b);
    x  = backsub(R, Rx);
end

function x = forwardsub(A,b)
    x = A\b;
end
function x = backsub(A,b)
    x = A\b;
end
