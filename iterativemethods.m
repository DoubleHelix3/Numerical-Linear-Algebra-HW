w = 1;
method = @(A,b,x) sor(A,b,w,x); % replace with whatever algorithm
h = 1/4;
m = 1/h; % number of points in discritization is m*m

f = zeros(m-1,m-1);
gL = zeros(1,m);
gR = zeros(1,m);
gU = zeros(1,m);
gD = zeros(1,m);
u0 = zeros(1,m-1);
solve_poisson(f,g,m,method);


function [u,iterations] = solve_poisson(f,gL,gR,gU,gD, u0, method)
    m = size(f,1)+1; % f is m-1 by m-1
    h = 1/m; % step size

    % Construct matrix equation from (m-1)*(m-1) scalar equations
    A = zeros(m-1);
    b = zeros(m-1,1);
    pack = @(i,j) j*(m-1)+i;
    I = eye(m-1);
    T = diag(4*ones(m-1,1)) + diag(-1*ones(m-2,1),-1) + diag(-1*ones(m-2,1),1);
    C = {};
    for i=1:m-1
        for j=1:m-1
            if i == j
                C(i,j) = {T};
            elseif i==j+1 | i==j-1
               C(i,j) = {-I};
            else
               C(i,j) = {zeros(m-1)};
            end
        end
     end
     A = cell2mat(C);

    method(A,b,u0);
end

% Utility function to decompose A = L + D + U
function [L,D,U] = sum_decomp(A)
    L = tril(A,-1);
    D = diag(diag(A)); % lmao
    U = triu(A,1);
end

% One step of Jacobi iteration 
function x_new = jacobi(A,b,x)
    [L,D,U] = sum_decomp(A);
    M = D;
    N = -(L + U);
    x_new = mn_iter(M,N,b,x);
end

% One step of Gauss-Seidel iteration 
function x_new = gs(A,b,x)
    [L,D,U] = sum_decomp(A);
    M = D + L;
    N = -U;
    x_new = mn_iter(M,N,b,x);
end

% One step of SOR iteration 
function x_new = sor(A,b,w,x)
    [L,D,U] = sum_decomp(A);
    M = D + w*L;
    N = ((1-w)*D - w*U);
    x_new = mn_iter(M,N,w*b,x); % note: b changed to w*b
end

% One step of MN iteration
function x_new = mn_iter(M,N,b,x)
    % Solve M*x_new = N*x + b
    x_new = M\(N*x + b);
end
