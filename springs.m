n = 20;
k = ones(n+1,1)
S = springMatrix(k);

f = zeros(n,1);
f(5) = 1;
f(16) = -1;

x = solveSprings(S,f)

function S = springMatrix(k)
    n = size(k,1) - 1;
    % Creating an nxn tridiagonal matrix 
    % -(k(i) + k(i+1)) along the main
    % k(i+1) along the sub and super diagonals
    S = diag(-(k(1:n) + k(2:n+1))) + diag(k(2:n),1) + diag(k(2:n),-1);
end

function x = solveSprings(S, f)
    % Sx = -f
    x = S\(-f);
end
