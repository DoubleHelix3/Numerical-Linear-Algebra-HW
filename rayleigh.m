
n = 100;
A = randn(n);
% Upper Hessenberg form 
[Q,B] = hess(A);

% Rayleigh iteration with random initial vector
q = randn(n,1) + i*randn(n,1);
q = q/norm(q);
rayquo = 0;
preverror = 1;
iterations = 1;


while iterations < 100
    rayquo = q'*B*q; % rayleigh quotient
    % check if matrix is close to singular    
    if rcond(B - rayquo*eye(n)) < 1e-14
        break
    end
    q = (B - rayquo*eye(n))\q;
    q = q/norm(q);


    % eigenvector and associated eigenvalue of A
    v = Q*q;
    lambda = rayquo;

    % To see rate of convergence
    error = norm(A*v - lambda*v);
    S = error/(preverror*preverror) % quadratic convergence
    preverror = error;

    iterations = iterations + 1;
end

iterations
error
