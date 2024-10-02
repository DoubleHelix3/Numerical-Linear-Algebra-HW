% 
% ns = [4,8,12,16];
% for i=1:length(ns)
%     n = ns(i);
%     z = ones(n,1);
% 
%     A = hilb(n);
%     b = A*z;
%     % the theoretical solution of Ax = b is extremely close to z because
%     % matrix multiplication is just additions and multiplications
%     % no catastrophic cancelation because all terms are positive
% 
%     xhat = A\b; % approximate solution
% 
%     n
%     condition_number = cond(A);
%     absolute_error = norm(xhat - z);
%     residual_error = norm(b - A*xhat);
% end


% 2.7.25

n = 160;
A = randn(n);
A(1,1) = 50*eps*A(1,1);

[L_hat,U_hat] = no_partial_pivoting(A);

L_hat_norm = norm(L_hat,inf)
U_hat_norm = norm(U_hat,inf)

E_hat = L_hat*U_hat - A;
E_hat_norm = norm(E_hat,inf)

[L,U] = lu(A);
L_norm = norm(L, inf)
U_norm = norm(U, inf)
E = L*U - A;
E_norm = norm(E,inf)

% Gaussean Elimination without pivoting
function [L,U] = no_partial_pivoting(A)
    n = size(A,1);
    assert(n == size(A,1)); % ensure A is square

    L = zeros(n);
    U = A;

    for k=1:(n-1)
        for i=(k+1):n
            l = U(i,k)/U(k,k);
            L(i,k) = l;
            U(i,:) = U(i,:) - l*U(k,:);
        end
    end

    L = L + eye(n);
end

