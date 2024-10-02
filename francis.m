rng("default");

n = 10;
md = randn(1,n); % random main diagonal
sd = randn(1,n-1); % random subdiagonal
a = diag(md) + diag(sd,-1) + diag(sd,1);

iterations = 0;
error = inf;
preverror = 1;
alpha = 1
while error > 1e-16
    a = francisIter(a, 4.964);

    error = abs(a(n,n-1));
    S = error/(preverror^alpha)
    preverror = error;
    iterations = iterations + 1;
end

iterations

% computes francis iteration given a shift
function a = francisIter(a, shift)
    [n,n] = size(a);

    % Create the bulge
    cs = a(1,1) - shift; sn = a(2,1);
    r = norm([cs sn]);
    cs = cs/r; sn = sn/r;
    % normalizing to make csA2 + snA2 = 1
    q0 = [ cs -sn; sn cs]; % Givens rotator
    a(1:2,:) = q0'*a(1:2,:); % left multiplication
    a(:,1:2) = a(:,1:2)*q0; % right multiplication

    % Chase the bulge
    for ii = 1:n-2
        % Chase the bulge from position (ii+2,ii).
        cs = a(ii+1,ii); sn = a(ii+2,ii);
        r = norm([cs sn]);
        cs = cs/r; sn = sn/r;
        a(ii+1,ii) = r; a(ii+2,ii) = 0;
        qi = [cs -sn; sn cs];
        % Givens rotator to chase the bulge,
        a(ii+1:ii+2, ii+1:n) = qi'*a(ii+1 : ii+2, ii+1:n);
        a (:, ii+1 : ii+2) = a(:, ii+1 : ii+2)*qi;
    end
end

% computes wilkinson shift of matrix a
function shift = wilkinson(a)
    [n,n] = size(a);

    htr = 0.5 * (a (n-1, n-1) + a(n,n)); % half a 2x2 trace
    dscr = sqrt((.5*(a(n-1,n-1)-a(n,n)))^2 + a(n,n-1)^2);
    % discriminant
    if htr < 0, dscr = -dscr; end % to avoid cancellation
        root1 = htr + dscr; % quadratic formula
    if root1 == 0 % almost never happens
        root2 = 0;
    else % almost always happens
        det = a (n-1,n-1)*a (n,n) - a(n,n-1)^2;
        % 2x2 determinant = product of roots
        root2 = det/root1;
    end
    if abs(a(n,n)-root1) < abs(a(n,n)-root2)
        shift = root1;
    else
        shift = root2;
    end
end
