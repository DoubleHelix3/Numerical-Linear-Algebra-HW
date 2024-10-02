dt = 0.00001;
T = 20;

k = 0;
x = solve(k, dt, T);

t = [0:dt:T];
x1 = x(1,:);

plot(t,x1)

function x = solve(k, dt, T)
    A = [0 1; -16 -k];
    b = [0; -12];
    x0 = [0;0];

    x = ode(A, b, x0, dt, T);
end

function x = ode(A,b,x0,dt,T)
    n = ceil(T/dt) + 1;
    % x values at various times in interval from 0 to T
    x = zeros(2,n);
    x(:,1) = x0;
    for k=2:n
        dx = (A*x(:,k-1) - b)*dt;
        x(:,k) = x(:,k-1) + dx;
    end
end
