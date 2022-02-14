% generates random LP data
% with (x,s) = (e,e) feasible
n = 1000; % number of variables
m = 500; % number of constraints
A = randn(m,n);
b = sum(A,2);
c = ones(n,1) - A'*randn(m,1);
ga = 1000;

% solve LP with long-step method along fixed combination of centering and
% affine scaling direction
tic
NT_iter_count = 0;
T = zeros(1,100);
M = zeros(1,100);
x = ones(n,1);
s = ones(n,1);
mu = 1;
while mu > 10^(-8)
    [dx,ds,px,ps] = NTdirections(A,b,c,x,s);
    tmax = 2;
    tmin = 0;
    while tmax - tmin > 10^(-4)
        t = (tmax + tmin)/2;
        xn = x + t*(0.9*px + 0.1*dx);
        sn = s + t*(0.9*ps + 0.1*ds);
        y = xn.*sn;
        if (min(xn) > 0) && (min(sn) > 0) && (max(y)/min(y) < ga)
            tmin = t;
        else
            tmax = t;
        end
    end
    x = x + tmin*(0.9*px + 0.1*dx);
    s = s + tmin*(0.9*ps + 0.1*ds);
    y = x.*s;
    mu = sum(y)/n;
    NT_iter_count = NT_iter_count + 1;
    T(NT_iter_count) = tmin;
    M(NT_iter_count) = log(mu);
end
T = T(1:NT_iter_count);
M = M(1:NT_iter_count);
NT_optval = c'*x;
NT_time = toc;

tic
cvx_begin
    variable x(n) nonnegative
    minimize(c'*x)
    A*x == b
cvx_end
cvx_time = toc;

format long
NT_iter_count
NT_optval
cvx_optval
NT_time
cvx_time
