% sparse recovery
% dimensional parameters
m = 100; % number of observations
n = 500; % dimension of decision vector
p = 0.05; % sparsity factor
s = 0.001; % noise level
% generate data
A = randn(m,n);
x = zeros(n,1);
fprintf("Nonzero components of true vector:\n");
count = 0;
for k = 1:n
    if rand < p
        x(k) = randn;
        count = count + 1;
        fprintf("x_%d = %5e\n",k,x(k));
    end
end
fprintf("%d non-zero components\n",count);
%y = A*x + s*randn(m,1);
y = A*x + s*(2*rand(m,1) - 1);
% solve problem
d = 0.001;
cvx_begin quiet
    variable z(n,1)
    variable t(n,1)
    minimize( ones(1,n)*t )
    A*z - y <= d*ones(m,1)
    -d*ones(m,1) <= A*z - y
    z <= t
    -t <= z
cvx_end
fprintf("Nonzero components of recovered vector:\n");
count = 0;
for k = 1:n
    if abs(z(k)) > 0.05
        fprintf("z_%d = %5e\n",k,z(k));
        count = count + 1;
    end
end
fprintf("%d non-zero components\n",count);
close all
plot(x)
hold on
plot(z)
ylabel('x,z')
