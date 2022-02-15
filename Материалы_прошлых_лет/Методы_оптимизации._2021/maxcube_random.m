% maximum of quadratic function over hypercube
close all
clear all

n = 30; % matrix size
N = 100; % number of generated random vertices

% generate cost matrix
W = randn(n);
Q = W*W';

% solve semi-definite relaxation
cvx_begin
    variable X(n,n) semidefinite
    diag(X) == 1
    maximize(trace(Q*X))
cvx_end
upper_bound = cvx_optval;
X = X - diag(diag(X)) + eye(n);
expectation = 2/pi*trace(Q*asin(X));

% generating cuts
average = 0;
C = zeros(n,N);
vertex_value = zeros(1,N);
count = 0;
[U,D] = eig(full(X));
D = max(D,0);
Xsq = U*sqrt(D);
for k = 1:N,
    x = Xsq*randn(n,1);
    x = x./abs(x);
    x = x*x(1);
    value = x'*Q*x;
    average = average + value;
    Dx = C(:,1:count) - x*ones(1,count);
    if (count == 0) || (min(sum(abs(Dx))) > 0),
        count = count + 1;
        C(:,count) = x;
        vertex_value(count) = value;
    end,
end,
average = average/N;
C = C(:,1:count);
vertex_value = vertex_value(1:count);
[cv,ind] = sort(vertex_value);
C = C(:,ind);
figure
hold on
title('Values of random vertices')
xlabel('Number of vertex')
ylabel('Value of vertex')
plot(1:count,upper_bound*ones(1,count),'r:')
plot(1:count,average*ones(1,count),'b:')
plot(1:count,expectation*ones(1,count),'c:')
plot(1:count,2/pi*upper_bound*ones(1,count),'g:')
plot(1:count,cv,'k.')
legend('upper bound','average','expectation','2/\pi * upper bound')
