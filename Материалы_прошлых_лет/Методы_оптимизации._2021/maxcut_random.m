% MaxCut relaxation on random instances
close all
clear all

n = 30; % number of vertices in the graph
p = 0.5; % probability of an edge
N = 100; % number of generated random cuts
m = 4; % number of best cuts plotted

% generate graph
W = zeros(n);
for k1 = 1:n-1,
    for k2 = k1+1:n,
        if rand < p,
            W(k1,k2) = 1;
            W(k2,k1) = 1;
        end,
    end,
end,

% plot graph
figure
hold on
phi = (1:n)*2*pi/n;
plot(cos(phi),sin(phi),'k.')
for k1 = 1:n-1,
    for k2 = k1+1:n,
        if W(k1,k2),
            plot(cos(phi([k1 k2])),sin(phi([k1 k2])),'k')
        end,
    end,
end,
axis([-1.1 1.1 -1.1 1.1])

% solve semi-definite relaxation
cvx_begin
    variable X(n,n) semidefinite
    diag(X) == 1
    maximize(1/4*trace(W*(1 - X)))
cvx_end
upper_bound = cvx_optval;
expectation = 1/(2*pi)*trace(W*acos(X));
title(['Upper bound = ',num2str(upper_bound)])

% generating cuts
average = 0;
C = zeros(n,N);
cut_value = zeros(1,N);
count = 0;
[U,D] = eig(full(X));
D = max(D,0);
Xsq = U*sqrt(D);
for k = 1:N,
    x = Xsq*randn(n,1);
    x = x./abs(x);
    x = x*x(1);
    value = 1/4*(sum(sum(W))-x'*W*x);
    average = average + value;
    Dx = C(:,1:count) - x*ones(1,count);
    if (count == 0) || (min(sum(abs(Dx))) > 0),
        count = count + 1;
        C(:,count) = x;
        cut_value(count) = value;
    end,
end,
average = average/N;
C = C(:,1:count);
cut_value = cut_value(1:count);
[cv,ind] = sort(cut_value);
C = C(:,ind);
figure
hold on
title('Values of random cuts')
xlabel('Number of cut')
ylabel('Value of cut')
plot(1:count,upper_bound*ones(1,count),'r:')
plot(1:count,average*ones(1,count),'b:')
plot(1:count,expectation*ones(1,count),'c:')
plot(1:count,0.87856*upper_bound*ones(1,count),'g:')
plot(1:count,cv,'k.')
legend('upper bound','average','expectation','\alpha * upper bound')

% plotting best cuts
for k = 1:min(m,count),
    x = C(:,end-k+1);
    figure
    hold on
    for l = 1:n,
        if x(l) == 1,
            plot(cos(phi(l)),sin(phi(l)),'b*')
        else
            plot(cos(phi(l)),sin(phi(l)),'r*')
        end,
    end,
    for k1 = 1:n-1,
        for k2 = k1+1:n,
            if W(k1,k2) == 1,
                if x(k1)*x(k2) == -1,
                    plot(cos(phi([k1 k2])),sin(phi([k1 k2])),'k')
                else
                    plot(cos(phi([k1 k2])),sin(phi([k1 k2])),'k:')
                end,
            end,
        end,
    end,
    axis([-1.1 1.1 -1.1 1.1])
    title(['Value of cut = ',num2str(cv(end-k+1))])
end,
