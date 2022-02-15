% minimizes random monic univariate polynomial

d = 10; % degree of polynomial is 2d
p = randn(1,2*d+1); % random coefficient vector
p(end) = 1; % coefficient at largest power equals 1

I = eye(2*d+1);
cvx_begin sdp
    variable t
    variable P(d+1,d+1) semidefinite
    sk = zeros(1,2*d + 1);
    for l1 = 0:d
        for l2 = 0:d
            sk = sk + I(l1+l2+1,:)*P(l1+1,l2+1);
        end
    end
    sk = sk + t*[1, zeros(1,2*d)];
    sk == p;
    maximize(t)
cvx_end

% plot polynomial
close all
x = -3:2^(-9):3;
pval = zeros(1,length(x));
for k = 0:2*d,
    pval = pval + p(k+1)*x.^k;
end,
figure
hold on
plot(x,pval)
plot(x,t*ones(1,length(x)),'k:')
xlabel('x')
ylabel('p(x)')
legend('value of p','minimum')
axis([min(x), max(x), t-0.5, t+1.5])
figure
plot(x,log(pval-t))
xlabel('x')
ylabel('log(p(x) - p_{min}')