% truss topology design
% knots arranged in rectangle
% left-most knots fixed
% force (0,-1) applied to central one of the right-most knots 
n = 5; % number of knots in vertical direction, odd
m = 3; % number of knots in horizontal direction
f = zeros(2*(m-1)*n,1);
f(2*(m-2)*n+n+1) = -1;
v = zeros(2*m*n,1); % position of knots
count = 0;
for k = 1:m
    for l = 1:n
        v(2*count+(1:2)) = [k; l];
        count = count + 1;
    end
end
count = 0;
ind = zeros((m-1)*n^2,2);
for k = 1:m-1
    for l1 = 1:n
        for l2 = 1:n
            if abs(l1-l2) <= 2
                count = count + 1;
                ind(count,1) = (k-1)*n+l1;
                ind(count,2) = k*n+l2;
            end
        end
    end
end
ind = ind(1:count,:); % indices of starting and ending knot for bars
B = zeros(2*m*n,count);
for k = 1:count
    b = v(2*ind(k,1)+(-1:0)) - v(2*ind(k,2)+(-1:0));
    b = b/(b'*b);
    B(2*ind(k,1)+(-1:0),k) = b;
    B(2*ind(k,2)+(-1:0),k) = -b;
end
cvx_begin
    variable ma(count) nonnegative
    variable M(2*(m-1)*n+1,2*(m-1)*n+1) semidefinite
    variable t
    minimize(t)
    M(end,end) == t
    M(end,1:2*(m-1)*n) == -f';
    M(1:2*(m-1)*n,end) == -f;
    M(1:2*(m-1)*n,1:2*(m-1)*n) == B(2*n+(1:2*(m-1)*n),:)*diag(ma)*B(2*n+(1:2*(m-1)*n),:)';
    sum(ma) == 1
cvx_end
close all
figure
hold on
for k = 1:m
    for l = 1:n
        plot(k,l,'k*')
    end
end
plot([1 1],[1 n],'k-')
for k = 1:count
    if ma(k) > exp(-15)
        v1 = v(2*ind(k,1)+(-1:0));
        v2 = v(2*ind(k,2)+(-1:0));
        plot([v1(1) v2(1)],[v1(2) v2(2)],'b-')
    end
end

