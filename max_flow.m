% % max flow
n = 100;
W = zeros(n);
for k = 1:n-1
    for l = k+1:n
        if rand < 0.5
            W(k,l) = rand;
        end
    end
end
W = W + W';
tic,
cvx_begin
    variable F(n,n) skew_symmetric
    maximize( F(1,:)*ones(n,1) )
    F <= W
    for j = 2:n-1
        ones(1,n)*F(:,j) == 0
    end
cvx_end
toc,

% finding a minimum cut
W2 = W - F;
checked = 0;
reached = 1;
not_reached = 2:n;
while checked < length(reached)
    checked = checked + 1;
    ind = find(W2(reached(checked),not_reached) > 10^(-9));
    if ~isempty(ind)
        reached = [reached, not_reached(ind)];
        not_reached(ind) = [];
    end
end
reached,
