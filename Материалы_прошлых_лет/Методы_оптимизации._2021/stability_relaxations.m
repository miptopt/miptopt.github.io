% Lovasz and SOS bound on stability number

n = 50; % number of vertices
p = 0.4; % probability of edge
A = zeros(n);
for k = 1:n-1,
    for l = k+1:n,
        if rand < p,
            A(k,l) = 1;
            A(l,k) = 1;
        end,
    end,
end,

% Lovasz
cvx_begin
    variable X(n,n) semidefinite
    trace(X) == 1
    for k = 1:n-1,
        for l = k+1:n,
            if A(k,l),
                X(k,l) == 0;
            end,
        end,
    end,
    maximize(sum(sum(X)))
cvx_end
lovasz = cvx_optval;

% Lovasz dual
cvx_begin
    variable Y(n,n)
    variable lam
    for k = 1:n,
        for l = k:n,
            if ~A(k,l),
                Y(k,l) == 0;
            end,
        end,
    end,
    lam*eye(n) - Y - 1 == semidefinite(n)
    minimize(lam)
cvx_end
lovasz_dual = cvx_optval;

% SOS
cvx_begin
    variable Z(n,n) semidefinite
    variable la
    minimize(la)
    Z <= la*(eye(n)+A) - 1
cvx_end
sos = cvx_optval;

disp(['Lovasz primal: ',num2str(lovasz),';   Lovasz dual: ',num2str(lovasz_dual),';   SOS relaxation: ',num2str(cvx_optval)])
