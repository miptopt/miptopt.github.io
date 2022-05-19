% Nemirovskis approximation of the unit disc by linear constraints
n = 10; % degree of approximation
N = 20; % number of trials
results = zeros(1,N);
M = zeros(2,2,n+1);
for i = 1:n
    M(:,:,i) = [cos(pi/2^(i+1)), sin(pi/2^(i+1)); -sin(pi/2^(i+1)), cos(pi/2^(i+1))];
end
M(:,:,n+1) = [1, 0; -tan(pi/2^(n+2)), 1];
for k = 1:N
    phi = rand*2*pi;
    c = [cos(phi), sin(phi), zeros(1,4*n+2)];
    cvx_begin
        variable x(4*n+4,1)
        maximize ( c*x - 1 )
        for i = 1:n
            x(2*(n+1+i) + (1:2)) == squeeze(M(:,:,i))*x(2*i + (1:2))
            x(2*i+3) == x(2*(n+i)+3)
            x(2*i+4) >= x(2*(n+i)+4)
            x(2*i+4) >= -x(2*(n+i)+4)
        end
        x(3) >= x(1)
        x(3) >= -x(1)
        x(4) >= x(2)
        x(4) >= -x(2)
        squeeze(M(:,:,n+1))*x(2*(n+1) + (1:2)) <= [1; 0]
    cvx_end
    results(k) = cvx_optval;
end
plot(results)