% ressource allocation
%
% generating random problem
K = 100; % number of types of raw material
n = 50; % number of products
p = abs(randn(1,n));
r = 100*n*abs(randn(K,1));
A = rand(K,n);
% solving continuous LP
cvx_begin
    variable x(n,1) nonnegative
    maximize( p*x )
    A*x <= r
cvx_end
% output solution
close all
plot(1:n,x,'o')
xlabel('Number of product')
ylabel('Quantity of product')
grid on
fprintf("Continuous solution:\n");
for k = 1:n
    if x(k) > 10^(-6)
        fprintf("Product %d: produce %5e units\n",k,x(k));
    end
end
% creating mps file containing the MILP
milpRA = BuildMPS(A, r, [], [], -p, zeros(n,1), ones(n,1)*Inf, [], 'Integer', 1:n);
SaveMPS('milpRA.mps', milpRA);
% commands for running in gurobi interactive shell
% m = read('C:/Users/Hildebrand/Desktop/teaching/EMO/mcsi2021/matlab/milpRA') 
% m.optimize()
% m.printAttr('X')