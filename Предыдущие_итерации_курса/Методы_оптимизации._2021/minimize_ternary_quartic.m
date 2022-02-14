% minimizes random monic ternary quartic

p = randn(5); % random coefficient vector
p(1,5) = 1; % coefficients at x^4 and y^4 equal 1
p(5,1) = 1;

cvx_begin sdp
    variable t
    variable P(6,6) semidefinite
    p(1,1) == P(3,3) + t
    p(1,2) == P(3,5) + P(5,3)
    p(1,3) == P(1,3) + P(3,1) + P(5,5)
    p(1,4) == P(1,5) + P(5,1)
    p(1,5) == P(1,1)
    p(2,1) == P(3,4) + P(4,3)
    p(2,2) == P(3,6) + P(6,3) + P(4,5) + P(5,4)
    p(2,3) == P(1,4) + P(4,1) + P(5,6) + P(6,5)
    p(2,4) == P(1,6) + P(6,1)
    p(3,1) == P(2,3) + P(3,2) + P(4,4)
    p(3,2) == P(2,5) + P(5,2) + P(4,6) + P(6,4)
    p(3,3) == P(1,2) + P(2,1) + P(6,6)
    p(4,1) == P(2,4) + P(4,2)
    p(4,2) == P(2,6) + P(6,2)
    p(5,1) == P(2,2)
    maximize(t)
cvx_end

% plot polynomial
close all
x = -2:2^(-7):2;
y = -2:2^(-7):2;
[xx,yy] = meshgrid(x,y);
pval = p(1,1)*ones(length(x),length(y)) + p(1,2)*xx + p(1,3)*xx.^2 + p(1,4)*xx.^3 + p(1,5)*xx.^4 + p(2,1)*yy + p(2,2)*xx.*yy + p(2,3)*xx.^2.*yy + p(2,4)*xx.^3.*yy + p(3,1)*yy.^2 + p(3,2)*xx.*yy.^2 + p(3,3)*xx.^2.*yy.^2 + p(4,1)*yy.^3 + p(4,2)*xx.*yy.^3 + p(5,1)*yy.^4;
figure
hold on
mesh(x,y,t*ones(length(x),length(y)))
mesh(x,y,pval)
xlabel('x')
ylabel('y')
zlabel('p(x,y)')
figure
hold on
mesh(x,y,t*ones(length(x),length(y)))
zlim([t - 0.5, t + 5.5])
mesh(x,y,pval)
xlabel('x')
ylabel('y')
zlabel('p(x,y)')
