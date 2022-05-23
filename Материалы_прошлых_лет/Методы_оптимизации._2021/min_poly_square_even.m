% minimization of random even polynomial over the square

d = 10; % d is the degree of p and must be even, degree of relaxation is equal to degree of p 
p = randn(d+1); % coefficients of random polynomial

% compute polynomial and construct index
ind = zeros((d+1)*(d+2)/2,2);
count = 0;
x = -1:2^(-6):1;
y = -1:2^(-6):1;
[xx,yy] = meshgrid(x,y);
pval = zeros(length(x),length(y));
for k = 0:d,
    for l = 0:d,
        if k+l <= d,
            pval = pval + p(k+1,l+1)*xx.^k.*yy.^l;
            count = count + 1;
            ind(count,:) = [k,l];
        end,
    end,
end,

% index of sigma_0 basis
d0 = floor(d/2);
ind0 = zeros((d0+1)*(d0+2)/2,2);
count0 = 0;
for k = 0:d0,
    for l = 0:d0,
        if k+l <= d0,
            count0 = count0 + 1;
            ind0(count0,:) = [k l];
        end,
    end,
end,

% index of sigma_1,...,sigma_6 basis
d1 = floor(d/2-1);
ind1 = zeros((d1+1)*(d1+2)/2,2);
count1 = 0;
for k = 0:d1,
    for l = 0:d1,
        if k+l <= d1,
            count1 = count1 + 1;
            ind1(count1,:) = [k l];
        end,
    end,
end,

% coefficient matrices
cc0 = zeros(count0,count0,count);
cc1 = zeros(count1,count1,count);
cc2 = zeros(count1,count1,count);
cc3 = zeros(count1,count1,count);
cc4 = zeros(count1,count1,count);
cc5 = zeros(count1,count1,count);
cc6 = zeros(count1,count1,count);
for k = 1:count,
    for l1 = 1:count0,
        for l2 = 1:count0,
            if sum(abs(ind0(l1,:)+ind0(l2,:)-ind(k,:))) == 0,
                cc0(l1,l2,k) = cc0(l1,l2,k) + 1;
            end,
        end,
    end,
    for l1 = 1:count1,
        for l2 = 1:count1,
            if sum(abs(ind1(l1,:)+ind1(l2,:)-ind(k,:))) == 0, % coef 1
                cc1(l1,l2,k) = cc1(l1,l2,k) + 1;
                cc2(l1,l2,k) = cc2(l1,l2,k) + 1;
                cc3(l1,l2,k) = cc3(l1,l2,k) + 1;
                cc4(l1,l2,k) = cc4(l1,l2,k) + 1;
                cc5(l1,l2,k) = cc5(l1,l2,k) + 1;
                cc6(l1,l2,k) = cc6(l1,l2,k) + 1;
            end,
            if sum(abs(ind1(l1,:)+ind1(l2,:)+[1 0]-ind(k,:))) == 0, % coef x
                cc2(l1,l2,k) = cc2(l1,l2,k) - 1;
                cc3(l1,l2,k) = cc3(l1,l2,k) - 1;
                cc4(l1,l2,k) = cc4(l1,l2,k) + 1;
                cc5(l1,l2,k) = cc5(l1,l2,k) + 1;
            end,
            if sum(abs(ind1(l1,:)+ind1(l2,:)+[2 0]-ind(k,:))) == 0, % coef x^2
                cc1(l1,l2,k) = cc1(l1,l2,k) - 1;
            end,
            if sum(abs(ind1(l1,:)+ind1(l2,:)+[0 1]-ind(k,:))) == 0, % coef y
                cc2(l1,l2,k) = cc2(l1,l2,k) - 1;
                cc3(l1,l2,k) = cc3(l1,l2,k) + 1;
                cc4(l1,l2,k) = cc4(l1,l2,k) - 1;
                cc5(l1,l2,k) = cc5(l1,l2,k) + 1;
            end,
            if sum(abs(ind1(l1,:)+ind1(l2,:)+[0 2]-ind(k,:))) == 0, % coef y^2
                cc6(l1,l2,k) = cc6(l1,l2,k) - 1;
            end,
            if sum(abs(ind1(l1,:)+ind1(l2,:)+[1 1]-ind(k,:))) == 0, % coef x*y
                cc2(l1,l2,k) = cc2(l1,l2,k) + 1;
                cc3(l1,l2,k) = cc3(l1,l2,k) - 1;
                cc4(l1,l2,k) = cc4(l1,l2,k) - 1;
                cc5(l1,l2,k) = cc5(l1,l2,k) + 1;
            end,
        end,
    end,
end,
    
cvx_begin
    variable t
    variable C0(count0,count0) semidefinite
    variable C1(count1,count1) semidefinite
    variable C2(count1,count1) semidefinite
    variable C3(count1,count1) semidefinite
    variable C4(count1,count1) semidefinite
    variable C5(count1,count1) semidefinite
    variable C6(count1,count1) semidefinite
    maximize t
    p(ind(1,1)+1,ind(1,2)+1) - t == trace(C0*squeeze(cc0(:,:,1))) + trace(C1*squeeze(cc1(:,:,1))) + trace(C2*squeeze(cc2(:,:,1))) + trace(C3*squeeze(cc3(:,:,1))) + trace(C4*squeeze(cc4(:,:,1))) + trace(C5*squeeze(cc5(:,:,1))) + trace(C6*squeeze(cc6(:,:,1)));
    for k = 2:count,
        p(ind(k,1)+1,ind(k,2)+1) == trace(C0*squeeze(cc0(:,:,k))) + trace(C1*squeeze(cc1(:,:,k))) + trace(C2*squeeze(cc2(:,:,k))) + trace(C3*squeeze(cc3(:,:,k))) + trace(C4*squeeze(cc4(:,:,k))) + trace(C5*squeeze(cc5(:,:,k))) + trace(C6*squeeze(cc6(:,:,k)));
    end,
cvx_end

close all
figure
hold on
mesh(x,y,pval)
mesh(x,y,cvx_optval*ones(length(x),length(y)))
