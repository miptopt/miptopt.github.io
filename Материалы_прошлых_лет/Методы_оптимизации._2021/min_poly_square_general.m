% minimization of random polynomial over the square

dp = 15; % dp is the degree of p 
d = 20; % d is the degree of the relaxation, d even and d >= dp
p = randn(dp+1); % coefficients of random polynomial

% compute polynomial and construct index
ind = zeros((d+1)*(d+2)/2,2);
count = 0;
x = -1:2^(-6):1;
y = -1:2^(-6):1;
[xx,yy] = meshgrid(x,y);
pval = zeros(length(x),length(y));
for k = 0:d,
    for l = 0:d,
        if k+l <= dp,
            pval = pval + p(k+1,l+1)*xx.^k.*yy.^l;
        end,
        if k+l <= d,
            count = count + 1;
            ind(count,:) = [k,l];
        end,
    end,
end,

% index of sigma_0 basis
d0 = d/2;
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

% index of sigma_1,...,sigma_4 basis
d1 = d/2-1;
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
            if sum(abs(ind1(l1,:)+ind1(l2,:)-ind(k,:))) == 0,
                cc1(l1,l2,k) = cc1(l1,l2,k) + 1;
                cc2(l1,l2,k) = cc2(l1,l2,k) + 1;
                cc3(l1,l2,k) = cc3(l1,l2,k) + 1;
                cc4(l1,l2,k) = cc4(l1,l2,k) + 1;
            end,
            if sum(abs(ind1(l1,:)+ind1(l2,:)+[1 0]-ind(k,:))) == 0,
                cc1(l1,l2,k) = cc1(l1,l2,k) - 1;
                cc2(l1,l2,k) = cc2(l1,l2,k) + 1;
            end,
            if sum(abs(ind1(l1,:)+ind1(l2,:)+[0 1]-ind(k,:))) == 0,
                cc3(l1,l2,k) = cc3(l1,l2,k) - 1;
                cc4(l1,l2,k) = cc4(l1,l2,k) + 1;
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
    maximize t
    p(ind(1,1)+1,ind(1,2)+1) - t == trace(C0*squeeze(cc0(:,:,1))) + trace(C1*squeeze(cc1(:,:,1))) + trace(C2*squeeze(cc2(:,:,1))) + trace(C3*squeeze(cc3(:,:,1))) + trace(C4*squeeze(cc4(:,:,1)));
    for k = 2:count,
        if sum(ind(k,:)) <= dp,
            p(ind(k,1)+1,ind(k,2)+1) == trace(C0*squeeze(cc0(:,:,k))) + trace(C1*squeeze(cc1(:,:,k))) + trace(C2*squeeze(cc2(:,:,k))) + trace(C3*squeeze(cc3(:,:,k))) + trace(C4*squeeze(cc4(:,:,k)));
        else
            trace(C0*squeeze(cc0(:,:,k))) + trace(C1*squeeze(cc1(:,:,k))) + trace(C2*squeeze(cc2(:,:,k))) + trace(C3*squeeze(cc3(:,:,k))) + trace(C4*squeeze(cc4(:,:,k))) == 0;
        end,
    end,
cvx_end

close all
figure
hold on
mesh(x,y,pval)
mesh(x,y,cvx_optval*ones(length(x),length(y)))
xlabel('y')
ylabel('x')
