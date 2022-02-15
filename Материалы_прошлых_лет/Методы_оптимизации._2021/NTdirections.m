function [dx,ds,px,ps] = NTdirections(A,b,c,x,s)
% function [dx,ds,px,ps] = NTdirections(A,b,c,x,s)
%
% computes Nesterov-Todd directions for linear program with data (A,b,c)
% at primal-dual point (x,s)
% centering direction (dx,ds)
% affine scaling direction (px,ps)
n = size(A,2);
m = size(A,1);
Axs = A*diag(x./s);
H = inv(Axs*A');
invs = 1./s;
mu = (x.'*s)/n;
dx = -(eye(n)-Axs'*H*A)*(x-mu*invs);
px = -x+Axs'*H*b;
ds = -A'*H*b+mu*A'*H*A*invs;
ps = -A'*H*b;