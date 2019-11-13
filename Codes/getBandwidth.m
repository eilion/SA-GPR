function [param] = getBandwidth(reg,data,param)

X = data.X;
N = size(X,1);

K = ceil(param.k*N/100);

index = (reg.H==0);
XX = X(index,:);

M = size(XX,1);

Dist = abs(X-XX');

AA = sort(Dist,2,'ascend');

param.bandwidth = AA(:,min(K,M));


end