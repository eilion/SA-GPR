function [reg] = learnVar_NW(reg,data,param)

X = data.X;
Y = data.Y;

BW = param.bandwidth;

index = (reg.H==0);
XX = X(index,:);

Dist2 = (X-XX').^2;
WW = - 0.5*Dist2./(BW.^2) - log(BW);
WW = WW - max(WW,[],2);
WW = exp(WW);
LL = WW./sum(WW,2);

Y = Y(index)';
mu = reg.mu(index)';
nu = reg.nu(index)';

reg.lambda = sum(((Y-mu).^2+nu).*LL,2);


end