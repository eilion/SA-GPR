function [reg] = getInitialReg(reg,data,param)

index = (reg.H==0);

XX = data.X(index);
YY = data.Y(index);
RR = reg.lambda(index);

N = length(YY);

K = getCov(XX,XX,param.eta,param.xi) + param.lambda^2*eye(N) + param.lambda0^2*eye(N);

B0 = diag(RR) + K;
B2 = B0\eye(N);
B1 = B2*YY;

K1 = getCov(data.X,XX,param.eta,param.xi);
K2 = param.eta^2 + param.lambda^2 + param.lambda0^2;

reg.mu = K1*B1;
reg.nu = K2 - sum((K1*B2).*K1,2);


end