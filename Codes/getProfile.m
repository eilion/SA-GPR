function [Profile] = getProfile(reg,data,param,setting)

Profile = zeros(1000,3);

st = 0.5*((2+setting.margin*2)*data.X(1)-setting.margin*2*data.X(end));
ed = 0.5*((2+setting.margin*2)*data.X(end)-setting.margin*2*data.X(1));
% st = 0.5*(2.1*data.X(1)-0.1*data.X(end));
% ed = 0.5*(2.1*data.X(end)-0.1*data.X(1));
% st = 0.5*(3*data.X(1)-data.X(end));
% ed = 0.5*(3*data.X(end)-data.X(1));

X = linspace(st,ed,1000)';

index = (reg.H==0);

XX = data.X(index);
YY = data.Y(index);
RR = reg.lambda(index);

N = length(YY);

K = getCov(XX,XX,param.eta,param.xi) + param.lambda^2*eye(N) + param.lambda0^2*eye(N);

B0 = diag(RR) + K;
B2 = B0\eye(N);
B1 = B2*YY;

K1 = getCov(X,XX,param.eta,param.xi);
K2 = param.eta^2 + param.lambda^2 + param.lambda0^2;

MU = K1*B1;
NU = K2 - sum((K1*B2).*K1,2);

LA = interp1(data.X,reg.lambda,X);
LA(X<data.X(1)) = reg.lambda(1);
LA(X>data.X(end)) = reg.lambda(end);


Profile(:,1) = X;
Profile(:,2) = MU;
Profile(:,3) = sqrt(NU+LA);


end