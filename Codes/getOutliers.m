function [reg] = getOutliers(reg,data)

q = 0.01;
d = 3;

Y = data.Y;

MU = reg.mu;
NU = reg.nu + reg.lambda;

N = size(Y,1);

LL = zeros(N,2);
LL(:,1) = log(1-q) - 0.5*(Y-MU).^2./NU - 0.5*log(NU) - 0.5*log(2*pi);

AA = zeros(N,2);
AA(:,1) = - 0.5*(Y-MU-d*sqrt(NU)).^2./NU - 0.5*log(NU) - 0.5*log(2*pi);
AA(:,2) = - 0.5*(Y-MU+d*sqrt(NU)).^2./NU - 0.5*log(NU) - 0.5*log(2*pi);
amax = max(AA,[],2);
LL(:,2) = log(q/2) + amax + log(sum(exp(AA-amax),2));

amax = max(LL,[],2);
PT = exp(LL(:,2)-amax-log(sum(exp(LL-amax),2)));

% rand_seed = rand(N,1);

reg.H = (PT>0.5);


end