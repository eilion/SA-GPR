function [Z] = getOutliers(Y,MU,NU,RR,setting)

q = setting.delta;
d = setting.d;

N = size(Y,1);

LL = zeros(N,2);
LL(:,1) = log(1-q) - 0.5*(Y-MU).^2./(NU+RR) - 0.5*log(NU+RR) - 0.5*log(2*pi);

AA = zeros(N,2);
AA(:,1) = - 0.5*(Y-MU-d*sqrt(NU+RR)).^2./(NU+RR) - 0.5*log(NU+RR) - 0.5*log(2*pi);
AA(:,2) = - 0.5*(Y-MU+d*sqrt(NU+RR)).^2./(NU+RR) - 0.5*log(NU+RR) - 0.5*log(2*pi);
amax = max(AA,[],2);
LL(:,2) = log(q/2) + amax + log(sum(exp(AA-amax),2));

amax = max(LL,[],2);
PT = exp(LL(:,2)-amax-log(sum(exp(LL-amax),2)));

Z = (PT>0.5);


end