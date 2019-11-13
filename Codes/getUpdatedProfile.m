function [param,Profile] = getUpdatedProfile(data,param,Profile,setting,dataType)

N = size(Profile,1);
L = length(data);

% Learn kernel hyperparameters:
MU = zeros(N,L);
NU = zeros(N,L);
QQ = zeros(N,L);
parfor ll = 1:L
    [param(ll),MU(:,ll),NU(:,ll),QQ(:,ll)] = learnKernelParam(data(ll),param(ll),Profile,setting,dataType);
end

% Update the profile:
AA = MU./NU;
BB = 1./NU;

index = (sum(QQ,2)>0);

Profile(index,2) = sum(AA(index,:).*QQ(index,:),2)./sum(BB(index,:).*QQ(index,:),2);
Profile(index,3) = sqrt(sum(QQ(index,:),2)./sum(BB(index,:).*QQ(index,:),2));


end