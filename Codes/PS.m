function [ZZ] = PS(ZZ,data,Profile,setting)

nu = 6;
h = setting.h*(Profile(end,1)-Profile(1,1))/100;
nPaths = setting.nPaths;
M = setting.nParticles;
alpha = data.alpha;
beta = data.beta;

Profile(:,2) = Profile(:,2) + data.h;
% Profile(:,3) = 4/3*Profile(:,3);

X0 = data.X;
Y0 = data.Y;
Z0 = ZZ;
R = data.R;
N0 = size(X0,1);

QQ = linspace(X0(1),X0(end),100);
index = zeros(100,1);
for m = 1:100
    [~,index(m)] = min(abs(QQ(m)-X0));
end
index = unique(index);
N = length(index);

X = X0(index);
Y = Y0(index);
old_Z = Z0(index,:);

rand_seed = ceil(nPaths*rand(N,M));

PP = zeros(N,M);
WW = zeros(N,M);
% Forward Posteriors:
ZZ = (2*rand(M,5)-1)*h + old_Z(1,rand_seed(1,:))';
ZZ = ZZ(ZZ>=Profile(1,1)&ZZ<=Profile(end,1));
PP(1,:) = ZZ(1:M);
MM = interp1(Profile(:,1),Profile(:,2),PP(1,:));
NN = interp1(Profile(:,1),Profile(:,3),PP(1,:));
LL = -(nu+1)/2*log(1+((Y(1)-MM)./NN).^2/nu) - log(NN);
WW(1,:) = LL - max(LL);
for n = 1:N-1
    ZZ = (2*rand(M,5)-1)*h + old_Z(n+1,rand_seed(n+1,:))';
    ZZ = ZZ(ZZ>=Profile(1,1)&ZZ<=Profile(end,1));
    PP(n+1,:) = ZZ(1:M);
    MM = interp1(Profile(:,1),Profile(:,2),PP(n+1,:));
    NN = interp1(Profile(:,1),Profile(:,3),PP(n+1,:));
    LL = -(nu+1)/2*log(1+((Y(n+1)-MM)./NN).^2/nu) - log(NN);
    
    TT = (PP(n+1,:)-PP(n,:)')./(R.*(X(n+1)-X(n)));
    SS = TT.^(alpha-1).*exp(-TT*beta);
    SS(TT<=0) = 0;
    SS(TT>setting.max_AccRate) = 0;
    SS(TT<setting.min_AccRate) = 0;
    SS = LL + log(exp(WW(n,:))*SS);
    WW(n+1,:) = SS - max(SS);
end


Z = zeros(N,nPaths);
rand_seed = rand(N,nPaths);
% Backward Sampling:
W = exp(WW(N,:));
W = cumsum(W)/sum(W);
index = sum(rand_seed(N,:)'>W,2) + 1;
Z(N,:) = PP(N,index);
for nn = 1:N-1
    n = N - nn;
    
    TT = (Z(n+1,:)'-PP(n,:))./(R.*(X(n+1)-X(n)));
    TT = max(TT,0);
    W = (alpha-1)*log(TT) - TT*beta + WW(n,:);
    W(TT==0) = -inf;
    W(TT>setting.max_AccRate) = -inf;
    W(TT<setting.min_AccRate) = -inf;
    W = exp(W-max(W,[],2));
    W = cumsum(W,2)./sum(W,2);
    index = sum(rand_seed(n,:)'>W,2) + 1;
    Z(n,:) = PP(n,index);
end


% Interpolation:
Z0 = zeros(N0,nPaths);
for m = 1:nPaths
    Z0(:,m) = interp1(X,Z(:,m),X0);
end

ZZ = Z0;





end