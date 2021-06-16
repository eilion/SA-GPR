function [ZZZ] = initializeAlignment(data,Profile,setting)

M = 300;

shortProfile = zeros(M,3);
shortProfile(:,1) = linspace(Profile(1,1),Profile(end,1),M);
shortProfile(:,2) = interp1(Profile(:,1),Profile(:,2),shortProfile(:,1));
shortProfile(:,3) = interp1(Profile(:,1),Profile(:,3),shortProfile(:,1));

CC = shortProfile(:,1)';
AA = shortProfile(:,2)' + data.h;
% BB = 4/3*shortProfile(:,3)';
BB = shortProfile(:,3)';

nPaths = setting.nPaths;
alpha = data.alpha;
beta = data.beta;
nu = 6;

X = data.X;
Y = data.Y;
RR = data.R;

N = size(X,1);

QQ = linspace(X(1),X(end),100);
index = zeros(100,1);
for m = 1:100
    [~,index(m)] = min(abs(QQ(m)-X));
end
index = unique(index);
NN = length(index);

XX = X(index);
YY = Y(index);

EM = -(nu+1)/2*log(1+((YY-AA)./BB).^2/nu) - log(BB);
EM = exp(EM-max(EM,[],2));

FT = zeros(NN,M);
% Forward Algorithm:
LL = EM(1,:)/M;
FT(1,:) = LL/max(LL);
for n = 1:NN-1
    TT = (CC-CC')./(RR.*(XX(n+1)-XX(n)));
    SS = TT.^(alpha-1).*exp(-TT*beta);
    SS(CC-CC'<=0) = 0;
    % SS(M,M) = 1;
    % SS = SS./sum(SS,2);
    % SS = SS*0.99;
    SS(TT>setting.max_AccRate) = 0;
    SS = SS + eye(M)*0.01;
    
    LL = EM(n+1,:).*(FT(n,:)*SS);
    FT(n+1,:) = LL/max(LL);
end

ZZ = zeros(NN,nPaths);
INDEX = zeros(NN,nPaths);
rand_seed = rand(NN,nPaths);
% Backward Sampling:
W = cumsum(FT(NN,:))/sum(FT(NN,:));
index = sum(rand_seed(NN,:)'>W,2) + 1;
ZZ(NN,:) = CC(index);
INDEX(NN,:) = index;
for nn = 1:NN-1
    n = NN - nn;
    
    TT = (CC-CC')./(RR.*(XX(n+1)-XX(n)));
    SS = TT.^(alpha-1).*exp(-TT*beta);
    SS(CC-CC'<=0) = 0;
    % SS(M,M) = 1;
    % SS = SS./sum(SS,2);
    % SS = SS*0.99;
    SS(TT>setting.max_AccRate) = 0;
    SS = SS + eye(M)*0.01;
    
    QQ = SS(:,INDEX(n+1,:))'.*FT(n,:);
    W = cumsum(QQ,2)./sum(QQ,2);
    index = sum(rand_seed(n,:)'>W,2) + 1;
    ZZ(n,:) = CC(index);
    INDEX(n,:) = index;
end


% Interpolation:
Z = zeros(N,nPaths);
for m = 1:nPaths
    Z(:,m) = interp1(XX,ZZ(:,m),X);
end

ZZZ = Z;


end