function [data] = shortHMM(data,Profile,setting)

M = size(Profile,1);

CC = Profile(:,1)';
AA = Profile(:,2)';
BB = Profile(:,3)';

nu = setting.nu;
nPaths = setting.nPaths;
alpha = setting.alpha;
beta = setting.beta;

X = data.X;
Y = data.Y;
R = data.R;

N = size(X,1);

QQ = linspace(X(1),X(end),100);
index = zeros(100,1);
for m = 1:100
    [~,index(m)] = min(abs(QQ(m)-X));
end
index = unique(index);
NN = length(index);

% index = ceil(N*rand(100,1));
% index = sort(index,'ascend');
% index(1) = 1;
% index(100) = N;
% index = unique(index);
% NN = length(index);

XX = X(index);
YY = Y(index);
RR = interp1(R(:,1),R(:,2),CC');


EM = -(nu+1)/2*log(1+((YY-AA)./BB).^2/nu) - log(BB);
EM = exp(EM-max(EM,[],2));

FT = zeros(NN,M);
% Forward Algorithm:
LL = EM(1,:)/M;
FT(1,:) = LL/max(LL);
for n = 1:NN-1
    TT = (CC-CC')./(RR.*(XX(n+1)-XX(n)));
    SS = TT.^(alpha-1).*exp(-TT/beta)./(RR.*(X(n+1)-X(n)));
    SS(CC-CC'<=0) = 0;
    % SS(M,M) = 1;
    % SS = SS./sum(SS,2);
    % SS = SS*0.99;
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
    SS = TT.^(alpha-1).*exp(-TT/beta)./(RR.*(X(n+1)-X(n)));
    SS(CC-CC'<=0) = 0;
    % SS(M,M) = 1;
    % SS = SS./sum(SS,2);
    % SS = SS*0.99;
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

data.Z = Z;


end