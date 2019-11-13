function [data] = MCMC(data,Profile,setting)

nu = setting.nu;
alpha = setting.alpha;
beta = setting.beta;
MH_Iters = setting.MH_Iters;

[N,M] = size(data.Z);
X = data.X;
Y = data.Y;
Z = data.Z;
R = data.R;

index0 = 2:2:N-1;
index1 = 3:2:N-1;

index0p = 1:2:N-2;
index0n = 3:2:N;
index1p = 2:2:N-2;
index1n = 4:2:N;

for r = 1:MH_Iters
    RAND = log(rand(N,M));
    
    % n = 1:
    ZZ_old = Z(1,:);
    ZZ_new = normrnd(0,1,[1,M]).*(Z(2,:)-Z(1,:))/8 + Z(1,:);
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_new);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_new);
    loglik_new = -(nu+1)/2*log(1+((Y(1)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),ZZ_new);
    TT = (Z(2,:)-ZZ_new)./(RR.*(X(2)-X(1)));
    TT = max(TT,0);
    loglik_new = loglik_new + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(2)-X(1)));
    
    QQ = max((Z(2,:)-ZZ_new)/8,0);
    loglik_new = loglik_new - (ZZ_new-ZZ_old).^2.*((Z(2,:)-ZZ_new)/8).^2 - log(QQ);
    loglik_new(isnan(loglik_new)==1) = -inf;
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_old);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_old);
    loglik_old = -(nu+1)/2*log(1+((Y(1)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),ZZ_old);
    TT = (Z(2,:)-ZZ_old)./(RR.*(X(2)-X(1)));
    TT = max(TT,0);
    loglik_old = loglik_old + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(2)-X(1)));
    
    QQ = max((Z(2,:)-ZZ_old)/8,0);
    loglik_old = loglik_old - (ZZ_new-ZZ_old).^2.*((Z(2,:)-ZZ_old)/8).^2 - log(QQ);
    loglik_old(isnan(loglik_old)==1) = -inf;
    
    ALPHA = loglik_new - loglik_old;
    ALPHA(isnan(ALPHA)==1) = -inf;
    index = (ALPHA>RAND(1,:));
    Z(1,index) = ZZ_new(index);
    
    % n in index0:
    ZZ_old = Z(index0,:);
    ZZ_new = rand(length(index0),M).*(Z(index0n,:)-Z(index0p,:)) + Z(index0p,:);
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_new);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_new);
    loglik_new = -(nu+1)/2*log(1+((Y(index0)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),ZZ_new);
    TT = (Z(index0n,:)-ZZ_new)./(RR.*(X(index0n)-X(index0)));
    TT = max(TT,0);
    loglik_new = loglik_new + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index0n)-X(index0)));
    
    RR = interp1(R(:,1),R(:,2),Z(index0p,:));
    TT = (ZZ_new-Z(index0p,:))./(RR.*(X(index0)-X(index0p)));
    TT = max(TT,0);
    loglik_new = loglik_new + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index0)-X(index0p)));
    
    loglik_new(isnan(loglik_new)==1) = -inf;
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_old);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_old);
    loglik_old = -(nu+1)/2*log(1+((Y(index0)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),ZZ_old);
    TT = (Z(index0n,:)-ZZ_old)./(RR.*(X(index0n)-X(index0)));
    TT = max(TT,0);
    loglik_old = loglik_old + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index0n)-X(index0)));
    
    RR = interp1(R(:,1),R(:,2),Z(index0p,:));
    TT = (ZZ_old-Z(index0p,:))./(RR.*(X(index0)-X(index0p)));
    TT = max(TT,0);
    loglik_old = loglik_old + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index0)-X(index0p)));
    
    loglik_old(isnan(loglik_old)==1) = -inf;
    
    ALPHA = loglik_new - loglik_old;
    ALPHA(isnan(ALPHA)==1) = -inf;
    index = (ALPHA>RAND(index0,:));
    UU = Z(index0,:);
    UU(index) = ZZ_new(index);
    Z(index0,:) = UU;
    
    % n in index1:
    ZZ_old = Z(index1,:);
    ZZ_new = rand(length(index1),M).*(Z(index1n,:)-Z(index1p,:)) + Z(index1p,:);
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_new);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_new);
    loglik_new = -(nu+1)/2*log(1+((Y(index1)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),ZZ_new);
    TT = (Z(index1n,:)-ZZ_new)./(RR.*(X(index1n)-X(index1)));
    TT = max(TT,0);
    loglik_new = loglik_new + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index1n)-X(index1)));
    
    RR = interp1(R(:,1),R(:,2),Z(index1p,:));
    TT = (ZZ_new-Z(index1p,:))./(RR.*(X(index1)-X(index1p)));
    TT = max(TT,0);
    loglik_new = loglik_new + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index1)-X(index1p)));
    
    loglik_new(isnan(loglik_new)==1) = -inf;
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_old);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_old);
    loglik_old = -(nu+1)/2*log(1+((Y(index1)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),ZZ_old);
    TT = (Z(index1n,:)-ZZ_old)./(RR.*(X(index1n)-X(index1)));
    TT = max(TT,0);
    loglik_old = loglik_old + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index1n)-X(index1)));
    
    RR = interp1(R(:,1),R(:,2),Z(index1p,:));
    TT = (ZZ_old-Z(index1p,:))./(RR.*(X(index1)-X(index1p)));
    TT = max(TT,0);
    loglik_old = loglik_old + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(index1)-X(index1p)));
    
    loglik_old(isnan(loglik_old)==1) = -inf;
    
    ALPHA = loglik_new - loglik_old;
    ALPHA(isnan(ALPHA)==1) = -inf;
    index = (ALPHA>RAND(index1,:));
    UU = Z(index1,:);
    UU(index) = ZZ_new(index);
    Z(index1,:) = UU;
    
    % n = N:
    ZZ_old = Z(N,:);
    ZZ_new = normrnd(0,1,[1,M]).*(Z(N,:)-Z(N-1,:))/8 + Z(N,:);
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_new);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_new);
    loglik_new = -(nu+1)/2*log(1+((Y(N)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),Z(N-1,:));
    TT = (ZZ_new-Z(N-1,:))./(RR.*(X(N)-X(N-1)));
    TT = max(TT,0);
    loglik_new = loglik_new + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(N)-X(N-1)));
    
    QQ = max((ZZ_new-Z(N-1,:))/8,0);
    loglik_new = loglik_new - (ZZ_new-ZZ_old).^2.*((ZZ_new-Z(N-1,:))/8).^2 - log(QQ);
    loglik_new(isnan(loglik_new)==1) = -inf;
    
    MM = interp1(Profile(:,1),Profile(:,2),ZZ_old);
    NN = interp1(Profile(:,1),Profile(:,3),ZZ_old);
    loglik_old = -(nu+1)/2*log(1+((Y(N)-MM)./NN).^2/nu) - log(NN);
    
    RR = interp1(R(:,1),R(:,2),Z(N-1,:));
    TT = (ZZ_old-Z(N-1,:))./(RR.*(X(N)-X(N-1)));
    TT = max(TT,0);
    loglik_old = loglik_old + (alpha-1)*log(TT) - TT/beta - log(RR.*(X(N)-X(N-1)));
    
    QQ = max((ZZ_old-Z(N-1,:))/8,0);
    loglik_old = loglik_old - (ZZ_new-ZZ_old).^2.*((ZZ_old-Z(N-1,:))/8).^2 - log(QQ);
    loglik_old(isnan(loglik_old)==1) = -inf;
    
    ALPHA = loglik_new - loglik_old;
    ALPHA(isnan(ALPHA)==1) = -inf;
    index = (ALPHA>RAND(N,:));
    Z(N,index) = ZZ_new(index);
end

data.Z = Z;


end