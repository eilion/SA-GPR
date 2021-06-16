function [Profile] = initializeProfile(data,setting)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;

lambda0 = 1e-2;

max_iters = 30000;

Profile = zeros(1000,3);
Profile(:,1) = linspace(-1,1,1000)';

X = data(1).X;
Y = data(1).Y;

mean_Y = mean(Y);
stdv_Y = sqrt(var(Y));

Y = (Y-mean_Y)/stdv_Y;

N = size(X,1);
N0 = min(N,100);

Z = zeros(N,1);

BW = zeros(N,4);
DIST = abs(X-X');
DIST = sort(DIST,2,'ascend');
for k = 1:4
    K = ceil(5*k*N/100);
    BW(:,k) = (DIST(:,K)+DIST(:,K+1))/2;
end

ETA = 1;
GAMMA = 1;
if strcmp(setting.mode,'graph') == 1
    LAMBDA = 0.01;
else
    LAMBDA = 0.1;
end

if strcmp(setting.kernel,'SE[SE]') == 1
    DELTA = 1;
    
    Mw_DELTA = 0;
    Vw_DELTA = 0;
end

Mw = zeros(3,1);
Vw = zeros(3,1);

ST = ceil(max(0,(N-200))*rand(max_iters,1));

for r = 1:max_iters
    XX = X(ST(r)+1:ST(r)+N0);
    YY = Y(ST(r)+1:ST(r)+N0);
    ZZ = Z(ST(r)+1:ST(r)+N0);
    
    XX = XX(ZZ==0);
    YY = YY(ZZ==0);
    NN0 = sum(ZZ==0);
    
    if strcmp(setting.mode,'graph') == 1 || strcmp(setting.variance,'homoscedastic') == 1
        RR = LAMBDA^2*ones(NN0,1);
    elseif strcmp(setting.variance,'heteroscedastic') == 1
        if r == 1
            R = LAMBDA^2*ones(N,1);
        end
        RR = R(ST(r)+1:ST(r)+N0);
        RR = RR(ZZ==0);
    end
    
    if rem(r,1000) == 0
        if strcmp(setting.mode,'graph') == 0
            if strcmp(setting.variance,'homoscedastic') == 1
                if strcmp(setting.kernel,'SE[SE]') == 1
                    [MU,NU] = getReg(X,X,Y,Z,ETA,[GAMMA;DELTA],LAMBDA^2*ones(N,1),setting);
                else
                    [MU,NU] = getReg(X,X,Y,Z,ETA,GAMMA,LAMBDA^2*ones(N,1),setting);
                end
                Z = getOutliers(Y,MU,NU,LAMBDA^2*ones(N,1),setting);
            elseif strcmp(setting.variance,'heteroscedastic') == 1
                if strcmp(setting.kernel,'SE[SE]') == 1
                    [MU,NU] = getReg(X,X,Y,Z,ETA,[GAMMA;DELTA],R,setting);
                else
                    [MU,NU] = getReg(X,X,Y,Z,ETA,GAMMA,R,setting);
                end
                Z = getOutliers(Y,MU,NU,R,setting);
                
                LOGLIK = zeros(4,1);
                for k = 1:4
                    LL = - 0.5*((X-X')./BW(:,k)').^2 - log(BW(:,k)');
                    LL = LL - diag(inf*ones(N,1));
                    LL = LL(:,Z==0);
                    LL = LL - max(LL,[],2);
                    LL = exp(LL);
                    
                    R = sum(((Y(Z==0)'-MU(Z==0)').^2+NU(Z==0)').*LL,2)./sum(LL,2);
                    
                    LOGLIK(k) = sum(-0.5*(Y(Z==0)-MU(Z==0)).^2./(NU(Z==0)+R(Z==0))-0.5*log(NU(Z==0)+R(Z==0))-0.5*log(2*pi));
                end
                [~,K] = max(LOGLIK);
                
                LL = - 0.5*((X-X')./BW(:,K)').^2 - log(BW(:,K)');
                LL = LL(:,Z==0);
                LL = LL - max(LL,[],2);
                LL = exp(LL);
                
                R = sum(((Y(Z==0)'-MU(Z==0)').^2+NU(Z==0)').*LL,2)./sum(LL,2);
            end
        end
    end
    
    PDEV = zeros(3,1);
    
    if strcmp(setting.kernel,'SE') == 1
        K_XX = exp(-GAMMA^2*(XX-XX').^2);
    elseif strcmp(setting.kernel,'OU') == 1
        K_XX = exp(-GAMMA^2*abs(XX-XX'));
    elseif strcmp(setting.kernel,'SE[SE]') == 1
        K_XX = (1+DELTA^2*(1-exp(-GAMMA^2*(XX-XX').^2))).^(-1/2);
    end
    
    WW_inv = (ETA^2*K_XX+diag(RR)+lambda0^2*eye(NN0))\eye(NN0);
    ALPHA = WW_inv*YY;
    ALPHA = (ALPHA*ALPHA'-WW_inv)';
    
    % ETA:
    PK = 2*ETA*K_XX;
    PDEV(1) = 0.5*sum(sum(PK.*ALPHA,2));
    
    % GAMMA:
    if strcmp(setting.kernel,'SE') == 1
        PK = ETA^2*K_XX.*(-2*GAMMA*(XX-XX').^2);
    elseif strcmp(setting.kernel,'OU') == 1
        PK = ETA^2*K_XX.*(-2*GAMMA*abs(XX-XX'));
    elseif strcmp(setting.kernel,'SE[SE]') == 1
        PK = ETA^2*K_XX.^3.*(-DELTA*(1-exp(-GAMMA^2*(XX-XX').^2)));
        PDEV_DELTA = 0.5*sum(sum(PK.*ALPHA,2));
        
        PK = ETA^2*K_XX.^3.*(-GAMMA*DELTA^2.*(XX-XX').^2.*exp(-GAMMA^2*(XX-XX').^2));
    end
    PDEV(2) = 0.5*sum(sum(PK.*ALPHA,2));
    
    % LAMBDA:
    if strcmp(setting.mode,'graph') == 0
        PK = 2*LAMBDA*eye(NN0);
        PDEV(3) = 0.5*sum(sum(PK.*ALPHA,2));
    end
    
    Mw = beta1*Mw - (1-beta1)*PDEV;
    Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
    
    ETA = abs(ETA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon));
    GAMMA = abs(GAMMA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon));
    LAMBDA = abs(LAMBDA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(3)./(sqrt(Vw(3))+epsilon));
    
    if strcmp(setting.kernel,'SE[SE]') == 1
        Mw_DELTA = beta1*Mw_DELTA - (1-beta1)*PDEV_DELTA;
        Vw_DELTA = beta2*Vw_DELTA + (1-beta2)*PDEV_DELTA.*PDEV_DELTA;
        
        DELTA = abs(DELTA - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw_DELTA./(sqrt(Vw_DELTA)+epsilon));
    end
end


MARGIN = max(X) - min(X);
A = linspace(min(X)-setting.margin*MARGIN,max(X)+setting.margin*MARGIN,1000)';

if strcmp(setting.mode,'graph') == 1 || strcmp(setting.variance,'homoscedastic') == 1
    RR = LAMBDA^2*ones(N,1);
elseif strcmp(setting.variance,'heteroscedastic') == 1
    RR = R;
end

if strcmp(setting.mode,'graph') == 1 || strcmp(setting.variance,'homoscedastic') == 1
    R = LAMBDA^2*ones(size(A));
elseif strcmp(setting.variance,'heteroscedastic') == 1
    
    LL = - 0.5*((A-X')./BW(:,K)').^2 - log(BW(:,K)');
    LL = LL(:,Z==0);
    LL = LL - max(LL,[],2);
    LL = exp(LL);
    
    R = sum(((Y(Z==0)'-MU(Z==0)').^2+NU(Z==0)').*LL,2)./sum(LL,2);
end

if strcmp(setting.kernel,'SE[SE]') == 1
    [MU,NU] = getReg(A,X,Y,Z,ETA,[GAMMA;DELTA],RR,setting);
else
    [MU,NU] = getReg(A,X,Y,Z,ETA,GAMMA,RR,setting);
end

NU = NU + R;

Profile(:,2) = MU*stdv_Y + mean_Y;
Profile(:,3) = sqrt(NU)*stdv_Y;


end