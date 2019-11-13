function [param,MU,NU,QQ] = learnKernelParam(data0,param,Profile,setting,dataType)

data1 = struct('X0',cell(1,1),'Y0',cell(1,1),'X',cell(1,1),'Y',cell(1,1),'mu_X',cell(1,1),'mu_Y',cell(1,1),'sig_X',cell(1,1),'sig_Y',cell(1,1));
data1.X0 = data0.Z(:,1);
data1.Y0 = data0.Y;
data1 = regularizeData(data1);

data = struct('Z0',cell(1,1),'Y0',cell(1,1),'Z',cell(1,1),'Y',cell(1,1),'mu_Z',cell(1,1),'mu_Y',cell(1,1),'sig_Z',cell(1,1),'sig_Y',cell(1,1));
data.Z0 = data0.Z;
data.Y0 = data0.Y;
data.Z = (data0.Z-data1.mu_X)/data1.sig_X;
data.Y = data1.Y;
data.mu_Z = data1.mu_X;
data.mu_Y = data1.mu_Y;
data.sig_Z = data1.sig_X;
data.sig_Y = data1.sig_Y;

[N,M] = size(data.Z);

reg = struct('mu',cell(1,1),'nu',cell(1,1),'lambda',cell(1,1),'H',cell(1,1));
PARAM = struct('eta',cell(1,1),'xi',cell(1,1),'lambda',cell(1,1),'lambda0',cell(1,1),'bandwidth',cell(1,1),'k',cell(1,1));
PARAM.eta = param.eta;
PARAM.xi = param.xi;
PARAM.lambda = param.lambda;
PARAM.k = 1;

reg.H = zeros(N,M);
if strcmp(dataType,'graph') == 1
    reg.lambda = zeros(N,M);
elseif strcmp(dataType,'homoscedastic') == 1
    reg.lambda = zeros(N,M);
elseif strcmp(dataType,'heteroscedastic') == 1
    reg.lambda = 0.3*max(data.Y.^2)*ones(N,M);
end
PARAM.lambda0 = 1e-4;

reg = getReg(reg,data,PARAM);
if strcmp(dataType,'heteroscedastic') == 1
    [reg,PARAM] = getLambda(reg,data,PARAM);
end

for r = 1:25
    if r > 10
        reg = getProfileOutliers(reg,data,setting);
    end
    if strcmp(dataType,'heteroscedastic') == 1
        [reg,PARAM] = getLambda(reg,data,PARAM);
    end
    PARAM = learnProfileParam(reg,data,PARAM,dataType);
    reg = getReg(reg,data,PARAM);
end

param.eta = PARAM.eta;
param.xi = PARAM.xi;
param.lambda = PARAM.lambda;


Z = Profile(:,1);
Z = (Z-data.mu_Z)/data.sig_Z;
MU = zeros(length(Z),M);
NU = zeros(length(Z),M);

QQ = zeros(length(Z),1);

for m = 1:M
    index = (reg.H(:,m)==0);
    
    ZZ = data.Z(index,m);
    YY = data.Y(index);
    RR = reg.lambda(index,m);
    
    N = length(YY);
    
    K = getCov(ZZ,ZZ,param.eta,param.xi) + param.lambda^2*eye(N) + PARAM.lambda0^2*eye(N);
    
    B0 = diag(RR) + K;
    B2 = B0\eye(N);
    B1 = B2*YY;
    
    K1 = getCov(Z,ZZ,param.eta,param.xi);
    K2 = param.eta^2 + param.lambda^2 + PARAM.lambda0^2;
    
    mu = K1*B1;
    nu = K2 - sum((K1*B2).*K1,2);
    
    LA = interp1(data.Z(:,m),reg.lambda(:,m),Z);
    LA(Z<data.Z(1,m)) = reg.lambda(1,m);
    LA(Z>data.Z(end,m)) = reg.lambda(end,m);
    
    MU(:,m) = mu;
    NU(:,m) = nu + LA;
    
    index = (Z>=ZZ(1))&(Z<=ZZ(end));
    QQ(index) = 1;
end

AA = sum(MU,2)/M;
NU = sum(NU+(MU-AA).^2,2)/M;
MU = AA;

MU = MU*data.sig_Y + data.mu_Y;
NU = NU*data.sig_Y^2;


end