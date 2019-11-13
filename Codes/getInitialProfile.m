function [data0,Profile] = getInitialProfile(data0,dataType,setting)

data = struct('X0',cell(1,1),'Y0',cell(1,1),'X',cell(1,1),'Y',cell(1,1),'mu_X',cell(1,1),'mu_Y',cell(1,1),'sig_X',cell(1,1),'sig_Y',cell(1,1));
data.X0 = data0(1).X;
data.Y0 = data0(1).Y;
data = regularizeData(data);

reg = struct('mu',cell(1,1),'nu',cell(1,1),'lambda',cell(1,1),'H',cell(1,1));
param = struct('eta',cell(1,1),'xi',cell(1,1),'lambda',cell(1,1),'lambda0',cell(1,1),'bandwidth',cell(1,1),'k',cell(1,1));

param.k = 5;

N = size(data.X,1);
reg.H = zeros(N,1);
if strcmp(dataType,'graph') == 1
    reg.lambda = zeros(N,1);
    param.eta = 1;
    param.xi = 10;
    param.lambda = 0.1;
elseif strcmp(dataType,'homoscedastic') == 1
    reg.lambda = zeros(N,1);
    param.eta = 1;
    param.xi = 1;
    param.lambda = 0.01;
elseif strcmp(dataType,'heteroscedastic') == 1
    reg.lambda = 0.3*max(data.Y.^2)*ones(N,1);
    param.eta = 1;
    param.xi = 10;
    param.lambda = 0.01;
end
param.lambda0 = 1e-4;

reg = getInitialReg(reg,data,param);
if strcmp(dataType,'heteroscedastic') == 1
    param = getBandwidth(reg,data,param);
    reg = learnVar_NW(reg,data,param);
end

for r = 1:25
    if r > 10
        reg = getOutliers(reg,data);
    end
    if strcmp(dataType,'heteroscedastic') == 1
        param = getBandwidth(reg,data,param);
        reg = learnVar_NW(reg,data,param);
    end
    param = learnParam(reg,data,param,dataType);
    reg = getInitialReg(reg,data,param);
end

Profile = getProfile(reg,data,param,setting);
Profile(:,1) = Profile(:,1)*data.sig_X + data.mu_X;
Profile(:,2) = Profile(:,2)*data.sig_Y + data.mu_Y;
Profile(:,3) = Profile(:,3)*data.sig_Y;

for ll = 1:length(data0)
    N = size(Profile,1);
    data0(ll).R = zeros(N,2);
    data0(ll).R(:,1) = Profile(:,1);
    % data0(ll).R(:,2) = 2*ones(N,1)*(data0(1).X(end)-data0(1).X(1))/(data0(ll).X(end)-data0(ll).X(1));
    data0(ll).R(:,2) = ones(N,1)*(data0(1).X(end)-data0(1).X(1))/(data0(ll).X(end)-data0(ll).X(1));
end


end