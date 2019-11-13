function [param] = learnProfileParam(reg,data,param,dataType)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;

max_iters = 1000;
N = size(data.Z,1);
M = min(100,N);

SCH = ceil(size(reg.H,2)*rand(max_iters,1));

Mw = zeros(3,1);
Vw = zeros(3,1);

r = 0;
while r < max_iters
    r = r + 1;
    old_param = param;
    
    m = SCH(r);
    
    index = datasample(1:N,M,'Replace',false)';
    HH = reg.H(index,m);
    index = index(HH==0);
    
    ZZ = data.Z(index,m);
    YY = data.Y(index);
    RR = reg.lambda(index,m);
    
    PDEV = getPDEV(YY,ZZ,RR,old_param);
    
    Mw = beta1*Mw - (1-beta1)*PDEV;
    Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
 
    param.eta = abs(old_param.eta - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon));
    param.xi = abs(old_param.xi - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon));
    if strcmp(dataType,'graph') == 0
        param.lambda = abs(old_param.lambda - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(3)./(sqrt(Vw(3))+epsilon));
    end
end


end