function param = learnParam(reg,data,param,dataType)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma = 1e-3;

max_iters = 1000;
N = size(data.X,1);
M = min(200,N);

Mw = zeros(3,1);
Vw = zeros(3,1);

r = 0;
while r < max_iters
    r = r + 1;
    old_param = param;
    
    index = datasample(1:N,M,'Replace',false)';
    HH = reg.H(index);
    index = index(HH==0);
    
    XX = data.X(index);
    YY = data.Y(index);
    RR = reg.lambda(index);
    
    PDEV = getPDEV(YY,XX,RR,old_param);
    
    Mw = beta1*Mw - (1-beta1)*PDEV;
    Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
 
    param.eta = abs(old_param.eta - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon));
    param.xi = abs(old_param.xi - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon));
    if strcmp(dataType,'graph') == 0
        param.lambda = abs(old_param.lambda - (gamma*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(3)./(sqrt(Vw(3))+epsilon));
    end
end


end