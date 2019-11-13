function [reg] = getReg(reg,data,param)

M = size(reg.H,2);

for m = 1:M
    index = (reg.H(:,m)==0);
    
    ZZ = data.Z(index,m);
    YY = data.Y(index);
    RR = reg.lambda(index,m);
    
    N = length(YY);
    
    K = getCov(ZZ,ZZ,param.eta,param.xi) + param.lambda^2*eye(N) + param.lambda0^2*eye(N);
    
    B0 = diag(RR) + K;
    B2 = B0\eye(N);
    B1 = B2*YY;
    
    K1 = getCov(data.Z(:,m),ZZ,param.eta,param.xi);
    K2 = param.eta^2 + param.lambda^2 + param.lambda0^2;
    
    reg.mu(:,m) = K1*B1;
    reg.nu(:,m) = K2 - sum((K1*B2).*K1,2);
end


end