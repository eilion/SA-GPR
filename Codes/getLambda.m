function [reg,param] = getLambda(reg,data,param)

[N,M] = size(reg.H);
K = ceil(param.k*N/100);

for m = 1:M
    Z = data.Z(:,m);
    
    index = (reg.H(:,m)==0);
    ZZ = Z(index);

    Dist = abs(Z-ZZ');

    AA = sort(Dist,2,'ascend');
    param.bandwidth(:,m) = AA(:,min(K,sum(index)));
end

for m = 1:M
    Z = data.Z(:,m);
    Y = data.Y;
    
    BW = param.bandwidth(:,m);
    
    index = (reg.H(:,m)==0);
    ZZ = Z(index,:);
    
    Dist2 = (Z-ZZ').^2;
    WW = - 0.5*Dist2./(BW.^2) - log(BW);
    WW = WW - max(WW,[],2);
    WW = exp(WW);
    LL = WW./sum(WW,2);
    
    Y = Y(index)';
    mu = reg.mu(index,m)';
    nu = reg.nu(index,m)';
    
    reg.lambda(:,m) = sum(((Y-mu).^2+nu).*LL,2);
end


end