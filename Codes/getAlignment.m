function [Samples] = getAlignment(data,Samples,Profile,setting)

L = length(data);

ZZ = cell(L,1);
for ll = 1:L
    ZZ{ll} = Samples.Z(Samples.ID==ll,:);
end

for ll = 1:L
    ZZ{ll} = PS(ZZ{ll},data(ll),Profile,setting);
end

for ll = 1:L
    ZZ{ll} = MCMC(ZZ{ll},data(ll),Profile,setting);
end

for ll = 1:L
    Samples.Z(Samples.ID==ll,:) = ZZ{ll};
    Samples.Y(Samples.ID==ll) = data(ll).Y - data(ll).h;
end


end