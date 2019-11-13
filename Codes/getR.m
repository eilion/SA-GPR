function [data] = getR(data)

Z = data.Z;
X = data.X;

M = size(Z,2);

A = data.R(:,1);

T = Z(1:end-1,:);
Y = (Z(2:end,:)-Z(1:end-1,:))./(X(2:end)-X(1:end-1));

T = T(:);
Y = Y(:);

Dist = abs(A-T');
Dist = sort(Dist,2,'ascend');

BW = Dist(:,max(M,2));

Dist2 = (A-T').^2;
WW = - 0.5*Dist2./(BW.^2) - log(BW);
WW = WW - max(WW,[],2);
WW = exp(WW);
LL = WW./sum(WW,2);

data.R(:,2) = sum(Y'.*LL,2);


end