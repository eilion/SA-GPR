function [Samples] = InitAlignment(data,Profile,setting)

Samples = struct('Z',cell(1,1),'ID',cell(1,1),'Y',cell(1,1));
L = length(data);

N = 0;
for ll = 1:L
    N = N + size(data(ll).Y,1);
    data(ll).R = 1;
    data(ll).alpha = 2;
    data(ll).beta = 1;
end

Samples.Y = zeros(N,1);
Samples.ID = zeros(N,1);
Samples.Z = zeros(N,setting.nPaths);

ZZ = cell(L,1);
for ll = 1:L
    ZZ{ll} = initializeAlignment(data(ll),Profile,setting);
end

n = 0;
for ll = 1:L
    N0 = size(data(ll).Y,1);
    Samples.Z(n+1:n+N0,:) = ZZ{ll};
    Samples.Y(n+1:n+N0) = data(ll).Y - data(ll).h;
    Samples.ID(n+1:n+N0) = ll;
    n = n + N0;
end


end