function [data] = getInitialAlignment(data,Profile,setting)

L = length(data);

shortProfile = zeros(200,3);
shortProfile(:,1) = linspace(Profile(1,1),Profile(end,1),200);
shortProfile(:,2) = interp1(Profile(:,1),Profile(:,2),shortProfile(:,1));
shortProfile(:,3) = interp1(Profile(:,1),Profile(:,3),shortProfile(:,1));

parfor ll = 1:L
    data(ll) = shortHMM(data(ll),shortProfile,setting);
end


end