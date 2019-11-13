function [data] = getAlignment(data,Profile,setting)

L = length(data);

% Particle Smoothing:
parfor ll = 1:L
    data(ll) = ParticleSmoothing(data(ll),Profile,setting);
end

% Metropolis-Hastings:
parfor ll = 1:L
    data(ll) = MCMC(data(ll),Profile,setting);
end

%{
% Updating R:
parfor ll = 1:L
    data(ll) = getR(data(ll));
end
%}


end