function [results] = getResults(data,param,Profile)

L = length(data);

results = struct('data',cell(1,1),'samples',cell(1,1),'param',cell(1,1),'profile',cell(1,1));

results.data = struct('name',cell(L,1),'inputs',cell(L,1),'outputs',cell(L,1));
results.samples = struct('name',cell(L,1),'latents',cell(L,1));
results.param = struct('name',cell(L,1),'eta',cell(L,1),'xi',cell(L,1),'R',cell(L,1));
results.profile = struct('latents',cell(1,1),'means',cell(1,1),'stdvs',cell(1,1));

for ll = 1:L
    results.data(ll).name = data(ll).name;
    results.data(ll).inputs = data(ll).X;
    results.data(ll).outputs = data(ll).Y;
end

for ll = 1:L
    results.samples(ll).name = data(ll).name;
    results.samples(ll).latents = data(ll).Z;
end

for ll = 1:L
    results.param(ll).name = data(ll).name;
    results.param(ll).eta = param(ll).eta;
    results.param(ll).xi = param(ll).xi;
    results.param(ll).R = data(ll).R;
end

results.profile.latents = Profile(:,1);
results.profile.means = Profile(:,2);
results.profile.stdvs = Profile(:,3);


end