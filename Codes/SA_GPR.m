function [results] = SA_GPR(inputFile)

% load settings:
setting = getSetting(inputFile);

% load data:
[data,setting] = getData(inputFile,setting);

% initialize the profile:
Profile = initializeProfile(data,setting);

% main iteration:
for r = 1:setting.nIters
    Samples = InitAlignment(data,Profile,setting);
    for t = 1:5
        Samples = getAlignment(data,Samples,Profile,setting);
        data = learnParam(data,Samples,Profile,setting);
    end
    if r < setting.nIters
        Profile = learnProfile(Samples,setting);
    end
end
for ll = 1:length(data)
    data(ll).Z = Samples.Z(Samples.ID==ll,:);
    data(ll).X = data(ll).X*setting.X_stdv(ll) + setting.X_mean(ll);
    data(ll).Y = data(ll).Y*setting.Y_stdv(ll) + setting.Y_mean(ll);
    data(ll).h = data(ll).h*setting.Y_stdv(ll);
end

% store results:
results = struct('data',cell(1,1),'Samples',cell(1,1),'Profile',cell(1,1),'setting',cell(1,1));
results.data = data;
results.Samples = Samples;
results.Profile = Profile;
results.setting = setting;

inputFile = saveResults(inputFile,results);
saveFigure(inputFile);


end