function [results] = SAGPR(inputFile,dataType)

setting = getSetting();
data = getData(inputFile);
param = getParam(data,dataType);

% Initialization:
disp('## Profile initialization is running...');
[data,Profile] = getInitialProfile(data,dataType,setting);
disp('#  Done.');

% Iteration:
disp('## Alignment initialization is running...');
data = getInitialAlignment(data,Profile,setting);
disp('#  Done.');
for r = 1:5
    disp(['## Iteration ',num2str(r),':']);
    % Alignment:
    disp('#  Alignment algorithm is running...');
    data = getAlignment(data,Profile,setting);
    disp('   Done.');
    % Profile Construction:
    disp('#  Profile is updated...');
    [param,Profile] = getUpdatedProfile(data,param,Profile,setting,dataType);
    disp('   Done.');
end

results = getResults(data,param,Profile);

% Storing the results:
path = ['Outputs/',inputFile];
if exist(path,'dir') == 7
    n = 0;
    DET = 1;
    while DET == 1
        n = n + 1;
        path = ['Outputs/',inputFile,'(',num2str(n),')'];
        if exist(path,'dir') == 0
            DET = 0;
        end
    end
end
mkdir(path);

fileID = [path,'/results.mat'];
save(fileID,'results');


end