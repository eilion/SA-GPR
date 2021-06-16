function [inputFile] = saveResults(inputFile,results)

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

inputFile = path(9:end);


end