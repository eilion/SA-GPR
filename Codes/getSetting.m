function [setting] = getSetting(inputFile)

setting = struct('nIters',cell(1,1));

path = 'Defaults/setting.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

setting.mode = INFO{2}{strcmp(INFO{1},'mode:')==1};
setting.nIters = str2double(INFO{2}{strcmp(INFO{1},'number_of_Iterations:')==1});

setting.kernel = INFO{2}{strcmp(INFO{1},'kernel:')==1};
setting.variance = INFO{2}{strcmp(INFO{1},'variance:')==1};

setting.p = str2double(INFO{2}{strcmp(INFO{1},'p:')==1});
setting.q = str2double(INFO{2}{strcmp(INFO{1},'q:')==1});
setting.r = str2double(INFO{2}{strcmp(INFO{1},'r:')==1});
setting.s = str2double(INFO{2}{strcmp(INFO{1},'s:')==1});

setting.h = str2double(INFO{2}{strcmp(INFO{1},'h:')==1});

setting.d = str2double(INFO{2}{strcmp(INFO{1},'d:')==1});
setting.delta = str2double(INFO{2}{strcmp(INFO{1},'prob_of_outliers:')==1});

setting.nParticles = str2double(INFO{2}{strcmp(INFO{1},'nParticles:')==1});
setting.MH_Iters = str2double(INFO{2}{strcmp(INFO{1},'MH_Iters:')==1});
setting.nPaths = str2double(INFO{2}{strcmp(INFO{1},'nPaths:')==1});

setting.isLearn_shift = INFO{2}{strcmp(INFO{1},'isLearn_shift:')==1};
setting.isLearn_gamma = INFO{2}{strcmp(INFO{1},'isLearn_gamma:')==1};
setting.isLearn_R = INFO{2}{strcmp(INFO{1},'isLearn_R:')==1};

setting.max_AccRate = str2double(INFO{2}{strcmp(INFO{1},'max_AccRate:')==1});
setting.min_AccRate = str2double(INFO{2}{strcmp(INFO{1},'min_AccRate:')==1});

setting.margin = str2double(INFO{2}{strcmp(INFO{1},'margin:')==1});


path = ['Inputs/',inputFile,'/setting.txt'];
if exist(path,'file') == 2
    fileID = fopen(path);
    INFO = textscan(fileID,'%s %s');
    fclose(fileID);
    
    if sum(strcmp(INFO{1},'mode:')==1) == 1
        setting.mode = INFO{2}{strcmp(INFO{1},'mode:')==1};
    end
    
    if sum(strcmp(INFO{1},'number_of_Iterations:')==1) == 1
        setting.nIters = str2double(INFO{2}{strcmp(INFO{1},'number_of_Iterations:')==1});
    end
    
    if sum(strcmp(INFO{1},'kernel:')==1) == 1
        setting.kernel = INFO{2}{strcmp(INFO{1},'kernel:')==1};
    end
    
    if sum(strcmp(INFO{1},'variance:')==1) == 1
        setting.variance = INFO{2}{strcmp(INFO{1},'variance:')==1};
    end
    
    if sum(strcmp(INFO{1},'p:')==1) == 1
        setting.p = str2double(INFO{2}{strcmp(INFO{1},'p:')==1});
    end
    
    if sum(strcmp(INFO{1},'q:')==1) == 1
        setting.q = str2double(INFO{2}{strcmp(INFO{1},'q:')==1});
    end
    
    if sum(strcmp(INFO{1},'r:')==1) == 1
        setting.r = str2double(INFO{2}{strcmp(INFO{1},'r:')==1});
    end
    
    if sum(strcmp(INFO{1},'s:')==1) == 1
        setting.s = str2double(INFO{2}{strcmp(INFO{1},'s:')==1});
    end
    
    if sum(strcmp(INFO{1},'h:')==1) == 1
        setting.h = str2double(INFO{2}{strcmp(INFO{1},'h:')==1});
    end
    
    if sum(strcmp(INFO{1},'d:')==1) == 1
        setting.d = str2double(INFO{2}{strcmp(INFO{1},'d:')==1});
    end
    
    if sum(strcmp(INFO{1},'prob_of_outliers:')==1) == 1
        setting.delta = str2double(INFO{2}{strcmp(INFO{1},'prob_of_outliers:')==1});
    end
    
    if sum(strcmp(INFO{1},'nParticles:')==1) == 1
        setting.nParticles = str2double(INFO{2}{strcmp(INFO{1},'nParticles:')==1});
    end
    
    if sum(strcmp(INFO{1},'MH_Iters:')==1) == 1
        setting.MH_Iters = str2double(INFO{2}{strcmp(INFO{1},'MH_Iters:')==1});
    end
    
    if sum(strcmp(INFO{1},'nPaths:')==1) == 1
        setting.nPaths = str2double(INFO{2}{strcmp(INFO{1},'nPaths:')==1});
    end
    
    if sum(strcmp(INFO{1},'isLearn_shift:')==1) == 1
        setting.isLearn_shift = INFO{2}{strcmp(INFO{1},'isLearn_shift:')==1};
    end
    
    if sum(strcmp(INFO{1},'isLearn_gamma:')==1) == 1
        setting.isLearn_gamma = INFO{2}{strcmp(INFO{1},'isLearn_gamma:')==1};
    end
    
    if sum(strcmp(INFO{1},'isLearn_R:')==1) == 1
        setting.isLearn_R = INFO{2}{strcmp(INFO{1},'isLearn_R:')==1};
    end
    
    if sum(strcmp(INFO{1},'max_AccRate:')==1) == 1
        setting.max_AccRate = str2double(INFO{2}{strcmp(INFO{1},'max_AccRate:')==1});
    end
    
    if sum(strcmp(INFO{1},'min_AccRate:')==1) == 1
        setting.min_AccRate = str2double(INFO{2}{strcmp(INFO{1},'min_AccRate:')==1});
    end
    
    if sum(strcmp(INFO{1},'margin:')==1) == 1
        setting.margin = str2double(INFO{2}{strcmp(INFO{1},'margin:')==1});
    end
end


end