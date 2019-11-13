function [setting] = getSetting()

setting = struct('q',cell(1,1),'d',cell(1,1),'h',cell(1,1),'nParticles',cell(1,1),'alpha',cell(1,1),'beta',cell(1,1),'MH_Iters',cell(1,1),'nPaths',cell(1,1),'margin',cell(1,1));

path = 'Defaults/setting.txt';
fileID = fopen(path);
INFO = textscan(fileID,'%s %s');
fclose(fileID);

setting.nu = str2double(INFO{2}{strcmp(INFO{1},'nu:')==1});
setting.q = str2double(INFO{2}{strcmp(INFO{1},'q:')==1});
setting.d = str2double(INFO{2}{strcmp(INFO{1},'d:')==1});
setting.h = str2double(INFO{2}{strcmp(INFO{1},'h:')==1});
setting.nParticles = str2double(INFO{2}{strcmp(INFO{1},'nParticles:')==1});
setting.alpha = str2double(INFO{2}{strcmp(INFO{1},'alpha:')==1});
setting.beta = str2double(INFO{2}{strcmp(INFO{1},'beta:')==1});
setting.MH_Iters = str2double(INFO{2}{strcmp(INFO{1},'MH_Iters:')==1});
setting.nPaths = str2double(INFO{2}{strcmp(INFO{1},'nPaths:')==1});
setting.margin = str2double(INFO{2}{strcmp(INFO{1},'margin:')==1});


end