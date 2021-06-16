function [data,setting] = getData(inputFile,setting)

list = dir(['Inputs/',inputFile,'/Signals/*.txt']);

L = length(list);

setting.X_mean = zeros(1,L);
setting.X_stdv = zeros(1,L);

setting.Y_mean = zeros(1,L);
setting.Y_stdv = zeros(1,L);

data = struct('name',cell(L,1),'X',cell(L,1),'Y',cell(L,1),'Z',cell(L,1),'h',cell(L,1),'R',cell(L,1),'alpha',cell(L,1),'beta',cell(L,1));
for ll = 1:L
    data(ll).name = list(ll).name(1:end-4);
    
    path = ['Inputs/',inputFile,'/Signals/',list(ll).name];
    AA = load(path);
    data(ll).X = AA(:,1);
    data(ll).Y = AA(:,2);
    data(ll).h = 0;
    data(ll).R = 1;
    data(ll).alpha = 2;
    data(ll).beta = 1;
    
    setting.X_mean(ll) = mean(data(ll).X);
    setting.X_stdv(ll) = sqrt(var(data(ll).X));
    
    setting.Y_mean(ll) = mean(data(ll).Y);
    setting.Y_stdv(ll) = 1;
    
    data(ll).X = (data(ll).X-setting.X_mean(ll))/setting.X_stdv(ll);
    data(ll).Y = (data(ll).Y-setting.Y_mean(ll))/setting.Y_stdv(ll);
end


end