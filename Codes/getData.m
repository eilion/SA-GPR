function [data] = getData(inputFile)

list = dir(['Inputs/',inputFile,'/*.txt']);

L = length(list);

data = struct('name',cell(L,1),'X',cell(L,1),'Y',cell(L,1),'Z',cell(L,1),'R',cell(L,1));
for ll = 1:L
    data(ll).name = list(ll).name(1:end-4);
    
    path = ['Inputs/',inputFile,'/',list(ll).name];
    AA = load(path);
    data(ll).X = AA(:,1);
    data(ll).Y = AA(:,2);
end



end