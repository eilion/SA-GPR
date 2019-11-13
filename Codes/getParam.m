function [param] = getParam(data,dataType)

L = length(data);

param = struct('name',cell(L,1),'eta',cell(L,1),'xi',cell(L,1),'lambda',cell(L,1),'mu_Z',cell(L,1),'mu_Y',cell(L,1),'sig_Z',cell(L,1),'sig_Y',cell(L,1));

for ll = 1:L
    param(ll).name = data(ll).name;
    if strcmp(dataType,'graph') == 1
        param(ll).eta = 1;
        param(ll).xi = 10;
        param(ll).lambda = 0.1;
    elseif strcmp(dataType,'homoscedastic') == 1
        param(ll).eta = 1;
        param(ll).xi = 10;
        param(ll).lambda = 0.01;
    elseif strcmp(dataType,'heteroscedastic') == 1
        param(ll).eta = 1;
        param(ll).xi = 10;
        param(ll).lambda = 0.01;
    end
end


end

