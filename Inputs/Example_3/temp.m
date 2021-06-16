x = linspace(-12,12,10000);

X = cell(4,1);
Y = cell(4,1);

for ll = 1:4
    X{ll} = ceil(10000*rand(220,1));
    X{ll} = sort(X{ll},'ascend');
    X{ll} = unique(X{ll});
    X{ll} = X{ll}(1:200);
    X{ll} = x(X{ll})';
end

x = linspace(-12,12,200);
std = 0.01*ones(1,200);
std_1 = normpdf(x,-6,0.75);
std_2 = normpdf(x,6,1);
std = std + 0.5*std_1 + std_2;

for ll = 1:4
    mu = (X{ll}/10).^3 + (X{ll}/10).^2 - X{ll}/20;
    y = normrnd(0,1,[200,1]).*std' + mu;
    Y{ll} = y;
end

figure;
hold on;
for ll = 1:4
    plot(X{ll},Y{ll}+2*ll);
end

figure;
hold on;
for ll = 1:4
    plot(X{ll},Y{ll},'.');
end

QQ = cell(4,1);
for ll = 1:4
    QQ{ll} = [X{ll},Y{ll}];
end