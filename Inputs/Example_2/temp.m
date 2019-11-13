x = linspace(-10,10,10000)';

X = cell(4,1);
Y = cell(4,1);

for ll = 1:4
    X{ll} = ceil(10000*rand(220,1));
    X{ll} = sort(X{ll},'ascend');
    X{ll} = unique(X{ll});
    X{ll} = X{ll}(1:200);
    X{ll} = x(X{ll})';
end

mu = (x/10+1).^2.*(cos(x/2)+sin(x+1))/2;

figure;
hold on;
plot(x,mu);