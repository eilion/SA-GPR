x = linspace(-10,10,100);
y = cos(pi*x/18*7);

index1 = 1:2:length(x);
index2 = 2:2:length(x);

figure;
hold on;
plot(x(index1),y(index1),'*','LineWidth',2);
plot(x(index2),y(index2),'*','LineWidth',2);

QQ1 = [x(index1)',y(index1)'];
QQ2 = [x(index2)',y(index2)'];

[~,order] = sort(QQ1(:,1),'ascend');
QQ1 = QQ1(order,:);

[~,order] = sort(QQ2(:,1),'ascend');
QQ2 = QQ2(order,:);