x = linspace(-5,10,200);
y = cos(pi*x/2).*(x.^2/10);

index2 = 1:100;
index1 = 51:150;
index3 = 101:200;

figure;
hold on;
plot(x(index1),y(index1),'o','LineWidth',2);
plot(x(index2),y(index2),'.','LineWidth',2);
plot(x(index3),y(index3),'.','LineWidth',2);

QQ1 = [x(index1)',y(index1)'];
QQ2 = [x(index2)',y(index2)'];
QQ3 = [x(index3)',y(index3)'];

[~,order] = sort(QQ1(:,1),'ascend');
QQ1 = QQ1(order,:);

[~,order] = sort(QQ2(:,1),'ascend');
QQ2 = QQ2(order,:);

[~,order] = sort(QQ3(:,1),'ascend');
QQ3 = QQ3(order,:);