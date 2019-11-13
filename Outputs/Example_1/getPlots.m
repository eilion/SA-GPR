AA = load('results.mat');
results = AA.results;
clear AA;

Data = results.data;
Samples = results.samples;
Profile = results.profile;

[~,id1,id2] = dtw(Data(1).outputs,Data(2).outputs);
AL_DTW = [Data(1).inputs(id1),Data(2).inputs(id2)];

ZZ = [median(Samples(1).latents,2),median(Samples(2).latents,2)];
AL_GPR = [interp1(ZZ(:,1),Data(1).inputs,ZZ(:,2)),Data(2).inputs];

L = length(Data);

fig = figure;
subplot(3,1,1);
hold on;
x = linspace(-10.5,10.5,1000);
y = cos(pi*x/18*7);
plot(x,y,'k');
plot(Data(1).inputs,Data(1).outputs,'*r');
plot(Data(2).inputs,Data(2).outputs,'*b');
plot([-10.5,10.5],[-1.2,-1.2],'k','LineWidth',1);
plot([-10.5,10.5],[1.2,1.2],'k','LineWidth',1);
plot(Data(1).inputs,-1.2*ones(1,length(Data(1).inputs)),'^r');
plot(Data(2).inputs,1.2*ones(1,length(Data(2).inputs)),'vb');
for m = 1:length(Data(1).inputs)
    plot([Data(1).inputs(m),Data(1).inputs(m)],[Data(1).outputs(m),-1.2],':r');
end
for m = 1:length(Data(2).inputs)
    plot([Data(2).inputs(m),Data(2).inputs(m)],[Data(2).outputs(m),1.2],':b');
end
xlim([-10.5,10.5]);
ylim([-1.2,1.2]);
xlabel('(a)','FontSize',12);

subplot(3,1,[2;3]);
hold on;
plot(x,y-1.5,'k');
plot(x,y+1.5,'k');
plot(Data(1).inputs,Data(1).outputs-1.5,'*r');
plot(Data(2).inputs,Data(2).outputs+1.5,'*b');
plot([-10.5,10.5],[-3.2,-3.2],'k','LineWidth',1);
plot([-10.5,10.5],[3.2,3.2],'k','LineWidth',1);
plot(Data(1).inputs,-3.2*ones(1,length(Data(1).inputs)),'^r');
plot(Data(2).inputs,3.2*ones(1,length(Data(2).inputs)),'vb');
for m = 1:length(Data(1).inputs)
    plot([Data(1).inputs(m),Data(1).inputs(m)],[Data(1).outputs(m)-1.5,-3.2],':r');
end
for m = 1:length(Data(2).inputs)
    plot([Data(2).inputs(m),Data(2).inputs(m)],[Data(2).outputs(m)+1.5,3.2],':b');
end
xlim([-10.5,10.5]);
ylim([-3.2,3.2]);
xlabel('(b)','FontSize',12);
set(fig,'Position',[20 20 900 900]);
movegui(fig,'center');

fig = figure;
subplot(3,1,1);
hold on;
xx = [Profile.latents;flipud(Profile.latents)];
yy = [Profile.means-1.96*Profile.stdvs;flipud(Profile.means+1.96*Profile.stdvs)];
patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
plot(median(Samples(1).latents,2),Data(1).outputs,'*-r');
plot(median(Samples(2).latents,2),Data(2).outputs,'*-b');
xlabel('(a)','FontSize',12);
xlim([-10.5,10.5]);
ylim([-1.2,1.2]);

subplot(3,1,2);
hold on;
plot(Data(1).inputs(id1),Data(1).outputs(id1),'*-r');
plot(Data(1).inputs(id1),Data(2).outputs(id2),'*-b');
xlabel('(b)','FontSize',12);
xlim([-10.5,10.5]);
ylim([-1.2,1.2]);

subplot(3,1,3);
hold on;
plot([0,0],[Data(1).inputs(1),Data(1).inputs(end)],'k','LineWidth',2);
plot(AL_DTW(:,2)-AL_DTW(:,1),AL_DTW(:,1),'r','LineWidth',2);
plot(AL_GPR(:,2)-AL_GPR(:,1),AL_GPR(:,1),'g','LineWidth',2);
xlim([-0.3 0.3]);
ylim([-10.5,10.5]);
xlabel('(c)','FontSize',12);
ylabel('inputs of signal 1','FontSize',12);

set(fig,'Position',[20 20 900 900]);
movegui(fig,'center');