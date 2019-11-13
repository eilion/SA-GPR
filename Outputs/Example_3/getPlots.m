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

INT = 3;
fig = figure;
subplot(2,3,[1;4]);
hold on;
for ll = 1:L
    plot(Data(ll).inputs,INT*(ll-1)+Data(ll).outputs,'.');
end
xlim([-12.5,11]);
ylim([-2,12]);
xlabel('(a)','FontSize',12);
subplot(2,3,[2,3]);
hold on;
xx = [Profile.latents;flipud(Profile.latents)];
yy = [Profile.means-1.96*Profile.stdvs;flipud(Profile.means+1.96*Profile.stdvs)];
patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
xlabel('(b)','FontSize',12);
xlim([-12.5,11]);
ylim([-1,2]);
subplot(2,3,[5,6]);
hold on;
xx = [Profile.latents;flipud(Profile.latents)];
yy = [Profile.means-1.96*Profile.stdvs;flipud(Profile.means+1.96*Profile.stdvs)];
patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
for ll = 1:L
    plot(median(Samples(ll).latents,2),Data(ll).outputs,'.');
end
xlabel('(c)','FontSize',12);
xlim([-12.5,11]);
ylim([-1,2]);


set(fig,'Position',[20 20 900 400]);
movegui(fig,'center');