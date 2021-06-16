function saveFigure(inputFile)

path = ['Outputs/',inputFile,'/results.mat'];
results = load(path);
results = results.results;

data = results.data;
Samples = results.Samples;
Profile = results.Profile;
setting = results.setting;

L = size(data,1);


% Raw Signals:
for k = 1:ceil(L/4)
    fig = figure;
    for ll = (k-1)*4+1:min(4*k,L)
        subplot(2,2,ll-(k-1)*4);
        hold on;
        title(['Signal ',num2str(ll)],'FontSize',16);
        plot(data(ll).X,data(ll).Y,'*k');
    end
    
    set(fig,'Position',[20 20 1000 600]);
    movegui(fig,'center');
    
    path = ['Outputs/',inputFile,'/raw_signals (',num2str(k),').fig'];
    savefig(fig,path);
end

MIN = min(Samples.Y) - 0.05*(max(Samples.Y)-min(Samples.Y));
MAX = max(Samples.Y) + 0.05*(max(Samples.Y)-min(Samples.Y));


% Alignments:
cc = jet(L);
h = zeros(L,1);
fig = figure;
hold on;
title('Alignments','FontSize',16);
xx = Profile(:,1);
mu = Profile(:,2);
sig = Profile(:,3);
aa = [xx;flipud(xx)];
bb = [mu-1.96*sig;flipud(mu+1.96*sig)];
patch(aa,bb,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
for ll = 1:L
    h(ll) = plot(median(Samples.Z(Samples.ID==ll,:),2),Samples.Y(Samples.ID==ll),'*','Color',cc(ll,:));
end
xlim([-1 1]);
ylim([MIN MAX]);
LIST = cell(L,1);
for ll = 1:L
    LIST{ll} = data(ll).name;
    LIST{ll}(LIST{ll}=='_') = '-';
end
legend(h,LIST,'Location','SouthOutside');

set(fig,'Position',[20 20 1000 600]);
movegui(fig,'center');

path = ['Outputs/',inputFile,'/alignments.fig'];
savefig(fig,path);

for k = 1:ceil(L/4)
    fig = figure;
    for ll = (k-1)*4+1:min(4*k,L)
        subplot(2,2,ll-(k-1)*4);
        hold on;
        title(['Signal ',num2str(ll)],'FontSize',16);
        
        xx = Profile(:,1);
        mu = Profile(:,2);
        sig = Profile(:,3);
        aa = [xx;flipud(xx)];
        bb = [mu-1.96*sig;flipud(mu+1.96*sig)];
        patch(aa,bb,1,'FaceColor','k','FaceAlpha',0.1,'EdgeColor','none');
        
        plot(median(Samples.Z(Samples.ID==ll,:),2),Samples.Y(Samples.ID==ll),'*','Color',cc(ll,:));
        
        xlim([-1 1]);
        ylim([MIN MAX]);
    end
    
    set(fig,'Position',[20 20 1000 600]);
    movegui(fig,'center');
    
    path = ['Outputs/',inputFile,'/alignments (',num2str(k),').fig'];
    savefig(fig,path);
end


% Accumulation Rates:
for k = 1:ceil(L/4)
    fig = figure;
    for ll = (k-1)*4+1:min(4*k,L)
        subplot(2,2,ll-(k-1)*4);
        hold on;
        title(['Signal ',num2str(ll)],'FontSize',16);
        
        h = zeros(2,1);
        
        ZZ = Samples.Z(Samples.ID==ll,:);
        XX = (data(ll).X-setting.X_mean(ll))/setting.X_stdv(ll);
        TT = (ZZ(2:end,:)-ZZ(1:end-1,:))./(XX(2:end)-XX(1:end-1));
        TT = TT/data(ll).R;
        
        bins = 0:0.25:10;
        hh = histogram(TT,bins,'FaceColor','c');
        hh.Normalization = 'pdf';
        
        xx = 0.01:0.01:10;
        yy = gampdf(xx,data(ll).alpha,1/data(ll).beta);
        plot(xx,yy,'g','LineWidth',2);
        
        h(1) = plot(inf,inf,'Color',[1,1,1]);
        h(2) = plot(inf,inf,'Color',[1,1,1]);
        
        legend(h,{['shape = ',num2str(data(ll).alpha)],['rate = ',num2str(data(ll).beta)]},'Location','NorthEast');
        
        xlim([0 6]);
    end
    
    set(fig,'Position',[20 20 1000 600]);
    movegui(fig,'center');
    
    path = ['Outputs/',inputFile,'/accumulation_rates (',num2str(k),').fig'];
    savefig(fig,path);
end



end