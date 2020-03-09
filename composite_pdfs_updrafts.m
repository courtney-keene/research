load('updraft_objects_20190227.mat','upzc','upzx','upxc','upZi','upctime','CWIDTH','upspd','upwmax','upwstar','upwmean','uparea','upzc2','upzx2','upxc2','upZi2','CWIDTH2','upspd2','upwmax2','upwstar2','upctime2','upwmean2','uparea2','xideal','zideal','ximat','zimat');
load('updraft_wnorm_20190227.mat','Wnorm','Wnorm2','xideal','zideal');

%% Compute upper 3/4 updraft location
upz_norm=(0.5.*(upzc+upzx))./upZi; %CBL height normalized 3/4 height. 
upzxn=upzx./upZi;
upz_norm2=(0.5.*(upzc2+upzx2))./upZi2;
upzxn2=upzx2./upZi2;
qtiles=[.25 .45 .65 .85 1.05 1.25] %normalized CBL height bins
% qtiles=[0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25] %normalized CBL height bins
% qtiles2=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]

for ii=1:(length(qtiles)-1)
    if ii<length(qtiles)
        idx(ii).locs=find(upz_norm>=qtiles(ii) & upz_norm<qtiles(ii+1)); %indices for updrafts in a given height range (e.g., between .25 and .45 Zi)
        idx2(ii).locs=find(upz_norm2>=qtiles(ii) & upz_norm2<qtiles(ii+1));
    end
end

%%
f=figure(3);clf;set(f,'color','w');
CCHORD=CWIDTH(:).*upspd;%(nidx);
CCHORD2=CWIDTH2(:).*upspd(2);
CCHORD_ALL=cat(1,CCHORD,CCHORD2);
UPWMAX=cat(2,upwmax,upwmax2);
for ii=1:(length(qtiles)-1)
    sp=subplot(3,length(qtiles)-1,length(qtiles)-1+ii);
    cla(sp);
    idxnow=[idx(ii).locs; idx2(ii).locs+53727];
    mchord(ii)=round(nanmean(CCHORD_ALL(idxnow)));
    mdchord(ii)=round(nanmedian(CCHORD_ALL(idxnow)));
    histogram(CCHORD_ALL(idxnow),[0:100:3000],'Normalization','probability','facecolor','b');
    pd = makedist(CCHORD_ALL(idxnow),[0:100:3000])
    y = pdf(pd,[0:100:3000])
    plot([0:100,3000],y)
    ylim([0 .45]);
    xlim([0 2000]);
    grid on; box on;
    if ismember(ii,[2 3 4 5])
        set(gca,'YTickLabel',[]');
    elseif ii==6
        set(gca,'YAxisLocation','right');
        ylabel('Prob.');
    elseif ii==1
        ylabel('Prob.');
    end
    xlabel('Chord Length [m]');
    text(800,.3,strcat('Mean=',num2str(mchord(ii))),'fontsize',17,'fontweight','bold');
    text(800,.25,strcat('Median=',num2str(mdchord(ii))),'fontsize',17,'fontweight','bold');
    set(gca,'fontsize',15,'fontweight','bold','linewidth',2,'layer','top');
    
    sp=subplot(3,length(qtiles)-1,2*(length(qtiles)-1)+ii);
    cla(sp);
    
    WM(ii)=(nanmean(UPWMAX(idxnow)));
    WMD(ii)=(nanmedian(UPWMAX(idxnow)));
    histogram(UPWMAX(idxnow),[0:.5:9],'Normalization','probability','facecolor','b');
    ylim([0 .4]);
    %xlim([0 2500]);
    grid on; box on;
    if ismember(ii,[2 3 4 5])
        set(gca,'YTickLabel',[]');
    elseif ii==6
        set(gca,'YAxisLocation','right');
        ylabel('Prob.');
    elseif ii==1
        ylabel('Prob.');
    end
    xlabel('max updraft [m s^{-1}]');
    text(1,.35,strcat('Mean=',num2str(WM(ii))),'fontsize',17,'fontweight','bold');
    text(1,.3,strcat('Median=',num2str(WMD(ii))),'fontsize',17,'fontweight','bold');
    set(gca,'fontsize',15,'fontweight','bold','linewidth',2,'layer','top');
end


