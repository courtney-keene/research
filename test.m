load('updraft_objects_20190227.mat','upzc','upzx','upxc','upZi','upctime','CWIDTH','upspd','upwmax','upwstar','upwmean','uparea','upzc2','upzx2','upxc2','upZi2','CWIDTH2','upspd2','upwmax2','upwstar2','upctime2','upwmean2','uparea2','xideal','zideal','ximat','zimat');
load('updraft_wnorm_20190227.mat','Wnorm','Wnorm2','xideal','zideal');


%% Compute upper 3/4 updraft location
upz_norm=(0.5.*(upzc+upzx))./upZi; %CBL height normalized 3/4 height. 
upzxn=upzx./upZi;
upz_norm2=(0.5.*(upzc2+upzx2))./upZi2;
upzxn2=upzx2./upZi2;
% qtiles=[.25 .45 .65 .85 1.05 1.25] %normalized CBL height bins
qtiles=[0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25] %normalized CBL height bins
qtiles2=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]

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

% M = [qtiles2;mchord; WM]
% csvwrite('test.csv',M)

%% Normalized

for ii=1:length(qtiles)-1
    figure(3);
    idxnow=[idx(ii).locs; idx2(ii).locs+53727];
    cmap=rbmapper_coarse(1.25,-.3);
    colormap(cmap);
    subplot(3,length(qtiles)-1,ii);
    WBAR1=squeeze(nanmean(Wnorm(idx(ii).locs,:,:),1));
    WBAR2=squeeze(nanmean(Wnorm2(idx2(ii).locs,:,:),1));
    WBAR=(WBAR1+WBAR2)./2;
    pcolor(ximat.*mdchord(ii)/2,zimat,WBAR);shading flat; 
    caxis([-.3 1.25]);
    xlabel('Median chord length [m]')
    ylabel('Z/Z_i')
    %contourf(ximat.*mdchord(ii)/2,zimat,WBAR,(-5:.1:5),'linestyle','none');shading flat; caxis([-.3 1.25]);

    hold on;
    %     contour(ximat.*mdchord(ii)/2,zimat,squeeze(nanmean(Wnorm(idx(ii).locs,:,:),1)),[0.5 .5],'k');shading flat; caxis([-.5 1.5]);
    
    ylim([.2 1.5]);
    %         xlim([-2 2]);
    xlim([-150 150]);
    grid on; box on; set(gca,'linewidth',2,'layer','top');
    text(-140,1.3,strcat('N=',num2str(length(idxnow))),'fontsize',17,'fontweight','bold')
        set(gca,'fontsize',15,'fontweight','bold','linewidth',2,'layer','top');
        title(strcat(num2str(qtiles(ii)),'<Z/Z_i<',num2str(qtiles(ii+1))))
    if ii==(length(qtiles)-1)
        cbh=colorbar
        cpos=get(cbh,'pos')
        set(cbh,'pos',[cpos(1)+.07 cpos(2) cpos(3) cpos(4)])
        ylabel(cbh,'m s^{-1}')
    end
end



%% Now examine size/intensity relationship for thermals in the upper CBL
f=figure(500);clf;set(f,'color','w');
nidx=[idx(3).locs; idx(4).locs; idx2(3).locs; idx2(4).locs];
UPWMAX=cat(1,upwmax,upwmax2);
UPWSTAR=cat(1,upwstar,upwstar2);
upwmaxn=UPWMAX(nidx)./UPWSTAR(nidx);
UPZI=cat(1,upZi,upZi2);
CCHORDN=(CCHORD_ALL(nidx))./UPZI(nidx);
histogram2(CCHORDN,upwmaxn,150,'normalization','probability','displaystyle','tile');%10,[.7 .7 .7],'filled');
hold on;
cbh=colorbar;ylabel(cbh,'prob','fontsize',12,'fontweight','bold');

xlim([0 4500]/3000);
ylim([0 4]);
cbins=[50:200:5000]/3000;
grid on; box on; set(gca,'linewidth',2);
clrs={'c','b','g','k','m','r'};

for cc=1:length(cbins)
    cidx=find(CCHORDN>= cbins(cc) & CCHORDN<cbins(cc)+.133/2);
    if numel(cidx)>50
        wbar=nanmean(upwmaxn(cidx));
        sc=scatter(cbins(cc)+.0667/2,wbar,80,'m','filled','markeredgecolor','k');
    end
end

ylim([0 3]);
xlabel('Chord Length/Z_i');
ylabel('Max Updraft/w^*');
set(gca,'fontsize',12,'fontweight','bold');
title('CBL updraft intensity vs size');
set(gca,'layer','top');
%%
f=figure(501);clf;set(f,'color','w');
nidx=[idx(5).locs; idx2(5).locs];
UPWMAX=cat(1,upwmax,upwmax2);
UPWSTAR=cat(1,upwstar,upwstar2);
upwmaxn=UPWMAX(nidx)./UPWSTAR(nidx);
UPZI=cat(1,upZi,upZi2);
CCHORDN=(CCHORD_ALL(nidx))./UPZI(nidx);
histogram2(CCHORDN,upwmaxn,200,'normalization','probability','displaystyle','tile');%10,[.7 .7 .7],'filled');
hold on;
cbh=colorbar;ylabel(cbh,'prob','fontsize',12,'fontweight','bold');

xlim([0 6000]/3000);
ylim([0 4]);
cbins=[0:400:5000]/3000;
grid on; box on; set(gca,'linewidth',2);
clrs={'c','b','g','k','m','r'};

for cc=1:length(cbins)
    cidx=find(CCHORDN>= cbins(cc) & CCHORDN<cbins(cc)+.133);
    if numel(cidx)>10
        wbar=nanmedian(upwmaxn(cidx));
        sc=scatter(cbins(cc)+.0667,wbar,70,'m','filled','markeredgecolor','k');
    end
end
hold on;
cvec=0:.1:1.2;
fitter=.6+(1.67*(cvec))-.83*(cvec.^2);
plot(cvec,fitter,'--k','linewidth',2);

ylim([0 4]);
xlabel('Chord Length/Z_i');
ylabel('Max Updraft/w^*');
set(gca,'fontsize',12,'fontweight','bold');
title('CBL updraft intensity vs size');

%%
% figure(502);clf;
% clrs=jet(10);
% for ii=1:length(qtiles)-1
%     hold on;
%     contour(ximat.*mdchord(ii)/2,zimat,squeeze(nanmedian(Wnorm(idx(ii).locs,:,:),1)),[.5 .5],'color',clrs(ii,:),'linewidth',3);
%     %         xlim([-2 2]);
%     
% end
% %
% ylim([.2 1.2])
% xlim([-150 150])
% grid on; box on; set(gca,'linewidth',2,'layer','top');
% ylabel('Z/Z_i','fontsize',18,'fontweight','bold');
% xlabel('Width [m]','fontsize',18,'fontweight','bold');

%%
f=figure(4);clf;set(f,'color','w');
clrs=jet(9);
qtiles2=[.85 .9 .95 1 1.05 1.1 1.15 1.2 1.25 1.3]
hold on;
for ii=1:(length(qtiles2)-1)
    if ii<length(qtiles2)
        idx3(ii).locs=find(upz_norm>=qtiles2(ii) & upz_norm<qtiles2(ii+1));
        idx4(ii).locs=find(upz_norm2>=qtiles2(ii) & upz_norm2<qtiles2(ii+1));
        idxnow=[idx3(ii).locs; idx4(ii).locs];
        mchord(ii)=round(nanmean(CCHORD_ALL(idxnow)));
        [~,zidx2]=min(abs(zideal-(qtiles2(ii))));
        WX1=squeeze(nanmax((nanmedian(Wnorm(idx3(ii).locs,:,:),1)),[],2));
        WX2=squeeze(nanmax((nanmedian(Wnorm2(idx4(ii).locs,:,:),1)),[],2));
        WX=0.5.*(WX1+WX2);
        a(ii)=plot(xideal*mchord(ii)/2,WX,'color',clrs(ii,:),'linewidth',2);
    end
    
end
xlim([-200 200]);
grid on; box on; set(gca,'linewidth',2,'fontsize',15,'fontweight','bold');
xlabel('Median Width [m]','fontsize',15,'fontweight','bold');
ylabel('W [m s^-1]', 'fontsize',15,'fontweight','bold');
legend(num2str(qtiles2(1)),num2str(qtiles2(2)),num2str(qtiles2(3)),num2str(qtiles2(4)),num2str(qtiles2(5)),num2str(qtiles2(6)),num2str(qtiles2(7)),num2str(qtiles2(8)));
return