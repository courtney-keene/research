%script to plot composite CBL updrafts at SGP

close all
clear all;

% save updraft object output data
load('updraft_objects_20190227.mat','z34','z342','upwmax_top','upzbot','upzbot2','upwmax_top2','upzc','upzx','upxc','upctime','upZi','CWIDTH','upspd','upwmax','upwstar','upwmean','uparea','upzc2','upzx2','upxc2','upZi2','CWIDTH2','upspd2','upwmax2','upctime2','upwstar2','upwmean2','uparea2','xideal','zideal');
load('updraft_wnorm_20190227.mat','Wnorm','Wnorm2','xideal','zideal');
[ximat,zimat]=meshgrid(xideal,zideal);

%% Compute upper 3/4 updraft location
upz_norm=(0.5.*(upzc+upzx))./upZi;
upzxn=upzx./upZi;
upz_norm2=(0.5.*(upzc2+upzx2))./upZi2;
upzxn2=upzx2./upZi2;

updnorm=(upzx-upzbot)./upZi;
updnorm2=(upzx2-upzbot2)./upZi2;

upzbot_norm=upzbot./upZi;
upzbot_norm2=upzbot2./upZi2;

upzx_norm=upzx./upZi;
upzx_norm2=upzx2./upZi2;

%qtiles=[.25 .45 .65 .85 1.05 1.25];

% qtiles=[.1 .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2 1.3 1.4]
% qtiles_=[.15 .25 .35 .45 .55 .65 .75 .85 .95 1.5 1.15 1.25 1.35]
% index=[0 1 2 3 4 5 6 7 8 9 10 11 12]

qtiles=[.1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9 .95 1.0 1.05 1.1 1.15 1.2 1.25 1.3];
qtiles_=[.125 .175 .225 .275 .325 .375 .425 .475 .525 .575 .625 .675 .725 .775 .825 .875 .925 .975 1.025 1.075 1.125 1.175 1.225 1.275];
index=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];

for ii=1:(length(qtiles)-1)
    if ii<length(qtiles)
        idx(ii).locs=find(upz_norm>=qtiles(ii) & upz_norm<qtiles(ii+1));
        idx2(ii).locs=find(upz_norm2>=qtiles(ii) & upz_norm2<qtiles(ii+1));
    end
end

%%
f=figure(3);clf;set(f,'color','w');
CCHORD=CWIDTH(:).*upspd;%(nidx);
CCHORD2=CWIDTH2(:).*upspd(2);
CCHORD_ALL=cat(1,CCHORD,CCHORD2);
UPWMAX=cat(2,upwmax_top,upwmax_top2);
for ii=1:(length(qtiles)-1)
    sp=subplot(3,length(qtiles)-1,length(qtiles)-1+ii);
    cla(sp);
    idxnow=[idx(ii).locs; idx2(ii).locs+53727];
    mchord(ii)=round(nanmean(CCHORD_ALL(idxnow)));
    mdchord(ii)=round(nanmedian(CCHORD_ALL(idxnow)));
    stdchord(ii)=round(nanstd(CCHORD_ALL(idxnow)));
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
    text(300,.37,strcat('Mean=',num2str(mchord(ii)),{' m'}),'fontsize',12,'fontweight','bold');
    text(300,.32,strcat('Median=',num2str(mdchord(ii)),{' m'}),'fontsize',12,'fontweight','bold');
    set(gca,'fontsize',12,'fontweight','bold','linewidth',2,'layer','top');
    
    sp=subplot(3,length(qtiles)-1,2*(length(qtiles)-1)+ii);
    cla(sp);
    
    WM(ii)=(nanmean(UPWMAX(idxnow)));
    WMD(ii)=(nanmedian(UPWMAX(idxnow)));
    Wstd(ii)=(nanstd(UPWMAX(idxnow)));
    N(ii)=length(idxnow);
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
    text(1.2,.37,strcat('Mean=',num2str(WM(ii)),{' m s^{-1}'}),'fontsize',12,'fontweight','bold');
    text(1.2,.32,strcat('Median=',num2str(WMD(ii)),{' m s^{-1}'}),'fontsize',12,'fontweight','bold');
    set(gca,'fontsize',12,'fontweight','bold','linewidth',2,'layer','top');

end

%%

% 
% M=[index;N;qtiles_;mchord;mdchord;stdchord;WM;WMD;Wstd]
% M_=transpose(M)
% row=1
% col=0
% csvwrite('coh_up_plum_plotter.csv',M_,row,col)

%% Normalized
tstr(1).txt='.25<=Z/Zi<.45';
tstr(2).txt='.45<=Z/Zi<.65';
tstr(3).txt='.65<=Z/Zi<.85';
tstr(4).txt='.85<=Z/Zi<1.05';
tstr(5).txt='1.05<=Z/Zi<1.25';


for ii=1:(length(qtiles)-1)
    figure(3);
    idxnow=[idx(ii).locs; idx2(ii).locs+53727];
    cmap=rbmapper_coarse(1.25,-.35);
    colormap(cmap);
    subplot(3,length(qtiles)-1,ii);
    WBAR1=squeeze(nanmean(Wnorm(idx(ii).locs,:,:),1));
    WBAR2=squeeze(nanmean(Wnorm2(idx2(ii).locs,:,:),1));
    WBAR=(WBAR1+WBAR2)./2;
    %     pcolor(ximat.*mdchord(ii)/2,zimat,WBAR);shading flat; caxis([-.3 1.25]);
    contourf(ximat.*mdchord(ii)/2,zimat,WBAR,[-5:.1:5],'linestyle','none');shading flat; caxis([-.35 1.25]); 
    if ii==(length(qtiles)-1)
        cbh=colorbar; 
        cpos=get(cbh,'pos');
        set(cbh,'pos',[cpos(1)+.05 cpos(2) cpos(3) cpos(4)]);
        ylabel(cbh,'m s^-1','fontsize',15,'fontweight','bold');
    end
    if ii==1
        ylabel('Z/Z_i','fontsize',15,'fontweight','bold');
    end
    xlabel('Median Chord Length','fontsize',15,'fontweight','bold');
    %title(tstr(ii).txt,'fontsize',15,'fontweight','bold');
    hold on;
    %     contour(ximat.*mdchord(ii)/2,zimat,squeeze(nanmean(Wnorm(idx(ii).locs,:,:),1)),[0.5 .5],'k');shading flat; caxis([-.5 1.5]);
    
    ylim([.2 1.5])
    %         xlim([-2 2]);
    xlim([-150 150])
    grid on; box on; set(gca,'linewidth',2,'layer','top');
    text(-140,1.3,strcat('N=',num2str(length(idxnow))),'fontsize',15,'fontweight','bold')
    set(gca,'fontsize',12,'fontweight','bold','linewidth',2,'layer','top');
    
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
        sc=scatter(cbins(cc)+.0667/2,wbar,70,'m','filled','markeredgecolor','k');
    end
end

ylim([0 3]);
xlabel('Chord Length/Z_i');
ylabel('Max Updraft/w^*');


%% 

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

% %% LOAD CIN DATA FOR 1700-1800 hrs
% load('CIN500_cbl.mat','CIN','CINtime');
% CIN(CIN<-2500)=-2500; 
% CIN2=sqrt(-CIN);
% 
% [~,~,~,HH1,~]=datevec(upctime);
% [~,~,~,HH2,~]=datevec(upctime2);
% 
% tidx1=find(HH1==17);% & HH1<20);
% tidx2=find(HH2==17);% & HH2<20);
% upz1=upz_norm(tidx1);
% upz2=upz_norm2(tidx2);
% mcall=cat(1,CCHORD(tidx1),CCHORD2(tidx2));
% ziall=cat(1,upZi(tidx1),upZi2(tidx2));
% W1=squeeze(Wnorm(tidx1,:,:));
% W2=squeeze(Wnorm2(tidx2,:,:));
% 
% UPW=cat(1,upwstar(tidx1),upwstar2(tidx2));
% upwmat=repmat(UPW,1,61,161);
% Wall=cat(1,W1,W2);
% wnorm=Wall./upwmat; 
% upzall=cat(1,upz1,upz2);
% timeall=cat(1,upctime(tidx1),upctime2(tidx2));
% draft_cin=nan(size(timeall));
% for tt=1:length(timeall);
%     %find the nearest sounding
%     [hdst,ttidx]=min(abs((timeall(tt)-(CINtime))));
%     if ttidx && hdst<1
%         draft_cin(tt)=CIN2(ttidx);
%     end
%     
% end
% %%
% figure(500);clf;hold on;
% clrs=jet(8);
% qtiles2=[.8 .9 1. 1.1 1.2];
% qtiles3=[-200 -100 0 100 200 300];
% 
% 
% upz341=(0.5.*(upzc+upzx));
% upz342=(0.5.*(upzc2+upzx2));
% upz34=cat(1,upz341,upz342);
% %% Create CIN bins (pun intended)
% clear ccc;
% cinbin=prctile(draft_cin,[50]); %set the percentile bins for CIN at the CBL top 
% ccc(1).locs=find(draft_cin<cinbin(1)); %find indices for times experienin
% 
% % ccc(2).locs=find(draft_cin>=cinbin(1) & draft_cin<cinbin(2));
% ccc(2).locs=find(draft_cin>=cinbin(1));
% 
% %%
% figure(22);clf;
% upZi_all=cat(1,upZi,upZi2);
% clear upznow chordnow Wnow cinnow;
% for cb=1:2%length(cinbin);
%     upznow=upzall(ccc(cb).locs);
%     chordnow=mcall(ccc(cb).locs);
%     zinow=(upZi_all(ccc(cb).locs));
%     Wnow=wnorm(ccc(cb).locs,:,:);
%     cinnow=draft_cin(ccc(cb).locs);
%     for ii=1:(length(qtiles2))
%         if ii<length(qtiles2)
%             idxnow=find(upznow>=(qtiles3(ii)+zinow) & upznow<(zinow+qtiles3(ii+1)));
%             mcc(ii)=round(nanmean(chordnow(idxnow)));
%             WX1=squeeze(((nanmean(Wnow(idxnow,:,:),1))));
%             %%
%             %             WXCIN=WX1./sqrt(-1.*nanmean(cinnow(idxnow)));
%             subplot(1,2,cb);hold on;
% %             WX1=inpaint_nans(WX1);
%             [bb,cc]=butter(5,1/5);
%             contour(xideal*mcc(ii)/2,zideal.*zinow,WX1,[0.1 0.1],'color',clrs(ii,:));hold on;
%             hold on; 
%             plot([-500 500],[zinow zinow],'--k','linewidth',3);
% %             plot(xideal*mcc(ii)/2,WX1,'color',clrs(ii,:),'linewidth',2);
%              xlim([-200 200]);
%             grid on; box on;
% %             ylim([-.3 .5]);
%         end
%         
%     end
% end
% 
% %%
% allarea=[uparea; uparea2];
% allarea(allarea==0)=nan;
% figure;
% histogram2(allarea,UPWMAX,50,'normalization','probability','displaystyle','tile');
% % xlim([