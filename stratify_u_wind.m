close all
clear all

load('updraft_objects_20190227.mat','u_top','u_bot','v_top','v_bot','u_top2','u_bot2','v_top2','v_bot2','z34','z342','upwmax_top','upwmax_top2','upzbot','upzc','upzx','upxc','upZi','upctime','CWIDTH','upspd','upwmax','upwstar','upwmean','uparea','upzbot2','upzc2','upzx2','upxc2','upZi2','CWIDTH2','upspd2','upwmax2','upwstar2','upctime2','upwmean2','uparea2','xideal','zideal','ximat','zimat');
load('updraft_wnorm_20190227.mat','Wnorm','Wnorm2','xideal','zideal');

%% calculate wind shear
u_shear=abs(u_top-u_bot);
v_shear=abs(v_top-v_bot);

u_shear2=abs(u_top2-u_bot2);
v_shear2=abs(v_top2-v_bot2);

shear=u_shear+v_shear;
shear2=u_shear2+v_shear2;

%% split wspd into upper 1/3, middle 1/3, upper 1/3
normalized = shear %change between 1. upwstar, 2. upspd, 3. upwstar./upspd 4. shear
normalized2 = shear2 %change between 1. upwstar2, 2. upspd2, 3. upwstar2./upspd2 4. shear

lower = prctile(normalized,33);
middle = prctile(normalized,66);
higher = max(normalized)

windindx=nan(size(upspd))

windindx(normalized<=lower)=1
windindx(normalized>=lower & normalized<=middle)=2
windindx(normalized>=middle)=3

windindx2=nan(size(upspd2))

windindx2(normalized2<=lower)=1
windindx2(normalized2>=lower & normalized2<=middle)=2
windindx2(normalized2>=middle)=3

%% Compute upper 3/4 updraft location
upz_norm=(0.5.*(upzc+upzx))./upZi; %CBL height normalized 3/4 height.
upzxn=upzx./upZi;
upz_norm2=(0.5.*(upzc2+upzx2))./upZi2;
upzxn2=upzx2./upZi2;
% qtiles=[.25 .45 .65 .85 1.05 1.25] %normalized CBL height bins

% qtiles=[0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25] %normalized CBL height bins
% qtiles2=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]

qtiles=[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3]

CCHORD=CWIDTH(:).*upspd;%(nidx);
CCHORD2=CWIDTH2(:).*upspd2;
CCHORD_ALL=cat(1,CCHORD,CCHORD2);
UPWMAX=cat(2,upwmax,upwmax2);
WM=nan(3,length(qtiles));
mchord=nan(3,length(qtiles));
for ww=1:3
    clear idx idx2

    for ii=1:(length(qtiles)-1)
        if ii<length(qtiles)
            idx(ii).locs=find((upz_norm>=qtiles(ii) & upz_norm<qtiles(ii+1)) & windindx==ww); %indices for updrafts in a given height range (e.g., between .25 and .45 Zi)
            idx2(ii).locs=find((upz_norm2>=qtiles(ii) & upz_norm2<qtiles(ii+1)) & windindx2==ww);
        end
    end
    
    %%
    f=figure(3+ww);clf;set(f,'color','w');
    for ii=1:(length(qtiles)-1)
        sp=subplot(3,length(qtiles)-1,length(qtiles)-1+ii);
        cla(sp);
        idxnow=[idx(ii).locs; idx2(ii).locs+53727];
        mchord(ww,ii)=round(nanmean(CCHORD_ALL(idxnow)));
        mdchord(ww,ii)=round(nanmedian(CCHORD_ALL(idxnow)));
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
        text(800,.3,strcat('Mean=',num2str(mchord(ww,ii))),'fontsize',17,'fontweight','bold');
        text(800,.25,strcat('Median=',num2str(mdchord(ww,ii))),'fontsize',17,'fontweight','bold');
        set(gca,'fontsize',15,'fontweight','bold','linewidth',2,'layer','top');
        
        sp=subplot(3,length(qtiles)-1,2*(length(qtiles)-1)+ii);
        cla(sp);
        
        WM(ww,ii)=(nanmean(UPWMAX(idxnow)));
        WMD(ww,ii)=(nanmedian(UPWMAX(idxnow)));
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
        text(1,.35,strcat('Mean=',num2str(WM(ww,ii))),'fontsize',17,'fontweight','bold');
        text(1,.3,strcat('Median=',num2str(WMD(ww,ii))),'fontsize',17,'fontweight','bold');
        set(gca,'fontsize',15,'fontweight','bold','linewidth',2,'layer','top');
    end
end

%%
% 1. normalized by wstar
% % ----------------------------------------------------------------------
% %max updraft
% figure(15);clf;
% clrs={'r','b','g'}
% hold on;
% for ww=1:3
%     plot(WM(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
%     xlabel('Maximum updraft [m/s]','fontsize',18);
%     ylabel('Z/Z_i','fontsize',18);
% end
% legend('low','middle','high','fontsize',16)
% 
% axes('pos',[.65 .15 .25 .25]);
% histogram(upwstar,40,'Normalization','probability','facecolor','b');
% xline(lower,'linewidth',2);
% % text(lower, .067,'1/3','fontsize',10);
% xline(middle,'linewidth',2);
% % text(middle,.067,'2/3','fontsize',10);
% % xlabel('upwstar');
% ylabel('Prob.','fontsize',14);
% saveas(gcf,'powerpoint_plots/wstar_updraft.png')
% 
% %chord length
% figure(16);clf;
% clrs={'r','b','g'}
% hold on;
% for ww=1:3
%     plot(mchord(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
%     xlabel('Chord length [m]','fontsize',18);
%     ylabel('Z/Z_i','fontsize',18);
% end
% legend('low','middle','high','fontsize',16);
% 
% axes('pos',[.65 .15 .25 .25]);
% histogram(upwstar,40,'Normalization','probability','facecolor','b');
% xline(lower,'linewidth',2);
% % text(lower, .068,'1/3','fontsize',10);
% xline(middle,'linewidth',2);
% % text(middle,.068,'2/3','fontsize',10);
% % xlabel('upwstar');
% ylabel('Prob.','fontsize',14);
% saveas(gcf,'powerpoint_plots/wstar_chordlength.png')

% 2. normalized by upspeed
% % ----------------------------------------------------------------------
% %max updraft
% figure(17);clf;
% clrs={'r','b','g'}
% hold on;
% for ww=1:3
%     plot(WM(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
%     xlabel('Maximum updraft [m/s]','fontsize',18);
%     ylabel('Z/Z_i','fontsize',18);
% end
% legend('low','middle','high','fontsize',16);
% 
% axes('pos',[.65 .15 .25 .25]);
% histogram(upspd,40,'Normalization','probability','facecolor','b');
% xline(lower,'linewidth',2);
% % text(lower, .068,'1/3','fontsize',10');
% xline(middle,'linewidth',2);
% % text(middle,.068,'2/3','fontsize',10)
% % xlabel('upspd');
% ylabel('Prob.','fontsize',14);
% saveas(gcf,'powerpoint_plots/upspd_updraft.png')
% 
% %chord length
% figure(18);clf;
% clrs={'r','b','g'}
% hold on;
% for ww=1:3
%     plot(mchord(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
%     xlabel('Chord length [m]','fontsize',18);
%     ylabel('Z/Z_i','fontsize',18);
% end
% legend('low','middle','high','fontsize',16);
% 
% axes('pos',[.65 .15 .25 .25]);
% histogram(upspd,40,'Normalization','probability','facecolor','b');
% xline(lower,'linewidth',2);
% % text(lower, .068,'1/3','fontsize',10);
% xline(middle,'linewidth',2);
% % text(middle,.068,'2/3','fontsize',10)
% % xlabel('upspd');
% ylabel('Prob.','fontsize',14);
% saveas(gcf,'powerpoint_plots/upspd_chordlength.png')

% % % 3. normalized by upwstar/upspeed
% % % ----------------------------------------------------------------------
% %max updraft
% figure(17);clf;
% clrs={'r','b','g'}
% hold on;
% for ww=1:3
%     plot(WM(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
%     xlabel('Maximum updraft [m/s]','fontsize',18);
%     ylabel('Z/Z_i','fontsize',18);
% end
% legend('low','middle','high','fontsize',16);
% 
% axes('pos',[.65 .15 .25 .25]);
% histogram(upwstar./upspd,40,'Normalization','probability','facecolor','b');
% xline(lower,'linewidth',2);
% % text(lower, .11,'1/3','fontsize',10);
% xline(middle,'linewidth',2);
% % text(middle,.11,'2/3','fontsize',10)
% % xlabel('upwstar/upspd');
% ylabel('Prob.','fontsize',14);
% saveas(gcf,'powerpoint_plots/upwstar_upspd_updraft.png')
% 
% %chord length
% figure(18);clf;
% clrs={'r','b','g'}
% hold on;
% for ww=1:3
%     plot(mchord(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
%     xlabel('Chord length [m]','fontsize',18);
%     ylabel('Z/Z_i','fontsize',18);
% end
% legend('low','middle','high','fontsize',16);
% 
% axes('pos',[.65 .15 .25 .25]);
% histogram(upwstar./upspd,40,'Normalization','probability','facecolor','b');
% xline(lower,'linewidth',2);
% % text(lower, .11,'1/3','fontsize',10);
% xline(middle,'linewidth',2);
% % text(middle,.11,'2/3','fontsize',10)
% % xlabel('upwstar/upspd');
% ylabel('Prob.','fontsize',14);
% saveas(gcf,'powerpoint_plots/upwstar_upspd_chordlength.png')


% % % 4. normalized by shear
% % % ----------------------------------------------------------------------
% %max updraft
figure(17);clf;
clrs={'r','b','g'}
hold on;
for ww=1:3
    plot(WM(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
    xlabel('Maximum updraft [m/s]','fontsize',18);
    ylabel('Z/Z_i','fontsize',18);
end
legend('low','middle','high','fontsize',16);

axes('pos',[.65 .15 .25 .25]);
histogram(shear,40,'Normalization','probability','facecolor','b');
xline(lower,'linewidth',2);
% text(lower, .11,'1/3','fontsize',10);
xline(middle,'linewidth',2);
% text(middle,.11,'2/3','fontsize',10)
% xlabel('upwstar/upspd');
ylabel('Prob.','fontsize',14);
saveas(gcf,'powerpoint_plots/shear_updraft.png')

%chord length
figure(18);clf;
clrs={'r','b','g'}
hold on;
for ww=1:3
    plot(mchord(ww,:),qtiles,'color',clrs{ww},'linewidth',2)
    xlabel('Chord length [m]','fontsize',18);
    ylabel('Z/Z_i','fontsize',18);
end
legend('low','middle','high','fontsize',16);

axes('pos',[.65 .15 .25 .25]);
histogram(shear,40,'Normalization','probability','facecolor','b');
xline(lower,'linewidth',2);
% text(lower, .11,'1/3','fontsize',10);
xline(middle,'linewidth',2);
% text(middle,.11,'2/3','fontsize',10)
% xlabel('upwstar/upspd');
ylabel('Prob.','fontsize',14);
saveas(gcf,'powerpoint_plots/shear_chordlength.png')