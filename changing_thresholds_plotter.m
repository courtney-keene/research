%script to plot composite CBL updrafts at SGP

close all
clear all;

% save updraft object output data
load('updraft_objects_thresholds.mat','Z34','upzbot','upzc','upzx','upxc','upZi','upctime','CWIDTH','upspd','upwmax','upwstar','upwmean','uparea','z34_','upwmax_top','hgt_wmax','xideal','zideal','ximat','zimat');
load('updraft_objects_thresholds23.mat','Z342','Z343','upzbot2','upzc2','upzx2','upxc2','upZi2','upctime2','CWIDTH2','upspd2','upwmax2','upwstar2','upwmean2','uparea2','z342','upwmax_top2','hgt_wmax2','upzbot3','upzc3','upzx3','upxc3','upZi3','upctime3','CWIDTH3','upspd3','upwmax3','upwstar3','upwmean3','uparea3','z343','upwmax_top3','hgt_wmax3');
load('updraft_objects_thresholds4.mat','Z344','upzbot4','upzc4','upzx4','upxc4','upZi4','upctime4','CWIDTH4','upspd4','upwmax4','upwstar4','upwmean4','uparea4','z344','upwmax_top4','hgt_wmax4');
load('updraft_wnorm_thresholds.mat','Wnorm','Wnorm2','xideal','zideal');

[ximat,zimat]=meshgrid(xideal,zideal);


%% Compute upper 3/4 updraft location
upz_norm=(0.5.*(upzc+upzx))./upZi;
upzxn=upzx./upZi;


upz_norm2=(0.5.*(upzc2+upzx2))./upZi2;
upzxn2=upzx2./upZi2;

upz_norm3=(0.5.*(upzc3+upzx3))./upZi3;
upzxn3=upzx3./upZi3;


upz_norm4=(0.5.*(upzc4+upzx4))./upZi4;
upzxn4=upzx4./upZi4;

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
        idx3(ii).locs=find(upz_norm3>=qtiles(ii) & upz_norm3<qtiles(ii+1));
        idx4(ii).locs=find(upz_norm4>=qtiles(ii) & upz_norm4<qtiles(ii+1));
    end
end

%%
f=figure(3);clf;set(f,'color','w');
CCHORD=CWIDTH(:).*upspd;%(nidx);
CCHORD2=CWIDTH2(:).*upspd(2);
CCHORD3=CWIDTH3(:).*upspd(3);
CCHORD4=CWIDTH4(:).*upspd(4);

CCHORD_ALL=cat(1,CCHORD,CCHORD2,CCHORD3,CCHORD4);
UPWMAX=cat(1,upwmax_top,upwmax_top2,upwmax_top3,upwmax_top4);
upzxs=cat(1,upzx,upzx2,upzx3,upzx4)
upbothgt = cat(1,upzbot,upzbot2,upzbot3,upzbot4);
Zi_s = cat(1,upZi,upZi2,upZi3,upZi4);
z34s= cat(1,Z34,Z342,Z343,Z344);

% for ii=1:(length(upbothgt));
%     if upbothgt(ii)<Zi_s(ii)*0.8;
%         cnt = cnt+1;
%     end

% end
        
        
for ii=1:(length(qtiles)-1)
    idxnow=[idx(ii).locs; idx2(ii).locs+53727];
    mchord(ii)=round(nanmean(CCHORD_ALL(idxnow)));
    mdchord(ii)=round(nanmedian(CCHORD_ALL(idxnow)));
    stdchord(ii)=round(nanstd(CCHORD_ALL(idxnow)));
    WM(ii)=(nanmean(UPWMAX(idxnow)));
    WMD(ii)=(nanmedian(UPWMAX(idxnow)));
    Wstd(ii)=(nanstd(UPWMAX(idxnow)));
    N(ii)=length(idxnow);
end

%%


M=[index;N;qtiles_;mchord;mdchord;stdchord;WM;WMD;Wstd]
M_=transpose(M)
row=1
col=0
csvwrite('coh_up_plum_counter_nv7.csv',M_,row,col)

% M=[CCHORD_ALL,UPWMAX,upzxs,upbothgt,Zi_s,z34s];
% % M_=transpose(M)
% row=1
% col=0
% csvwrite('updrafts_widths_hgts.csv',M,row,col)