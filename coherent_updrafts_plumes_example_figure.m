%% Script to compute updraft properties at different heights in the CBL
close all
clear all

%addpath('~/Dropbox/MATLAB/');

daysec=24*60*60; %number of seconds in a day
pflag=1; % flag to toggle plotting on/off on=1
cmap=[gray(12);flipud(gray(12))]
cmap(1:12,1)=1
cmap(13:end,3)=1
target_years={'2015' '2016'}; %select the year of data you want to process

%% Load locally forced convection dates to consider
% test=dlmread('~/Dropbox/Shallow_to_Deep/Specturm_of_depths.csv','/'); %this list has various categories of convective development
% dates=datenum(test(:,3)+2000,test(:,1),test(:,2)); %convert to matlab time
% cudays_all=unique(dates); %make sure the list is unique

%% List of Lidar Sites at SGP
sites={'E32', 'E37', 'E39', 'E41', 'C1'};

% %% idealized time/height vectors for use later in compositing updrafts
% idealtime=-1.75:.005:1.75; % a cloud centered time vector at arbitrary resolution scaled to time of cloud overpass
% hideal=0:.005:.5; % a cloud base relative height vector at arbitrary resolution
% [tmati,hmati]=meshgrid(idealtime,hideal); %grid the interpolation values


%% Initialize some empty matrices for compositing of draft properties
%size is arbitrarily large to avoid space constraints later
upctime=nan(53727,1); %updraft center times
updur=nan(53727,1); %updraft duration
upchord=nan(53727,1); %updraft chord length
upwmax=nan(53727,1); %updraft maximum velocity
upwmean=nan(53727,1); %updraft mean velocity
upspd=nan(53727,1); %mean CBL wind speed for each updraft
upwstar=nan(53727,1); %CBL wstar value
upz=nan(53727,1); %updraft centroid height;
uparea=nan(53727,1);%updraft area (m^2 time-> space converted using wind speed)
draft_count=0; % keep track of the numer of updrafts in our sample
upZi=nan(53727,1); %CBL height for each updraft.
% uptilt=nan(53727,1); %orientation of updraft object
upxc=nan(53727,1); %updraft "x" center
upzc=nan(53727,1); %updraft "z" center
upzx=nan(53727,1); %updraft "z" corresponding to maximum updraft
% wsnapshot=nan(53727,100,501); %this is for parsing out the updraft scene for each updraft
Wnorm=nan(53727,61,161); %height/time normalized updraft scene
CWIDTH=nan(53727,1); %updraft width matrix
wxz34=nan(53727,1);%max updraft at the updraft 3/4 height
z34=nan(53727,1); %3/4 height
upwmax_top=nan(53727,1); %max updraft in the top quarter of updraft
hgt_wmax=nan(53727,1);
top_spd=nan(53727,1);
bot_spd=nan(53727,1);
theta_top=nan(53727,1);
theta_bot=nan(53727,1);
u_top=nan(53727,1);
u_bot=nan(53727,1);
v_top=nan(53727,1);
v_bot=nan(53727,1);

upctime2=nan(53727,1); %updraft center times
updur2=nan(53727,1); %updraft duration
upchord2=nan(53727,1); %updraft chord length
upwmax2=nan(53727,1); %updraft maximum velocity
upwmean2=nan(53727,1); %updraft mean velocity
upspd2=nan(53727,1); %mean CBL wind speed for each updraft
upwstar2=nan(53727,1); %CBL wstar value
upz2=nan(53727,1); %updraft centroid height;
uparea2=nan(53727,1);%updraft area (m^2 time-> space converted using wind speed)
upZi2=nan(53727,1); %CBL height for each updraft.
% uptilt2=nan(53727,1); %orientation of updraft object
upxc2=nan(53727,1); %updraft "x" center
upzc2=nan(53727,1); %updraft "z" center
upzx2=nan(53727,1); %updraft "z" corresponding to maximum updraft
% wsnapshot2=nan(53727,100,501); %this is for parsing out the updraft scene for each updraft
Wnorm2=nan(53727,61,161); %height/time normalized updraft scene
CWIDTH2=nan(53727,1); %updraft width matrix
wxz342=nan(53727,1);%max updraft at the updraft 3/4 height
z342=nan(53727,1); %3/4 height
upwmax_top2=nan(53727,1); %max updraft in the top quarter of updraft
hgt_wmax2=nan(53727,1);
top_spd2=nan(53727,1);
bot_spd2=nan(53727,1);
theta_top2=nan(53727,1);
theta_bot2=nan(53727,1);
u_top2=nan(53727,1);
u_bot2=nan(53727,1);
v_top2=nan(53727,1);
v_bot2=nan(53727,1);

files=dir('/Users/courtneykeene/Desktop/Research/lidar_data/SGP_LIDAR_DATA/RAW_DATA/C1/*.201*.cdf');


%% load at CBL heights
load('/Users/courtneykeene/Desktop/Research/C1_2011_2018_CBLH.mat','Zi_all','Time_all','Wstar_all','CBH25_all','CF_all','bhgts');


%% ***** PLEASE NOTE THESE ARE 1 Sec DATA AND PROCESSING MORE THAN A FEW DAYS CAN GET CUMBERSOME ON A PERSONAL MACHINE **************
% index=dlmread('shallow_cumulus_days_2016_2017_lidar_check.csv',',',1,0); %SHCU DAYS DURING PERIOD OF LIDAR OVERLAP
% dates=datenum(index(:,1),index(:,2),index(:,3)); %CONVERT TO MATLAB DATE
% gb=index(:,end); %this is an integer 0-5 indicating the number of lidar locations showing favorable CBL characteristics
% dates(gb<5)=[]; %REMOVE DAYS WHERE NOT ALL LIDARS ARE FAVORABLE.
% cudays=(datenum(dates)); %CREATE A LIST OF DAYS OF LIDAR DATA TO PROCESS
% datestr(cudays(25))


%startday=datenum(2013,7,24);

%  startday=datenum(2016,6,11);%datenum(2013,7,21);
%startday=datenum(2017,6,1);%datenum(2013,7,21);
startday=datenum(2016,6,11);%datenum(2013,7,21);
endday= datenum(2016,6,11);%startday+(26/24);%datenum(2015,4,3);
daylist=startday:1:endday;
thr=02;
%startday=datenum(2015,9,06);
%endday= datenum(2015,9,07);
cmap=[gray(12);flipud(gray(12))];
cmap(1:12,3)=1;
cmap(13:end,1)=1;
%load('vert_vel_map');



%% PROCESS LIDAR DATA
%% load the lidar boundary layer wind data (produced in "lidar_wind.m")
fin=strcat('/Users/courtneykeene/Desktop/Research/lidar_data/SGP_LIDAR_DATA/WINDS/','lidar_winds_',sites{5},'.mat');
load(fin,'wspdmaster','wdmaster','timemaster','hgts')
whgts=hgts; clear hgts; %rename some of the variables
Wall=[];
timeall=[];
f=figure(1);clf;set(f,'color','w','pos',[100 100 1500 700]);
for ff=1:length(files)
    
    %% read in one hour of raw lidar data 
    fname=strcat('/Users/courtneykeene/Desktop/Research/lidar_data/SGP_LIDAR_DATA/RAW_DATA/C1/',files(ff).name);
    yyyy=str2double(fname(end-18:end-15));
    month=str2double(fname(end-14:end-13));
    day=str2double(fname(end-12:end-11));
    hour=str2double(fname(end-9:end-8));
    fileday=datenum(yyyy,month,day);
    filetime=datenum(yyyy,month,day,hour,0,0);
    
    if filetime>=(startday+(19/24)) && filetime<endday+(21/24) %&& hour>=15 && hour<23 % thr-5  && hour<22%thr+4
        ncdisp(fname);
        
        B=double(ncread(fname,'attenuated_backscatter'));
        W=double(ncread(fname,'radial_velocity'));
        if size(W,2)>1
            %filter W using a range and/or standard deviation filter
            %WR=rangefilt(W,ones(3,3));
            %W(del2(W)>.05)=nan;
            W(1:2,:)=nan;
            SNR=double(ncread(fname,'intensity'));
            W(SNR<1.01)=nan;
            
            if yyyy<=2012
                time=ncread(fname,'time_offset');
            else
                time=ncread(fname,'time');
            end
            basetime=double(ncread(fname,'base_time'));
            %convert to matlab serial time
            basetime=(basetime./(60*60*24))+datenum(1970,1,1);
            time=(time./(60*60*24))+basetime;
            time_ind=find(timemaster>=time(1) & timemaster<time(end));
            if numel(time_ind)<1
                speed_bar= 4;
            else
                speed_bar=nanmean(nanmean(wspdmaster(1:10,time_ind),1),2);
            end
            timeall=cat(1,timeall,time);
            W(:,gradient(time)>(1.5.*mode(gradient(time))))=nan;
            %B(:,gradient(time)>(1.2.*mode(gradient(time))))=nan;
            bb=ones(3,3); bb=bb./numel(bb);
%             W=filter2(bb,W);
            W=medfilt2(W);
            Wall=cat(2,Wall,W);
            
            
            dz=30; %gate length
            base_height=318.000000; % this can be parsed from the global attributes
            hgts=base_height+ dz./2 + dz.*(1:size(B,1)); %constructs a height vector
            
            sp1=subplot(2,1,1);
            hold on;
            %colormap(vert_vel_map);
%             cmap=jet(12);
            %pcolor(time,hgts,medfilt2(W,[1 1]));shading flat; caxis([-5 5]);%caxis([-5 5]);
            
            pp=pcolor(time,hgts,W);shading flat; caxis([-2 2]);%caxis([-5 5]);
            colormap(sp1,cmap)
            %alpha(pp,.75);
            hold on;
            [tmat,hmat]=meshgrid(time,hgts);
            if yyyy<2015
                tmat(B<.25.*10^-4)=[];hmat(B<.25*10^-4)=[];
            else
                tmat(B<(.2*10^-3))=[];hmat(B<(.2*10^-3))=[];
            end
            
            scatter(tmat(:),hmat(:),5,'ok','linewidth',2);%,'filled','markeredgecolor','k');
            %             xlim([startday+1.25/24 startday+3/24]);
            hold on;
            ylim([base_height 3500]);
            %             plot(Time_all(zi_idx),Zi_all(zi_idx),'--k','linewidth',2);
            
            
            %             subplot(2,1,1);
            %             hold on;
            %             scatter(tmat(:),hmat(:),.5,'*k');%,'filled','markeredgecolor','k');
            
        end
    end
end
%%
xlim([startday+19/24 startday+21/24]);
datetick('x','keeplimits');
ylim([300 2500]);
xlabel('Time [UTC]')
ylabel('Height [m MSL]')
grid on; box on; set(gca,'layer','top','linewidth',2,'fontsize',15,'fontweight','bold');
cbh=colorbar; ylabel(cbh,'m s^{-1}','fontsize',18,'fontweight','bold');

%% USE IMAGE PROCESSING TO LOCATE TEMPORALLY AND SPATIALLY COHERENT UPDRAFT OBJECTS


%### create a binary mask for vertical velocity greater
%than some threshold vertical velocity:
wthreshold=0.25; %define the updraft threshold value
W=Wall;
time=timeall;
hgt=hgts;

Wmask=W; Wmask(W<wthreshold)=false; Wmask(W>=wthreshold)=true; %create true/false mask for points on either side of the threshold
%set edges to false
Wmask(1,:) = false;
Wmask(end,:) = false;
Wmask(:,1) = false;
Wmask(:,end) = false;
Wmask(isnan(Wmask))=false;
Wmask=logical(Wmask); %make sure the whole array is logical true/false, which is the same as white/black

if pflag==1
    %     figure(100);clf;
    %     sp1=subplot(2,1,1);
    %     rbm=rbmapper(5,-5);
    %     colormap(sp1,rbm);
    %     pcolor(time,hgts,medfilt2(W,[1 1]));shading flat; caxis([-3 3]);
    %     hold on;
    %     ylim([0 3500]);
    %     %                         sp2=subplot(3,1,2);
    %
    %                         pcolor(time,hgts,Wmask);shading flat;
end
horiz_threshold = round(100/speed_bar); %number of points to reach 100 meters based on the windspeed
size_threshold=12*horiz_threshold; 
upperbound=.1*(round(2000/speed_bar)^2); %figure out appropriate upperbound
Wmaskfilt=bwareafilt(Wmask,[size_threshold upperbound]); %filter on an area of "size_threshold" pixels (what remains we will call "coherent" updrafts)
updrafts=bwlabel(Wmaskfilt,4); %label each independent connected updraft region
% upfrac(cuidx,hour,:)=nansum(Wmaskfilt,2)./size(Wmaskfilt,2); %this is the coherent updraft fraction for the hour as a function of height
% cblhr(cuidx,hour)=cblhnow; %this stores the CBL height for each hour
Wcopy=W; Wcopy(Wmaskfilt==0)=nan; %now use the updraft mask to blank downdraft values for

% Leverage MATLAB regionprop algorithm to compute
% properties of each updraft object (kind of a black
% box)
regions=regionprops(updrafts,'Area','Extrema','BoundingBox','Centroid','PixelList','PixelIdxList','Orientation','convexhull','Eccentricity','MajorAxisLength','MinorAxisLength')
%weight_cent = regionprops(updrafts,I,{'WeightedCentroid'})
%## Optional plotting to show the unique updraft
%objects
if pflag==1
    sp3=subplot(2,1,2);
    colormap(sp3,lines(100));
    upplot=updrafts;upplot(upplot==0)=nan;
    pcolor(time,hgts,upplot);shading flat;%
    ylim([0 3500]);
end

%Parse out various properties
centers=[regions.Centroid];
areas=[regions.Area];

% loop over each updraft object (contained in regions)
if numel(regions)>0
    for rr=1:length(regions)
        if regions(rr).Area>size_threshold %must be at least 100 pixels in the area
            draft_count=draft_count+1; %count the number of coherent updrafts
            idxnow=regions(rr).PixelList; %list of pixels
            %scatter(idxnow(:,1),idxnow(:,2),'*k');
            %regions(rr).Wmax=nanmax(W(regions(rr).PixelIdxList));
            
            [upwmax(draft_count),upidx]=nanmax(W(regions(rr).PixelIdxList)); %find the maximum vertical velocity value in the updraft
            [upi,upj]=ind2sub(size(W),regions(rr).PixelIdxList(upidx))
            hgt_wmax(draft_count)=hgt(upi)
            upwmean(draft_count)=nanmean(W(regions(rr).PixelIdxList)); %compute the mean vertical velocity in the updraft
            %                 upwstar(draft_count)=wstarnow; %assign each updraft a convective velocity scale based on the hourly values computed previously (this is just for reference and mormalization and not a property of the updraft)
            %                 upspd(draft_count)=spdbar; %mean CBL wind speed
            %                 pixel_area=(spdbar*1).*mode(gradient(hgts));   %this is the "size" of a given pixek = speed*1 second (meters) x gate length (30 m)
            %                 uparea(draft_count)=regions(rr).Area.*pixel_area; %the total updraft area is therefore the #of pixels (regions(rr).Area) multiplied by the "pixel size"
            %                 upZi(draft_count)=cblhnow; % cbl height corresponding to each updraft (not locally computed).
            upz(draft_count)=hgts(round(regions(rr).Centroid(2))); %the height of the center of the updraft
            hold on;
            [ztop,zidx]=max(regions(rr).ConvexHull(:,2)); % ztop is the "pixel" top, but not in whole pixels... strange, but workable
            zcenter=round(regions(rr).Centroid(2)); %must be rounded because it can be non integer number (whereas actual data is integer spacing)
            xcenter=round(regions(rr).Centroid(1));
            width=max(regions(rr).Extrema(:,1))-min(regions(rr).Extrema(:,1)); % maximum width (but not necessarily along a single height)
            zpoints=regions(rr).PixelList(:,2);
            xpoints=regions(rr).PixelList(:,1);
            cidx=find(zpoints==(round((zcenter+ztop)/2))); %this is the location half way between the upper edge of the updraft and the updraft center (upper quarter locaton)
            LC=min(xpoints(cidx)); %left edge at the upper quarter location
            RC=max(xpoints(cidx)); %right edge at the upper quarter location
            
            CWIDTH(draft_count)=RC-LC; %width of the updraft at the upper quarter location
            z34(draft_count)=round((zcenter+ztop)/2); %three quarter height of the updraft
            upxc(draft_count)=LC+(RC-LC)/2;%xcenter; %xcenter of the updraft at 3/4 height
            upzc(draft_count)=hgts(zcenter); %zcenter of the updraft
            upzx(draft_count)=hgts(round(ztop)); %top of the updraft (x here indicates maximum height)
            upctime(draft_count)=time(round(upxc(draft_count)));
            
            upwmax_top(draft_count)=nanmax(W(regions(rr).PixelIdxList(zpoints>=z34(draft_count))));
            %% normalized snap shot (normalized in time and in CBL height (Zi))
            if isfinite(CWIDTH(draft_count)) && CWIDTH(draft_count)>0
                zvec=(hgts(1:100));
                if isfinite(upZi(draft_count))
                    znorm=zvec./upZi(draft_count);
                else
                    znorm=zvec./2000;
                end
                xnorm=((1:length(time))-xcenter)./CWIDTH(draft_count); %this is a time (x-data) normalization centered on the updraft center period
                xideal=-2:.025:2;
                zideal=0:.025:1.5;
                [xmat,zmat]=meshgrid(xnorm,znorm);
                [ximat,zimat]=meshgrid(xideal,zideal);
                Wnorm(draft_count,:,:)=interp2(xmat,zmat,W(1:100,:),ximat,zimat);
            end
            
            %%
%             z34=round((ztop+zcenter)/2);
            phi = linspace(0,2*pi,50);
            cosphi = cos(phi);
            sinphi = sin(phi);
            if pflag==1
                hold on;
                plot(time(round(regions(rr).ConvexHull(:,1))),hgt(round(regions(rr).ConvexHull(:,2))),'k','linewidth',2);
%                 plot(time(round(regions(rr).Centroid(1))), hgt(round(regions(rr).Centroid(2))),'*w','markersize',5);
%                 plot(time(round([regions(rr).Centroid(1)-width/2 regions(rr).Centroid(1)+width/2 ])), hgt(round([regions(rr).Centroid(2) regions(rr).Centroid(2)])),'--k');
                plot(time(round([LC RC])),hgt(round([zcenter zcenter])),'--r');
                plot(time(round(upxc(draft_count))),hgt(round(z34(draft_count))),'*m','markersize',5,'linewidth',2);
                plot(time(round([LC RC])),hgt(round([z34(draft_count) z34(draft_count)])),'-k','linewidth',2);
                plot(regions(rr).ConvexHull(zidx,1), regions(rr).ConvexHull(zidx,2),'*r','markersize',15);
                %                     ylim([0 120]);
%                 plot(time(upj(draft_count)),hgt(upi(draft_count)),'*k')
                xbar = regions(rr).Centroid(1);
                ybar = regions(rr).Centroid(2);

                a = regions(rr).MajorAxisLength/2;
                b = regions(rr).MinorAxisLength/2;

                theta = pi*regions(rr).Orientation/180;
                R = [ cos(theta)   sin(theta)
                     -sin(theta)   cos(theta)];

                xy = [a*cosphi; b*sinphi];
                xy = R*xy;

                x = xy(1,:) + xbar;
                y = xy(2,:) + ybar;
                x(x<=0.5)=nan;
                y(y<=0.5)=nan;
                nanner=x+y;
                x(isnan(nanner))=[];
                y(isnan(nanner))=[];
%                 plot(time(round(x)+1),hgt(round(y)+1),'r','LineWidth',2);
%                 text(time(round(upxc(draft_count))),hgt(round(z34(draft_count))),num2str(regions(rr).MajorAxisLength),'fontsize',15,'color','k');
%                 pause(1);

            end
            %
            
            
            
        end %end of size threshold check
    end %end up updraft region loop
end

subplot(2,1,2);
xlim([startday+19/24 startday+21/24]);
datetick('x','keeplimits');
ylim([300 2500]);
xlabel('Time [UTC]')
ylabel('Height [m MSL]')
grid on; box on; set(gca,'layer','top','linewidth',2,'fontsize',15,'fontweight','bold');
cbh=colorbar; ylabel(cbh,'Draft Object ID','fontsize',18,'fontweight','bold');


% export_fig updraft_object_example.png -m2