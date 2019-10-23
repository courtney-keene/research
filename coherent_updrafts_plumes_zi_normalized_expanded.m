%% Script to compute updraft properties at different heights in the CBL
close all
clear all

daysec=24*60*60; %number of seconds in a day
load('~/Dropbox/MATLAB/Lidar_Codes/vert_vel_map'); %divergent colormap
pflag=1; % flag to toggle plotting on/off on=1

target_years={'2015' '2016'}; %select the year of data you want to process

%% Load locally forced convection dates to consider
test=dlmread('~/Dropbox/Shallow_to_Deep/Specturm_of_depths.csv','/'); %this list has various categories of convective development
dates=datenum(test(:,3)+2000,test(:,1),test(:,2)); %convert to matlab time
cudays_all=unique(dates); %make sure the list is unique

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
% % dnctime=nan(100000,1); %downdraft center times
% % dndur=nan(100000,1); %downdraft duration
% % dnchord=nan(100000,1); %downdraft chord length
% % dnwmin=nan(100000,1); %downdraft maximum velocity
% % dnwmean=nan(100000,1); %downdraft mean velocity
% % dnspd=nan(100000,1); %mean CBL wind speed for each downdraft
% % dnwstar=nan(100000,1); %CBL wstar value
% % dnz=nan(100000,1); %downdraft centroid height;
% % dnarea=nan(100000,1);%updraft area (m^2 time-> space converted using wind speed)
% % dn_draft_count=0; % keep track of the numer of updrafts in our sample
% % dnZi=nan(100000,1); %CBL height for each updraft.
% % dntilt=nan(100000,1); %orientation of updraft object

%% PROCESS LIDAR DATA

for ss=5%1:length(sites)  % loop over the 5 lidar sites at SGP
    
    [upnum,upfrac,dnfrac,upflux,Wmax,Wmin]=deal(nan(length(cudays_all),24,320)); %initialize some matrices
    cblhr=nan(length(cudays_all),24);
    for ty=1%:length(target_years); %right now this does nothing
        target_year=target_years{ty};
        
        %% load the lidar derived boundary layer horizontal wind profoile data (produced in "lidar_wind.m")
        fin=strcat('/Volumes/My Passport for Mac/SGP_LIDAR_DATA/WINDS/','lidar_winds_',sites{ss},'.mat');
        load(fin,'wspdmaster','wdmaster','timemaster','hgts')
        whgts=hgts; clear hgts; %rename some of the variables
        
        %% list the Lidar vertical velocity and backscatter data files
        dirstr=strcat('/Volumes/My Passport for Mac/SGP_LIDAR_DATA/RAW_DATA/',sites{ss},'/','*.cdf');
        files=dir(dirstr);
        
        %% idealized vectors for use later in compositing
        idealtime=-1.75:.005:1.75; % a cloud center relative time vector at arbitrary resolution
        hideal=0:.005:1.5; % a cloud base relative height vector at arbitrary resolution
        [tmati,hmati]=meshgrid(idealtime,hideal); %grid the interpolation values
        
        
        
        %     cloud_flag=nan(1500,1); %cloudy/not cloudy flag (1=cloudy)
        %     cloudB=nan(1500,1); %value of backscatter in the updraft
        %     Wmaster=nan(1500,length(hideal),length(idealtime)); %normalized cloud centered vertical velocity over the depth of the CBL
        %
        %% load day categories - > manually derived cumulus days from looking at modis and apply other criteria
        %     fin2=strcat('cu_days_manual_',target_year,'.mat');
        %     load(fin2);
        %load('cu_days_manual_2014.mat','cudays');
        
        %% load lidar derived CBL heights and Wstar values-> these can be useful for identifying when clouds are near the top of the CBL
        
        if ss==5
            load('~/Dropbox/MATLAB/SGP_LIDAR_CODES/C1_2011_2017_CBLH');
            timevec2=Time_all; %rename the time vector for the CBLH heights and Wstar data
            CBLH=Zi_all; %rename the Convective Boundary Layer height (Zi)
            WSTAR=Wstar_all; %rename the wstar data
            CBLH(abs(gradient(CBLH))>700)=nan; %eliminate large jumps in the CBL depth (spurious features)
        else
            load(strcat('/Volumes/My Passport for Mac/SGP_LIDAR_DATA/EXTEND_WSTATS_VAP/',sites{ss},'_CBLH.mat'));
            WSTAR=Wstar_all;
            CBLH=Zi_all;
            CBLH(abs(gradient(CBLH))>700)=nan; %eliminate large jumps in the CBL depth (spurious features)
            timevec2=Time_all; %rename the time vector for the CBLH heights and Wstar data
        end
        %% load and process the 1 hz lidar data
        
        for ff=1:length(files)
            
            fname=strcat(files(ff).folder,'/',files(ff).name);
            %Parse the date/time from the file name
            yyyy=str2double(fname(end-18:end-15));
            month=str2double(fname(end-14:end-13));
            day=str2double(fname(end-12:end-11));
            hour=str2double(fname(end-9:end-8));
            
            fileday=datenum(yyyy,month,day); %construct the fileday in matlab serial time
            filetime=datenum(yyyy,month,day,hour,0,0);
            cuidx=find(fileday==cudays_all);
            %if fileday==datenum(2016,6,25);
            if fileday>=datenum(2011,5,1) && ismember(fileday,cudays_all)%check to see if this day is in the list of cumulus days
                %ncdisp(fname);
                
                if hour>=17 && hour<23 % only process data for hours where ShCu are likely and the CBL growth is less extreme
                    %hour
                    B=double(ncread(fname,'attenuated_backscatter')); %attenuated backscatter coefficeint
                    B(1:3,:)=nan; %remove the first few range gate (bad data)
                    W=double(ncread(fname,'radial_velocity')); %read the Line of Site velocity (aka vertical velocity)
                    %filter W using a range and/or standard deviation filter
                    WR=rangefilt(W,ones(5,5));%determine the range of data at every grid point using a 5x5 tile of points
                    W(WR>8)=nan; %eliminate points characterized by large noise
                    SNR=double(ncread(fname,'intensity')); %signal to noise ratio
                    W(SNR<1.00)=nan; %remove W less than SNR
                    W=wiener2(W,[3 3]); %this is a noise filter in a 3x3 domain to remove "speckling" the algorithm is based on the local variance and mean
                    %W(abs(W)>7)=nan;%range check
                    if str2double(target_year)<=2012
                        time=ncread(fname,'time_offset');
                    else
                        time=ncread(fname,'time');
                    end
                    basetime=double(ncread(fname,'base_time')); %base time
                    %convert to matlab serial time
                    basetime=(basetime./(60*60*24))+datenum(1970,1,1); %add 1970 time offset and convert basetime to matlab time
                    time=(time./(60*60*24))+basetime; %combine time and basetime to create the full
                    
                    %remove points adjacent to the wind profiles (this is optional)
                    %W(:,gradient(time)>(1.2.*mode(gradient(time))))=nan;
                    
                    base_height=double(ncread(fname,'alt')); % this can be parsed from the global attributes
                    hgts=base_height+double(ncread(fname,'range')); %constructs a height vector
                    %for switch to XR system, ignore the extra range gates
                    %(changes from 9.6km to 12 km)
                    if length(hgts)>320
                        hgts=hgts(1:320);
                        W=W(1:320,:);
                        SNR=SNR(1:320,:);
                        B=B(1:320,:);
                    end
                    
                    %### find the CBL heights for this hour
                    cblhidx=find(timevec2>=datenum(yyyy,month,day,hour,5,0) & timevec2<=datenum(yyyy,month,day,hour,55,0));
                    cblhnow=nanmedian(CBLH(cblhidx));% compute the median CBL height for the hour, use this later to restrict which clouds are in the composites
                    wstarnow=nanmedian(WSTAR(cblhidx));
                    
                    %% find hourly-average CBL-averaged wind speed
                    if ~isempty(cblhnow) && ~isnan(cblhnow) %if the CBL height is defined
                        whidx=find(whgts<=cblhnow);
                    else %if now CBL height is defined assume a CBL height of 1800 m
                        whidx=find(whgts<=1800);
                    end
                    
                    wtidx=find(timemaster>=filetime & timemaster<(filetime+1/24)); %time index for the hour of interest
                    spdnow=wspdmaster(whidx,wtidx); %parse out the wind speeds
                    spdbar=nanmean(spdnow(:)); %flatten and take the mean of all values... this is admittedly a bit course but still useful
                    
                    %### create a meshed time/height grid from the lidar range
                    %gates and time stamps
                    [tmat,hmat]=meshgrid(time,hgts);
                    
                    %Define a cloud backscatter threshold (.25*10^-4 before 2015,
                    %.8*10^-4 for 2015 and later).
                    if yyyy<2015
                        cloud_threshold=.25.*10^-4;
                    else
                        cloud_threshold=.8*10^-4;
                    end
                    
                    %W(B>cloud_threshold)=nan;
                    
                    %% Optional plotting of the data for inspection purposes
                    if pflag==1
                        f=figure(1);set(f,'color','w','pos',[10 100 700 500]);clf
                        subplot(2,1,1)
                        hold on;
                        colormap(jet);
                        %pcolor(time,hgts,(real(log10(B))));shading flat; caxis([-6.5 -4.5]);
                        pcolor(time,hgts,B);shading flat;caxis([0 1]*10^-4);
                        ylim([base_height 3500]);
                        xlim([time(1) time(end)]); datetick('x','keeplimits');
                        cbh=colorbar;
                        ylabel(cbh,'m^{-1} sr^{-1}','fontsize',15,'fontweight','bold');
                        xlabel('Time [UTC]','fontsize',15,'fontweight','bold');
                        ylabel('Height [m]','fontsize',15,'fontweight','bold');
                        title(datestr(nanmean(time)));
                        sp2=subplot(2,1,2);
                        rbm=rbmapper(5,-5);
                        colormap(sp2,rbm);
                        pcolor(time,hgts,medfilt2(W,[1 1]));shading flat; caxis([-3 3]);
                        hold on;
                        plot(timevec2(cblhidx),CBLH(cblhidx),'--k','linewidth',2);
                        
                        plot(timevec2(cblhidx),CBLH(cblhidx)-100,'--k','linewidth',2);
                        plot([timevec2(cblhidx(1)) timevec2(cblhidx(end))],[cblhnow-500 cblhnow-500],'--c','linewidth',2);
                        %scatter(tmat(:),hmat(:),10,'ow','filled','markeredgecolor','k'); %this adds markers indicative of clouds
                        ylim([base_height 3500]);
                        xlim([time(1) time(end)]); datetick('x','keeplimits');
                        xlabel('Time [UTC]','fontsize',15,'fontweight','bold');
                        ylabel('Height [m]','fontsize',15,'fontweight','bold');
                        cbh=colorbar;
                        ylabel(cbh,'m s^{-1}','fontsize',15,'fontweight','bold');
                        
                    end
                    
                    %% USE IMAGE PROCESSING TO LOCATE TEMPORALLY AND SPATIALLY COHERENT UPDRAFT OBJECTS
                    %### create a binary mask for vertical velocity greater
                    %than some threshold vertical velocity:
                    wthreshold=.5; %define the updraft threshold value
                    
                    Wmask=W; Wmask(W<wthreshold)=false; Wmask(W>=wthreshold)=true; %create true/false mask for points on either side of the threshold
                    %set edges to false
                    Wmask(1,:) = false;
                    Wmask(end,:) = false;
                    Wmask(:,1) = false;
                    Wmask(:,end) = false;
                    Wmask(isnan(Wmask))=false;
                    Wmask=logical(Wmask); %make sure the whole array is logical true/false, which is the same as white/black
                    
                    if pflag==1
                        figure(100);clf;
                        sp1=subplot(3,1,1);
                        rbm=rbmapper(5,-5);
                        colormap(sp1,rbm);
                        pcolor(time,hgts,medfilt2(W,[1 1]));shading flat; caxis([-3 3]);
                        hold on;
%                         sp2=subplot(3,1,2);
%                         
%                         pcolor(time,hgts,Wmask);shading flat;
                    end
                    
                    size_threshold=50;
                    Wmaskfilt=bwareafilt(Wmask,[size_threshold inf]); %filter on an area of "size_threshold" pixels (what remains we will call "coherent" updrafts)
                    updrafts=bwlabel(Wmaskfilt); %label each independent connected updraft region
                    upfrac(cuidx,hour,:)=nansum(Wmaskfilt,2)./size(Wmaskfilt,2); %this is the coherent updraft fraction for the hour as a function of height
                    cblhr(cuidx,hour)=cblhnow; %this stores the CBL height for each hour
                    Wcopy=W; Wcopy(Wmaskfilt==0)=nan; %now use the updraft mask to blank downdraft values for
                    Wmax(cuidx,hour,:)=nanmax(Wcopy,[],2); %maximum vertical velocity for the hour;
                    
                    % Leverage MATLAB regionprop algorithm to compute
                    % properties of each updraft object (kind of a black
                    % box)
                    regions=regionprops(updrafts,'Area','Extrema','BoundingBox','Centroid','PixelList','PixelIdxList','Orientation','convexhull');
                    
                    %## Optional plotting to show the unique updraft
                    %objects
                    if pflag==1
                        sp3=subplot(2,1,2);
                        colormap(sp3,colorcube(200));
                        pcolor(time,hgts,updrafts(:,:));shading flat;%
                    end
                    return
                    
                    %Parse out various properties
                    centers=[regions.Centroid];
                    areas=[regions.Area];
                    
                    % loop over each updraft object (contained in regions)
                    if numel(regions)>0
                        for rr=1:length(regions)
                            if regions(rr).Area>size_threshold %must be at least 100 pixels in the area
                                draft_count=draft_count+1; %count the number of coherent updrafts
                                if draft_count<53727
                                    idxnow=regions(rr).PixelList; %list of pixels
                                    %scatter(idxnow(:,1),idxnow(:,2),'*k');
                                    %regions(rr).Wmax=nanmax(W(regions(rr).PixelIdxList));
                                    upwmax(draft_count)=nanmax(W(regions(rr).PixelIdxList)); %find the maximum vertical velocity value in the updraft
                                    upwmean(draft_count)=nanmean(W(regions(rr).PixelIdxList)); %compute the mean vertical velocity in the updraft
                                    upwstar(draft_count)=wstarnow; %assign each updraft a convective velocity scale based on the hourly values computed previously (this is just for reference and mormalization and not a property of the updraft)
                                    upspd(draft_count)=spdbar; %mean CBL wind speed
                                    pixel_area=(spdbar*1).*mode(gradient(hgts));   %this is the "size" of a given pixek = speed*1 second (meters) x gate length (30 m)
                                    uparea(draft_count)=regions(rr).Area.*pixel_area; %the total updraft area is therefore the #of pixels (regions(rr).Area) multiplied by the "pixel size"
                                    upZi(draft_count)=cblhnow; % cbl height corresponding to each updraft (not locally computed).
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
                                    
%                                     if xcenter>250 && xcenter<(size(W,2)-250) %this keeps it away from the edges of the time domain...
%                                         if (round((zcenter+ztop)/2))<100 %this restricts the scenes to updrafts confined to the CBL
%                                             wsnapshot(draft_count,:,:)=W(1:100,(xcenter-250):(xcenter+250)); %this pulls the data from 250 second before-> 250 seconds after the center of the updraft
%                                             wxz34(draft_count)=nanmax(wsnapshot(draft_count,round((zcenter+ztop)/2),:));
%                                         end
%                                     end
                                    
                                    upxc(draft_count)=LC+(RC-LC)/2;%xcenter; %xcenter of the updraft at 3/4 height
                                    upzc(draft_count)=hgts(zcenter); %zcenter of the updraft
                                    upzx(draft_count)=hgts(round(ztop)); %top of the updraft (x here indicates maximum height)
                                    upctime(draft_count)=time(round(upxc(draft_count)));
                                    
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
                                    z34=round((ztop+zcenter)/2);
                                    if pflag==1
                                        plot(regions(rr).ConvexHull(:,1),regions(rr).ConvexHull(:,2),'k');
                                        plot(regions(rr).Centroid(1), regions(rr).Centroid(2),'*w','markersize',5);
                                        plot([regions(rr).Centroid(1)-width/2 regions(rr).Centroid(1)+width/2 ], [regions(rr).Centroid(2) regions(rr).Centroid(2)],'--k');
                                        plot([LC RC],[zcenter zcenter],'--r');
                                        plot([LC RC],[z34 z34],'--c');
                                        plot(regions(rr).ConvexHull(zidx,1), regions(rr).ConvexHull(zidx,2),'*r','markersize',5);
                                        ylim([0 120]);
                                        pause(1);
                                    end
                                    %
                                elseif draft_count>53727
                                    idxnow=regions(rr).PixelList; %list of pixels
                                    dcn=draft_count-53727;
                                    %scatter(idxnow(:,1),idxnow(:,2),'*k');
                                    %regions(rr).Wmax=nanmax(W(regions(rr).PixelIdxList));
                                    upwmax2(dcn)=nanmax(W(regions(rr).PixelIdxList)); %find the maximum vertical velocity value in the updraft
                                    upwmean2(dcn)=nanmean(W(regions(rr).PixelIdxList)); %compute the mean vertical velocity in the updraft
                                    upwstar2(dcn)=wstarnow; %assign each updraft a convective velocity scale based on the hourly values computed previously (this is just for reference and mormalization and not a property of the updraft)
                                    upspd2(dcn)=spdbar; %mean CBL wind speed
                                    pixel_area=(spdbar*1).*mode(gradient(hgts));   %this is the "size" of a given pixek = speed*1 second (meters) x gate length (30 m)
                                    uparea2(dcn)=regions(rr).Area.*pixel_area; %the total updraft area is therefore the #of pixels (regions(rr).Area) multiplied by the "pixel size"
                                    upZi2(dcn)=cblhnow; % cbl height corresponding to each updraft (not locally computed).
                                    upz2(dcn)=hgts(round(regions(rr).Centroid(2))); %the height of the center of the updraft
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
                                    
                                    CWIDTH2(dcn)=RC-LC; %width of the updraft at the upper quarter location
                                    
                                    %                                     if xcenter>250 && xcenter<(size(W,2)-250) %this keeps it away from the edges of the time domain...
                                    %                                         if (round((zcenter+ztop)/2))<100 %this restricts the scenes to updrafts confined to the CBL
                                    %                                             wsnapshot(draft_count,:,:)=W(1:100,(xcenter-250):(xcenter+250)); %this pulls the data from 250 second before-> 250 seconds after the center of the updraft
                                    %                                             wxz34(draft_count)=nanmax(wsnapshot(draft_count,round((zcenter+ztop)/2),:));
                                    %                                         end
                                    %                                     end
                                    
                                    upxc2(dcn)=LC+(RC-LC)/2;%xcenter; %xcenter of the updraft at 3/4 height
                                    upzc2(dcn)=hgts(zcenter); %zcenter of the updraft
                                    upzx2(dcn)=hgts(round(ztop)); %top of the updraft (x here indicates maximum height)
                                    upctime2(dcn)=time(round(upxc2(dcn)));

                                    
                                    %% normalized snap shot (normalized in time and in CBL height (Zi))
                                    if isfinite(CWIDTH2(dcn)) && CWIDTH2(dcn)>0
                                        zvec=(hgts(1:100));
                                        if isfinite(upZi2(dcn))
                                            znorm=zvec./upZi2(dcn);
                                        else
                                            znorm=zvec./2000;
                                        end
                                        xnorm=((1:length(time))-xcenter)./CWIDTH2(dcn); %this is a time (x-data) normalization centered on the updraft center period
                                        xideal=-2:.025:2;
                                        zideal=0:.025:1.5;
                                        [xmat,zmat]=meshgrid(xnorm,znorm);
                                        [ximat,zimat]=meshgrid(xideal,zideal);
                                        Wnorm2(dcn,:,:)=interp2(xmat,zmat,W(1:100,:),ximat,zimat);
                                    end
                                    
                                    %%
                                    z34=round((ztop+zcenter)/2);
                                    if pflag==1
                                        plot(regions(rr).ConvexHull(:,1),regions(rr).ConvexHull(:,2),'k');
                                        plot(regions(rr).Centroid(1), regions(rr).Centroid(2),'*w','markersize',5);
                                        plot([regions(rr).Centroid(1)-width/2 regions(rr).Centroid(1)+width/2 ], [regions(rr).Centroid(2) regions(rr).Centroid(2)],'--k');
                                        plot([LC RC],[zcenter zcenter],'--r');
                                        plot([LC RC],[z34 z34],'--c');
                                        plot(regions(rr).ConvexHull(zidx,1), regions(rr).ConvexHull(zidx,2),'*r','markersize',5);
                                        ylim([0 120]);
                                        pause(1);
                                    end
                                end %end of excess updraft loop
                            end %end of size threshold check
                        end %end up updraft region loop
                        
                    end %end of region check
                    
                end %end hour loop
            end % end check against CU day list
            disp(strcat('Draft Count=',num2str(draft_count)));
            
        end %end file loop
        
    end %end year loop
end
return
%%
% save updraft object output data
save('updraft_objects_20190227.mat','upzc','upzx','upxc','upZi','upctime','CWIDTH','upspd','upwmax','upwstar','upwmean','uparea','upzc2','upzx2','upxc2','upZi2','CWIDTH2','upspd2','upwmax2','upwstar2','upctime2','upwmean2','uparea2','xideal','zideal','-v7.3');
save('updraft_wnorm_20190227.mat','Wnorm','Wnorm2','xideal','zideal','-v7.3');

%% Compute upper 3/4 updraft location
upz_norm=(0.5.*(upzc+upzx))./upZi;
upzxn=upzx./upZi;
upz_norm2=(0.5.*(upzc2+upzx2))./upZi2;
upzxn2=upzx2./upZi2;
qtiles=[.25 .45 .65 .85 1.05 1.25]

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
    text(800,.3,strcat('Mean=',num2str(mchord(ii))),'fontsize',12,'fontweight','bold');
    text(800,.25,strcat('Median=',num2str(mdchord(ii))),'fontsize',12,'fontweight','bold');
    set(gca,'fontsize',12,'fontweight','bold','linewidth',2,'layer','top');
    
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
    text(1,.35,strcat('Mean=',num2str(WM(ii))),'fontsize',12,'fontweight','bold');
    text(1,.3,strcat('Median=',num2str(WMD(ii))),'fontsize',12,'fontweight','bold');
    set(gca,'fontsize',12,'fontweight','bold','linewidth',2,'layer','top');
end

%% Normalized
for ii=1:length(qtiles)-1
    figure(3);
    idxnow=[idx(ii).locs; idx2(ii).locs+53727];
    cmap=rbmapper_coarse(1.25,-.3);
    colormap(cmap);
    subplot(3,length(qtiles)-1,ii);
    WBAR1=squeeze(nanmean(Wnorm(idx(ii).locs,:,:),1));
    WBAR2=squeeze(nanmean(Wnorm2(idx2(ii).locs,:,:),1))
    WBAR=(WBAR1+WBAR2)./2;
%     pcolor(ximat.*mdchord(ii)/2,zimat,WBAR);shading flat; caxis([-.3 1.25]);
    contourf(ximat.*mdchord(ii)/2,zimat,WBAR,[-5:.1:5],'linestyle','none');shading flat; caxis([-.3 1.25]);

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