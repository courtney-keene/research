%% Script to compute cloud statistics for individual cumulus clouds at ARM SGP
close all
clear all

daysec=24*60*60; %number of seconds in a day
% load('~/Dropbox/MATLAB/Lidar_Codes/vert_vel_map'); %divergent colormap
pflag=1; % flag to toggle plotting on/off on=1

target_years={'2015' '2016'}; %select the year of data you want to process

% index=dlmread('../shallow_cumulus_days_2016_2017_lidar_check.csv',',',1,0);
% dates=datenum(index(:,1),index(:,2),index(:,3));
% gb=index(:,end);
% dates(gb<4)=[];
% cudays_all=(datenum(dates));

% addpath('~/Dropbox/MATLAB/SGP_LIDAR_CODES/');
index=readtable('/Users/courtneykeene/Desktop/Research/QC_water_vapor_flux_days.txt'); %SHCU DAYS DURING PERIOD OF LIDAR OVERLAP
dates=index.Var1;
cudays=datenum(dates);
cudays_all=[cudays; datenum(2016,8,18)];

sites={'E32', 'E37', 'E39', 'E41', 'C1'};

%% idealized vectors for use later in compositing
idealtime=-1.75:.01:1.75; % a cloud center relative time vector at arbitrary resolution
hideal=0:.01:1.5; % a cloud base relative height vector at arbitrary resolution
[tmati,hmati]=meshgrid(idealtime,hideal); %grid the interpolation values
%%

for ss=1:5%:length(sites)
    [upnum,upfrac,upflux]=deal(nan(length(cudays_all),24,320));
    cblhr=nan(length(cudays_all),24);
    for ty=1%:length(target_years);
        target_year=target_years{ty};
        %% load the lidar boundary layer wind data (produced in "lidar_wind.m")
        fin=strcat('/Users/courtneykeene/Desktop/Research/lidar_data/SGP_LIDAR_DATA/WINDS/','lidar_winds_',sites{ss},'.mat');
        load(fin,'wspdmaster','wdmaster','timemaster','hgts')
        whgts=hgts; clear hgts; %rename some of the variables
        
        %% list the Lidar vertical velocity and backscatter data files
        dirstr=strcat('/Users/courtneykeene/Desktop/Research/lidar_data/SGP_LIDAR_DATA/RAW_DATA/',sites{ss},'/','*.cdf');
        files=dir(dirstr);
        
        %% idealized vectors for use later in compositing
        idealtime=-1.75:.005:1.75; % a cloud center relative time vector at arbitrary resolution
        hideal=0:.005:1.5; % a cloud base relative height vector at arbitrary resolution
        [tmati,hmati]=meshgrid(idealtime,hideal); %grid the interpolation values
        
        %% Initialize some empty matrices for compositing
        
        upctime=nan(10000,1); %updraft center times
        updur=nan(10000,1); %updraft duration
        upchord=nan(10000,1); %updraft chord length
        upwmax=nan(10000,1); %updraft maximum velocity
        upwmean=nan(10000,1); %updraft mean velocity
        upspd=nan(10000,1); %mean CBL wind speed for each updraft
        upwstar=nan(10000,1); %CBL wstar value
        upz=nan(10000,1); %updraft centroid height;
        uparea=nan(10000,1);%updraft area (m^2 time-> space converted using wind speed)
        draft_count=0; % keep track of the numer of updrafts in our sample
        upZi=nan(10000,1); %CBL height for each updraft.
        uptilt=nan(10000,1); %orientation of updraft object
        %     cloud_flag=nan(1500,1); %cloudy/not cloudy flag (1=cloudy)
        %     cloudB=nan(1500,1); %value of backscatter in the updraft
        %     Wmaster=nan(1500,length(hideal),length(idealtime)); %normalized cloud centered vertical velocity over the depth of the CBL
        %
        %% load day categories - > manually derived cumulus days from looking at modis and apply other criteria
        %     fin2=strcat('cu_days_manual_',target_year,'.mat');
        %     load(fin2);
        %load('cu_days_manual_2014.mat','cudays');
        
        %% load lidar derived CBL heights and Wstar values-> these can be useful for identifying when clouds are near the top of the CBL
        load(strcat('/Users/courtneykeene/Desktop/Research/lidar_data/SGP_LIDAR_DATA/EXTEND_WSTATS_VAP/',sites{ss},'_CBLH.mat'));
        timevec2=Time_all;CBLH=Zi_all;WSTAR=Wstar_all;
        CBLH(abs(gradient(CBLH))>700)=nan; %eliminate large jumps in the CBL depth (spurious features)
        
        
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
            if ismember(fileday,cudays_all)%check to see if this day is in the list of cumulus days
                %ncdisp(fname);
                
                if hour>=12 && hour<24 % only process data for hours where ShCu are likely
                    %hour
                    B=double(ncread(fname,'attenuated_backscatter'));
                    B(1:3,:)=nan; %remove the first few range gate (bad data)
                    W=double(ncread(fname,'radial_velocity'));
                    %filter W using a range and/or standard deviation filter
                    WR=rangefilt(W,ones(5,5));%determine the range of data at every grid point
                    W(WR>8)=nan; %eliminate points characterized by large noise
                    SNR=double(ncread(fname,'intensity'));
                    W(SNR<1.00)=nan;
                    W=wiener2(W,[3 3]);
                    %W(abs(W)>7)=nan;%range check
                    if str2double(target_year)<=2012
                        time=ncread(fname,'time_offset');
                    else
                        time=ncread(fname,'time');
                    end
                    basetime=double(ncread(fname,'base_time'));
                    %convert to matlab serial time
                    basetime=(basetime./(60*60*24))+datenum(1970,1,1);
                    time=(time./(60*60*24))+basetime;
                    
                    %remove points adjacent to the wind profiles (this is optional)
                    %W(:,gradient(time)>(1.2.*mode(gradient(time))))=nan;
                    
                    base_height=double(ncread(fname,'alt')); % this can be parsed from the global attributes
                    hgts=base_height+double(ncread(fname,'range')); %constructs a height vector
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
                    %% find hourly average CBL averaged wind speed
                    whidx=find(whgts<=1800);
                    wtidx=find(timemaster>=filetime & timemaster<(filetime+1/24));
                    spdnow=wspdmaster(whidx,wtidx);
                    spdbar=nanmean(spdnow(:));
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
                    
                    %% Optional plotting of the data
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
%                         rbm=rbmapper(5,-5);
                        colormap(sp2,jet(12));
                        pcolor(time,hgts,medfilt2(W,[1 1]));shading flat; caxis([-3 3]);
                        hold on;
%                         plot(timevec2(cblhidx),CBLH(cblhidx),'--k','linewidth',2);
                        
%                         plot(timevec2(cblhidx),CBLH(cblhidx)-100,'--k','linewidth',2);
%                         plot([timevec2(cblhidx(1)) timevec2(cblhidx(end))],[cblhnow-500 cblhnow-500],'--c','linewidth',2);
                        %scatter(tmat(:),hmat(:),10,'ow','filled','markeredgecolor','k'); %this adds markers indicative of clouds
                        ylim([base_height 3500]);
                        xlim([time(1) time(end)]); datetick('x','keeplimits');
                        xlabel('Time [UTC]','fontsize',15,'fontweight','bold');
                        ylabel('Height [m]','fontsize',15,'fontweight','bold');
                        cbh=colorbar;
                        ylabel(cbh,'m s^{-1}','fontsize',15,'fontweight','bold');
                        
                    end
                    
                    
                    
                    
                    %% Coherent updrafts
                    
                    %figure; sp1=subplot(2,1,1);pcolor(W);shading flat;
                    %colormap(sp1,vert_vel_map);
                    %caxis([-5 5]);
                    %Wmask=im2bw(W,.1);% create a binary BW image from the W data, threshold is 0.7 m/s for updrafts
                    
                    Wmask=W; Wmask(W<0.5)=false; Wmask(W>=0.5)=true;
                    %set edges to false
                    Wmask(1,:) = false;
                    Wmask(end,:) = false;
                    Wmask(:,1) = false;
                    Wmask(:,end) = false;
                    Wmask(isnan(Wmask))=false;
                    Wmask=logical(Wmask);
                    
                    %subplot(2,1,2);
                    %pcolor(Wmask);shading flat;
                    Wmaskfilt=bwareafilt(Wmask,[100 inf]); %filter on an area of 100 pixels (what remains we will call "coherent" updrafts)
                    updrafts=bwlabel(Wmaskfilt); %label each updraft
                    upfrac(cuidx,hour,:)=nansum(Wmaskfilt,2)./size(Wmaskfilt,2); %this is the coherent updraft fraction for the hour
                    cblhr(cuidx,hour)=cblhnow;
                    regions=regionprops(updrafts,'Area','Centroid','PixelList','PixelIdxList','Orientation')%,'WeightedCentroid');
                    %figure(80);clf; colormap(colorcube(200));
                    %pcolor(updrafts(:,:));shading flat;%
                    centers=[regions.Centroid];
                    
                    areas=[regions.Area];
                    %ylim([0 80]);
                    %hold on;
                    if numel(regions)>0
                        for rr=1:length(regions)
                            if regions(rr).Area>100 %must be at least 100 pixels in the area
                                draft_count=draft_count+1;
                                idxnow=regions(rr).PixelList;
                                %scatter(idxnow(:,1),idxnow(:,2),'*k');
                                %regions(rr).Wmax=nanmax(W(regions(rr).PixelIdxList));
                                upwmax(draft_count)=nanmax(W(regions(rr).PixelIdxList));
                                upwmean(draft_count)=nanmean(W(regions(rr).PixelIdxList));
                                upwstar(draft_count)=wstarnow;
                                upspd(draft_count)=spdbar; %mean CBL wind speed
                                pixel_area=(spdbar*1).*mode(gradient(hgts));   %speed*1 second (meters) x gate length (30 m)
                                uparea(draft_count)=regions(rr).Area.*pixel_area;
                                upZi(draft_count)=cblhnow;
                                uptilt(draft_count)=regions(rr).Orientation;
                                upz(draft_count)=hgts(round(regions(rr).Centroid(2)));
                            else
                                %regions(rr).Wmax=nan;
                            end
                        end
                        %hold on
                        
                        %                     wmax=[regions.Wmax];
                        %                     figure(88);hold on;
                        %                     scatter(areas,wmax);
                        %                     pause(.2);
                    end
                    
                    
                    
                    
                end %end hour loop
                
            end % end check against CU day list
        end %end file loop
        fout=strcat('lidar_updraft_objects_extended_',sites{ss},'.mat');
        fout2=strcat('lidar_updraft_fraction_',sites{ss},'.mat');
        save(fout,'upspd','upwmax','upwmean','uparea','upz','uptilt','upwstar');
        save(fout2,'upfrac','cblhr','hgts');
    end %end year loop
end
%%


%%
figure(200);clf;
tvec=datenum(yyyy,month,day):1/24:(datenum(yyyy,month,day,23,0,0));
AGL=[1:size(upfrac,3)].*30;
upfrac_mean=squeeze(nanmean(upfrac,1));
Zi=nanmean(cblhr,1);
pcolor(tvec,AGL,upfrac_mean');shading flat;
ylim([0 3500]);
hold on; 
plot(tvec,Zi,'--k','linewidth',2);
xlim([datenum(yyyy,month,day)+12/24 datenum(yyyy,month,day)+1]);
datetick('x','HH:MM','keeplimits');

%% bin updraft area and updraft vertical velocity
%%
fins=dir('/Users/courtneykeene/Desktop/Research/lidar_updraft_objects_extended*.mat');
UPAREA=[];
UPWX=[];
WS=[];
for fff=1:length(fins)
    load(strcat(fins(fff).folder,'/',fins(fff).name));
    UPAREA=cat(1,UPAREA,uparea);
    UPWX=cat(1,UPWX,upwmax);
    WS=cat(1,WS,upwstar);
end

bin_width=(.15*10^6);
bins=0:bin_width:(5*10^6);
f=figure(3000);clf;hold on;set(f,'color','w');
for bb=2:length(bins)
    wstarbin(bb)=nanmedian(WS(UPAREA>=bins(bb-1) & UPAREA<bins(bb)));
    wbin(bb)=nanmedian(UPWX(UPAREA>=bins(bb-1) & UPAREA<bins(bb)));
    wbin75(bb)=prctile(UPWX(UPAREA>=bins(bb-1) & UPAREA<bins(bb)),75);
    wbin95(bb)=prctile(UPWX(UPAREA>=bins(bb-1) & UPAREA<bins(bb)),95);
    wbin5(bb)=prctile(UPWX(UPAREA>=bins(bb-1) & UPAREA<bins(bb)),5);
    
    wbin25(bb)=prctile(UPWX(UPAREA>=bins(bb-1) & UPAREA<bins(bb)),25);
    %    wbin(bb)=nanmedian(upwmax(uparea>=bins(bb-1) & uparea<bins(bb)));
    %    wbin75(bb)=prctile(upwmax(uparea>=bins(bb-1) & uparea<bins(bb)),75);
    %    wbin95(bb)=prctile(upwmax(uparea>=bins(bb-1) & uparea<bins(bb)),95);
    %    wbin5(bb)=prctile(upwmax(uparea>=bins(bb-1) & uparea<bins(bb)),5);
    %
    %    wbin25(bb)=prctile(upwmax(uparea>=bins(bb-1) & uparea<bins(bb)),25);
    plot([bins(bb) bins(bb)],[wbin25(bb) wbin75(bb)]./wstarbin(bb),'-','linewidth',4,'color',[1 .7 .7]);
    plot([bins(bb) bins(bb)],[wbin5(bb) wbin95(bb)]./wstarbin(bb),'-','linewidth',1,'color',[1 .7 .7]);
    
    scatter(bins(bb),wbin(bb)./wstarbin(bb),30,'ok','filled');
end
xlabel('Updraft Area [km^2]','fontsize',15,'fontweight','bold');
grid on; box on; set(gca,'linewidth',2,'fontsize',15,'fontweight','bold','layer','top');
xlim([0 3]*10^6);
ticks=get(gca,'xtick');
set(gca,'xticklabel',ticks./10^6);
ylabel('updraft [m s^{-1}]','fontsize',15,'fontweight','bold');