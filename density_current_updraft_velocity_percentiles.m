%% Density Current Updraft percentiles

close all
clear all

addpath('~/Dropbox/MATLAB/');
daysec=24*60*60;
pflag=1
dates=csvread('~/Dropbox/density_current_list2.csv'); %this is a list of year, month, day, hour, minutes for density current passage times the last column also has the speed
gtimes=datenum(dates(:,1),dates(:,2),dates(:,3),dates(:,4),dates(:,5),0*dates(:,5)); %converts to matlab serial time
gspeed=dates(:,end); %parse out the speed of each density current

aa=0; %initialize a counter

for gg=1:length(gtimes);%Loop over the list of density current start times 
    aa=aa+1; %keep count
    target_time=gtimes(gg);%parsae out the current density current start time
    speed=gspeed(gg); %current density current speed
    
    %% Define start end times for gust front
    starttime=(target_time);% start time for density current
    endtime=starttime+(1)/24; %define end time (for now) as one hour later
    %% List the lidar files spanning the start and end times (note there should be one file for every hour
    f1=dir(strcat('/Volumes/My Passport for Mac/SGP_LIDAR_DATA/RAW_DATA/C1/sgpdlfptC1.b1.',datestr(target_time,'YYYYmmdd'),'.',datestr(target_time,'HH'),'*.cdf'));
    f2=dir(strcat('/Volumes/My Passport for Mac/SGP_LIDAR_DATA/RAW_DATA/C1/sgpdlfptC1.b1.',datestr(endtime,'YYYYmmdd'),'.',datestr(endtime,'HH'),'*.cdf'));
    files=[f1 f2]; %concatenate the list of lidar files
    fname=strcat(files(1).folder,'/',files(1).name); %create the path to the files 
%     yyyy=str2double(fname(end-18:end-15)); %parse out year
%     month=str2double(fname(end-14:end-13)); %parse out month
%     day=str2double(fname(end-12:end-11)); %parse out day
%     hour=str2double(fname(end-9:end-8));  
%     fileday=datenum(yyyy,month,day);
%     filetime=datenum(yyyy,month,day,hour,0,0);
    
    ncdisp(fname); %display file header contents
    B=double(ncread(fname,'attenuated_backscatter')); %read the attenuated backscatter
    W=double(ncread(fname,'radial_velocity')); %read the vertical velocity
    if size(W,2)>1 %provided there is some actual data in here (size check)
        
        SNR=double(ncread(fname,'intensity'));
        WR=rangefilt(W,ones(5,5));
        W=wiener2(W,[3 3]); %noise filter the velocity data
        
        W(WR>10)=nan; %remove large noisy points
        W(SNR<1.002)=nan; %remove low SNR data
        W(1:3,:)=nan; %blank the first 3 range gates which don't have good data

        %deal with times in the files
        if target_time<=datenum(2012,5,1)
            time=ncread(fname,'time_offset');
        else
            time=ncread(fname,'time');
        end
        basetime=double(ncread(fname,'base_time'));
        %convert to matlab serial time
        basetime=(basetime./(60*60*24))+datenum(1970,1,1);
        time=(time./(60*60*24))+basetime;
        dz=30; %gate length
        base_height=ncread(fname,'alt');%318.000000; % this can be parsed from the global attributes
        hgts=base_height+ dz./2 + dz.*(1:size(B,1)); %constructs a height vector
    end
    %% load second file and concatenate
    fname=strcat(files(2).folder,'/',files(2).name);
    ncdisp(fname);
    B2=double(ncread(fname,'attenuated_backscatter'));
    W2=double(ncread(fname,'radial_velocity'));
    if size(W2,2)>1
      
        SNR2=double(ncread(fname,'intensity'));
        WR=rangefilt(W2,ones(5,5));
        W2=wiener2(W2,[5 5]);
        
        W2(WR>10)=nan;

        W2(SNR2<1.002)=nan;
        W2(1:3,:)=nan;

        if target_time<=datenum(2012,5,1)
            time2=ncread(fname,'time_offset');
        else
            time2=ncread(fname,'time');
        end
        basetime2=double(ncread(fname,'base_time'));
        %convert to matlab serial time
        basetime2=(basetime2./(60*60*24))+datenum(1970,1,1);
        time2=(time2./(60*60*24))+basetime2;
        
        dz=30; %gate length
        base_height=318.000000; % this can be parsed from the global attributes
        hgts2=base_height+ dz./2 + dz.*(1:size(B2,1)); %constructs a height vector
    end
    
    %% concatenate all the lidar data
    Wall=cat(2,W,W2);
    Ball=cat(2,B,B2);
    timeall=cat(1,time,time2);
    hgtsall=hgts;
    %%
    if pflag==1
            f=figure(1);clf
            set(f,'color','w','pos',[200 200 1400 700])
            sp1=subplot(1,4,1:3);cla;
            hold on;
            cmap=rbmapper(8,-5);colormap(sp1,cmap);
            pcolor(timeall,hgtsall,(Wall));shading flat; caxis([-5 8]);cbh1=colorbar;
            ylabel(cbh1,'m s^{-1}','fontsize',15,'fontweight','bold');
            ylim([base_height 3500]);
            datetick('x')
            grid on;
            box on;
            set(gca,'linewidth',2,'layer','top','fontsize',15,'fontweight','bold');
            xlabel('Time [UTC]','fontsize',15,'fontweight','bold');ylabel('Height [m]','fontsize',15,'fontweight','bold');
            title(datestr(starttime));
            xlim([starttime-.5/24 starttime+1/24]);
            datetick('x','keeplimits')
        
%     
    end
    %% UPDRAFT PERCENTILES
    wcopy=Wall; %make a copy of the vertical velocity
    nanmask=Wall; %now make a mask to keep track of what percent of the data are nans
    nanmask(~isnan(nanmask))=1; %if it isn't a nan set to 1
    nanmask(isnan(nanmask))=0; %if it IS a nan set to 0
    nprctile=1-(nansum(nanmask,2)./size(nanmask,2)); %compute the perecnt of time at a given height that the data are populated by nans
    npmat=repmat(nprctile,1,size(nanmask,2));
    wcopy(npmat>.5)=nan; %remove data from the processing if it has more than 50% of the time sereies at that height as nans
    wvar=nanvar(wcopy,[],2); %this is the vertical velocity variance
    uptiles=prctile(wcopy,[5 25 50 75 95],2); %now compute the percentiles (as a function of height) of vertical velocity
    bigups(gg,:)=uptiles(:,5); %save the 95th percentile data (i.e., strong updraft profile) for each density current
    bigdowns(gg,:)=uptiles(:,1); %save the 5th percentil data (i.e., strong downdrafts).
    if pflag==1 %optional plotting for each density current
        figure(1);
        subplot(1,4,4);
        hold on;
        
        plot(uptiles(:,1),hgts2);hold on;
        plot(uptiles(:,2),hgts2);hold on;
        plot(uptiles(:,3),hgts2);hold on;
        plot(uptiles(:,4),hgts2);hold on;
        plot(uptiles(:,5),hgts2);hold on;
        ylim([base_height 2000]);
    end

end

%%now make a plot showing the distribution of the 95th and 5th percentile data (i.e. the 95th precentile of the 95th percentile data!)

figure(22);clf;hold on;
plot(prctile(bigups,50,1),hgts,'--k','linewidth',3)
plot(prctile(bigups,75,1),hgts,'--k','linewidth',3)
plot(prctile(bigups,90,1),hgts,'--k','linewidth',3)
plot(prctile(bigups,99,1),hgts,'--k','linewidth',3)

plot(prctile(bigdowns,50,1),hgts,'--b','linewidth',3)
plot(prctile(bigdowns,25,1),hgts,'--b','linewidth',3)

plot(prctile(bigdowns,10,1),hgts,'--b','linewidth',3)
plot(prctile(bigdowns,1,1),hgts,'--b','linewidth',3)
xlim([-12 12]);
ylim([base_height 3000]);
