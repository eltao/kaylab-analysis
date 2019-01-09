% Odor Presentation analysis
% Presented Hexanone, EMB, Geraniol, Limonene, Citronellal
% 10 presentations per odor block and each odor presented 1x various order
% (pre infusion) and 2x various order (post infusion)
% Channel map
% NOTE: RL93 uses a 9pin headstage
% -----------------------------%
%       MCS Channel arrangement       %Vdata
% L OB electrode :       8,7           6,5
% R OB electrode:        6,5           4,3
% PC:                    4,3           2,1
% -----------------------------%
% Author: Emily Tao
% This program adapts/revises analysis code by Shane Peace, 2018. 
% Changes: - Removed optogenetic analysis portions,
% - Created struct-based variable conventions for data storage and looping, 
% - General restructuring and organization


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars

% add analysis functions in path
addpath('Z:\Code.Repository\MATLAB\AnalysisFunctions')
addpath('Z:\Code.Repository\MATLAB\MultiChan import')
addpath(genpath('Z:\Code.Repository\MATLAB\chronux_2_10'))
addpath(genpath('C:\Program Files (x86)\Multi Channel Systems\MC_Rack\MCStreamSupport\Matlab\meatools\mcintfac'))

%% Load and initialize data
URL = 'Z:\Personal.Folders\Nora and Vivian\RL93\RL93_OdorOrder.xlsx'; %excel file for stim codes
uiimport(URL);

%% Load and initialize the rest of the data

% session info - change every time
Rname = 'RL93';
date = '121718';
drugCondition = 'HighAnt';
recPhase = 'Post';
% end session info

Nsess = [date '_' drugCondition]; % name session for saving files
dfile = [Rname '_' date '_odorpres_' drugCondition '_' recPhase '.mcd'];
direc = ['Z:\Personal.Folders\Nora and Vivian\' Rname '\' Nsess '\']; % directory where data is stored

chan_names = {'LOB1','LOB2','ROB1','ROB2','LPC1','LPC2'};
nchan = length(chan_names);

% filename
fname = [direc,dfile];

tic
[~, sf, Tdata, Vdata, mtrig, dtrig] = LoadMCD(fname); %opens data file
toc

save([direc,date,'_',drugCondition,'_data_',recPhase,'.mat'], 'Vdata','dtrig');

% specify channel for button press
all_trig1 = dtrig == 3;    %change this according to digital channel
% all_trig2 = dtrig == 2;
all_trig1_smooth = smoothdtrig(all_trig1);
% all_trig2_smooth = smoothdtrig(all_trig2);

% assign odor variables to odor name on spreadsheet
odor1 = Hex;
odor2 = EMB;
odor3 = Ger;
odor4 = Lim;
odor5 = Cit;

odorList = {'Hex','EMB','Ger','Lim','Cit'};

% design filter used for movement artefact removal
% lpFilt = designfilt('lowpassiir','FilterOrder',5, ...
%          'PassbandFrequency',5,'PassbandRipple',0.2, ...
%          'SampleRate',sf);
lpFilt = designfilt('lowpassfir','PassbandFrequency',2, ...
         'StopbandFrequency',6.5,'PassbandRipple',0.5, ...
         'StopbandAttenuation',65,'DesignMethod','kaiserwin','SampleRate',sf);
% fvtool(lpFilt)

% Note: I find that fir produces less edge effects than iir

maf = [1,1,1,1,1,1]; % movement artefact flag. If 1 then correct for 
                     % movement artefact by subtracting lowpass filtered signal

% usechan = [1,2,3,4,5,6]; % NEW SYSTEM
usechan = [8,7,6,5,4,3]; % channels to use for analysis OLD SYSTEM

%% Extract presentation windows

pre = 2; % time interval before presentation onset (s)
post = 5; % time interval after presentation onset (s)

skipdp = ceil(7*sf); % # of data points to skip durring dtrig scan

% get timepts of each odor onset
all_inds1 = scandtrig_while(all_trig1,skipdp);
%     all_inds2 = scandtrig_while(all_trig2,skipdp);
%     all_inds = sort([all_inds1;all_inds2]);
odor1_inds = all_inds1(odor1);
n1 = length(odor1_inds);
odor2_inds = all_inds1(odor2);
n2 = length(odor2_inds);
odor3_inds = all_inds1(odor3);
n3 = length(odor3_inds);
odor4_inds = all_inds1(odor4);
n4 = length(odor4_inds);
odor5_inds = all_inds1(odor5);
n5 = length(odor5_inds);

%Uncomment to check indent times
figure
plot(all_trig1);
hold on
plot((all_inds1),.5,'o');

%make matrices of all trials per odor
if sum(rm_offtrials) > 0
    loff_inds(rm_offtrials) = []; % remove selectd off trials
end
% NOTE trialalign requires input data vector to be dim ndp x 1
for ch = 1:length(usechan)
    odor1trials.(chan_names{ch}) = trialalign(Vdata(usechan(ch),:),odor1_inds,pre,post,sf);
end

if sum(rm_ontrials) > 0
    lon_inds(rm_ontrials) = []; % remove selected on trials
end
for ch = 1:length(usechan)
    odor2trials.(chan_names{ch}) = trialalign(Vdata(usechan(ch),:),odor2_inds,pre,post,sf);
end

if sum(rm_ontrials) > 0
    lon_inds(rm_ontrials) = []; % remove selected on trials
end
for ch = 1:length(usechan)
    odor3trials.(chan_names{ch}) = trialalign(Vdata(usechan(ch),:),odor3_inds,pre,post,sf);
end

if sum(rm_ontrials) > 0
    lon_inds(rm_ontrials) = []; % remove selected on trials
end
for ch = 1:length(usechan)
    odor4trials.(chan_names{ch}) = trialalign(Vdata(usechan(ch),:),odor4_inds,pre,post,sf);
end

 if sum(rm_ontrials) > 0
    lon_inds(rm_ontrials) = []; % remove selected on trials
end
for ch = 1:length(usechan)
    odor5trials.(chan_names{ch}) = trialalign(Vdata(usechan(ch),:),odor5_inds,pre,post,sf);
end

%% detrend trials and remove movement artefact
twin = -pre:1/sf:post; % time window (s)

for ch = 1:length(usechan)
     for tt = 1:n1
         odor1trials.(chan_names{ch})(:,tt) = detrend(odor1trials.(chan_names{ch})(:,tt));
     end
    for tt = 1:n2
        odor2trials.(chan_names{ch})(:,tt) = detrend(odor2trials.(chan_names{ch})(:,tt));
    end
    for tt = 1:n3
        odor3trials.(chan_names{ch})(:,tt) = detrend(odor3trials.(chan_names{ch})(:,tt));
    end
    for tt = 1:n4
        odor4trials.(chan_names{ch})(:,tt) = detrend(odor4trials.(chan_names{ch})(:,tt));
    end

     for tt = 1:n5
        odor5trials.(chan_names{ch})(:,tt) = detrend(odor5trials.(chan_names{ch})(:,tt));
    end

    if maf(ch) == 1
         for tt = 1:n1
             filt_temp = filtfilt(lpFilt,odor1trials.(chan_names{ch})(:,tt));
             odor1trials.(chan_names{ch})(:,tt) = odor1trials.(chan_names{ch})(:,tt) - filt_temp;
         end
        for tt = 1:n2
            filt_temp = filtfilt(lpFilt,odor2trials.(chan_names{ch})(:,tt));
            odor2trials.(chan_names{ch})(:,tt) = odor2trials.(chan_names{ch})(:,tt) - filt_temp;
        end
         for tt = 1:n3
            filt_temp = filtfilt(lpFilt,odor3trials.(chan_names{ch})(:,tt));
            odor3trials.(chan_names{ch})(:,tt) = odor3trials.(chan_names{ch})(:,tt) - filt_temp;
         end
         for tt = 1:n4
            filt_temp = filtfilt(lpFilt,odor4trials.(chan_names{ch})(:,tt));
            odor4trials.(chan_names{ch})(:,tt) = odor4trials.(chan_names{ch})(:,tt) - filt_temp;
         end
         for tt = 1:n5
            filt_temp = filtfilt(lpFilt,odor5trials.(chan_names{ch})(:,tt));
            odor5trials.(chan_names{ch})(:,tt) = odor5trials.(chan_names{ch})(:,tt) - filt_temp;
         end
    end
end

clear tt ch

%% Set parameters for plotting
sff = 1; % save figure flag. If 1 then save figure. 2=early trials removed

lfs = 10; % legend fontsize
afs = 20; % axis fontsize

remEarlyTrials = 0; %1=true 0=false. flag to exclude the first 3 trials

chanPairs = [1 2; 3 4; 5 6];

%% Plot all trials
% plot channel pairs on same axes for all channels in usechan
if ~exist([direc,'traces'],'dir')
    mkdir(direc, 'traces');
end

for electrode=1:size(chanPairs,1) %for each electrode pair
    First = chanPairs(electrode,1);
    Second = chanPairs(electrode,2);

    % Odor sessions 1
    figure
    set(gcf,'position',[100,200,1700,700])
    for ii = 1:n1
        subplot(n1,1,ii)
        hold on
        plot(twin,odor1trials.(chan_names{First})(:,ii),'k')
        plot(twin,odor1trials.(chan_names{Second})(:,ii),'r')
        line([0,0],[-1500,1500],'Color','g','LineWidth',3);
        axis off
        axis tight
    end
    axis on
    legend({chan_names{First},chan_names{Second}},'fontsize',lfs,'position',[0.84,0.2,0.01,0.06])
    set(gca,'YTick',[],'fontsize',afs); set(gca,'XTick',[],'fontsize',afs)
    xlabel([num2str(pre+post),'s'])
    xlim([-pre post])
    suptitle(['Hexanone ',recPhase])
    if sff == 1
        fname = [direc,'traces\',chan_names{First},'&',chan_names{Second},'_Hexanone_',recPhase,'.png'];
        saveas(gcf,fname)
    end


    % Odor sessions 2
    figure
    set(gcf,'position',[100,200,1700,700])
    for ii = 1:n2
        subplot(n2,1,ii)
        hold on
        plot(twin,odor2trials.(chan_names{First})(:,ii),'k')
        plot(twin,odor2trials.(chan_names{Second})(:,ii),'r')
        line([0,0],[-1500,1500],'Color','g','LineWidth',3);
        axis off
        axis tight
    end
    axis on
    legend({chan_names{First},chan_names{Second}},'fontsize',lfs,'position',[0.84,0.2,0.01,0.06])
    set(gca,'YTick',[],'fontsize',afs);set(gca,'XTick',[],'fontsize',afs)
    xlabel([num2str(pre+post),'s'])
    xlim([-pre post])
    suptitle(['EMB ',recPhase])
    if sff == 1
        fname = [direc,'traces\',chan_names{First},'&',chan_names{Second},'_EMB_',recPhase,'.png'];
        saveas(gcf,fname)
    end

     % Odor sessions 3
    figure
    set(gcf,'position',[100,200,1700,700])
    for ii = 1:n3
        subplot(n3,1,ii)
        hold on
        plot(twin,odor3trials.(chan_names{First})(:,ii),'k')
        plot(twin,odor3trials.(chan_names{Second})(:,ii),'r')
        line([0,0],[-1500,1500],'Color','g','LineWidth',3);
        axis off
        axis tight
    end
    axis on
    legend({chan_names{First},chan_names{Second}},'fontsize',lfs,'position',[0.84,0.2,0.01,0.06])
    set(gca,'YTick',[],'fontsize',afs);set(gca,'XTick',[],'fontsize',afs)
    xlabel([num2str(pre+post),'s'])
    xlim([-pre post])
    suptitle(['Geraniol ',recPhase])
    if sff == 1
        fname = [direc,'traces\',chan_names{First},'&',chan_names{Second},'_Geraniol_',recPhase,'.png'];
        saveas(gcf,fname)
    end

     % Odor sessions 4
    figure
    set(gcf,'position',[100,200,1700,700])
    for ii = 1:n4
        subplot(n4,1,ii)
        hold on
        plot(twin,odor4trials.(chan_names{First})(:,ii),'k')
        plot(twin,odor4trials.(chan_names{Second})(:,ii),'r')
        line([0,0],[-1500,1500],'Color','g','LineWidth',3);
        axis off
        axis tight
    end
    axis on
    legend({chan_names{First},chan_names{Second}},'fontsize',lfs,'position',[0.84,0.2,0.01,0.06])
    set(gca,'YTick',[],'fontsize',afs);set(gca,'XTick',[],'fontsize',afs)
    xlabel([num2str(pre+post),'s'])
    xlim([-pre post])
    suptitle(['Limonene ',recPhase])
    if sff == 1
        fname = [direc,'traces\',chan_names{First},'&',chan_names{Second},'_Limonene_',recPhase,'.png'];
        saveas(gcf,fname)
    end

%      Odor sessions 5
    figure
    set(gcf,'position',[100,200,1700,700])
    for ii = 1:n5
        subplot(n5,1,ii)
        hold on
        plot(twin,odor5trials.(chan_names{First})(:,ii),'k')
        plot(twin,odor5trials.(chan_names{Second})(:,ii),'r')
        line([0,0],[-1500,1500],'Color','g','LineWidth',3);
        axis off
        axis tight
    end
    axis on
    legend({chan_names{First},chan_names{Second}},'fontsize',lfs,'position',[0.84,0.2,0.01,0.06])
    set(gca,'YTick',[],'fontsize',afs);set(gca,'XTick',[],'fontsize',afs)
    xlabel([num2str(pre+post),'s'])
    xlim([-pre post])
    suptitle(['Citronellol ',recPhase])
    if sff == 1
        fname = [direc,'traces\',chan_names{First},'&',chan_names{Second},'_Citronellol_',recPhase,'.png'];
        saveas(gcf,fname)
    end
end

clear First Second
    
%% Power plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YLp = [0 12]; % Ylim for power plots
ylabp = 'log(Power)'; % ylabel for power plots
YLc = [0 3]; % Ylim for coh plots
XL = [0 125]; % Xlim (frequency Hz)

% define parameters
% params.tapers = [thbw_prod (2*thbw_prod - 1)]; % [thbw ntapers]
params.tapers = [6 11]; % [thbw ntapers] %from LK [4 7]
params.pad = 0; % zero padding
params.Fs = sf; % sampling fq
params.fpass = [1 150]; % fq range
params.err = [2 1]; % error codes
params.trialave = 0; % average over trials? yes(1) no(0)

movingwin = [post+pre post+pre]; % [window length, step size] in (s)

ROI = 2*sf:5*sf; %Was 1.2 and 2.75

% exclude early trials from power plots
% ** no longer necessary because we can do this easier later on!!! 
if remEarlyTrials==1 && strcmp(recPhase,'Pre')==1
    trialsToUse=4:10;
elseif remEarlyTrials==1 && strcmp(recPhase,'Pre')==0
    trialsToUse=[4:10 14:20];
elseif remEarlyTrials==0 && strcmp(recPhase,'Pre')==1
    trialsToUse=1:10;
elseif remEarlyTrials==0 && strcmp(recPhase,'Pre')==0
    trialsToUse=1:20;
end

% create temporary processing matrix from master data matrix
% dmat dimensions are samples x trials
for iCh = 1:length(usechan)
    chName = chan_names{iCh};
    dmat.(odorList{1}).(chName) = odor1trials.(chName)(ROI,trialsToUse);
    dmat.(odorList{2}).(chName) = odor2trials.(chName)(ROI,trialsToUse);
    dmat.(odorList{3}).(chName) = odor3trials.(chName)(ROI,trialsToUse);
    dmat.(odorList{4}).(chName) = odor4trials.(chName)(ROI,trialsToUse);
    dmat.(odorList{5}).(chName) = odor5trials.(chName)(ROI,trialsToUse);
end
   
%     Spacer for subtraction if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcuate coherence and PSD (power spectral density)
for iOd=1:5
    odName = odorList{iOd};
    [C.(odName).('LOB1_LPC1'),ph,s12,S.(odName).('LOB1'),S.(odName).('LPC1'),f,confc,phie,cerr]=coherencyc(dmat.(odName).('LOB1'),dmat.(odName).('LPC1'),params);
    [C.(odName).('LOB2_LPC2'),~ ,~ , S.(odName).('LOB2'),S.(odName).('LPC2'),~,~    ,~   ,~   ]=coherencyc(dmat.(odName).('LOB2'),dmat.(odName).('LPC2'),params);
    [C.(odName).('LOB1_LOB2'),~ ,~ ,          ~         ,         ~         ,~,~    ,~   ,~   ]=coherencyc(dmat.(odName).('LOB1'),dmat.(odName).('LOB2'),params);
    [C.(odName).('LPC1_LPC2'),~ ,~ ,          ~         ,         ~         ,~,~    ,~   ,~   ]=coherencyc(dmat.(odName).('LPC1'),dmat.(odName).('LPC2'),params);
    [C.(odName).('ROB1_ROB2'),~ ,~ , S.(odName).('ROB1'),S.(odName).('ROB2'),~,~    ,~   ,~   ]=coherencyc(dmat.(odName).('ROB1'),dmat.(odName).('ROB2'),params);    
    [C.(odName).('LOB1_ROB1'),~ ,~ ,          ~         ,          ~        ,~,~    ,~   ,~   ]=coherencyc(dmat.(odName).('LOB1'),dmat.(odName).('ROB1'),params);
    [C.(odName).('LOB2_ROB2'),~ ,~ ,          ~         ,          ~        ,~,~    ,~   ,~   ]=coherencyc(dmat.(odName).('LOB2'),dmat.(odName).('ROB2'),params);
end

% take log of power spectra
for iOd=1:length(odorList)
    odName = odorList{iOd};
    for iChan=1:length(usechan)
        chanName = chan_names{iChan};
        logpwr.(odName).(chanName) = log(S.(odName).(chanName));
    end
end
 
% save workspace data
if sff == 1
    save([direc,'workspace_',recPhase,'.mat'], ...
        'logpwr','remEarlyTrials','C','sf','post','pre',...
        'f','YLp','chan_names','ylabp','XL','all_inds1','Rname',...
        'Nsess','odorList','recPhase','direc','usechan',...
        'odor1trials','odor2trials','odor3trials','odor4trials','odor5trials');
end
    
%% plotting power figures

% make folders to save figures
if sff==1 && ~exist([direc,'individual power plots'],'dir')
    mkdir(direc, 'individual power plots');
elseif sff==2 && ~exist([direc,'individual power plots (RET)'],'dir')
    mkdir(direc, 'individual power plots (RET)');
end

% plot power spectra for each set of electrodes, for each odor. 
for iOd=1:5
    odName=odorList{iOd};
    for iChanPair=1:2:6
        chanPair = chan_names(iChanPair:iChanPair+1);
        figure;
    %     subplot(2,1,1)
        hold on
        H1=shadedErrorBar(f,mean(logpwr.(odName).(chanPair{1}),2),std(logpwr.(odName).(chanPair{1}),0,2),'k');
        H2=shadedErrorBar(f,mean(logpwr.(odName).(chanPair{2}),2),std(logpwr.(odName).(chanPair{1}),0,2),'r',0.1);
        hold off
        ylim(YLp);
        legend([H1.mainLine H2.mainLine],chanPair,'location','northeast')
        set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
        title([odorList{iOd},' ',recPhase]);
        ylabel(ylabp)
        xlim(XL)
        xlabel('Frequency (Hz)')
        if sff == 1
            fname = [direc,'individual power plots\',chanPair{1},'&',chanPair{2},'_' odorList{iOd} '_power_',recPhase,'.png'];
            saveas(gcf,fname)
        elseif sff == 2
            fname = [direc,'individual power plots (RET)\',chanPair{1},'&',chanPair{2},'_' odorList{iOd} '_power_',recPhase,'.png'];
            saveas(gcf,fname)
        end
    end
end
    
%%    %%%%%Plot coherence values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist([direc,'coherence plots'],'dir')
    mkdir(direc, 'coherence plots');
end

% plot coherence values LOB-LPC
for iOd=1:5
    odName=odorList{iOd};
    figure
    a = C.(odName).('LOB1_LPC1'); 
    b = C.(odName).('LOB2_LPC2'); 
    c = C.(odName).('LOB1_LOB2'); 
    d = C.(odName).('LPC1_LPC2');
    
    subplot(2,1,1)
    hold on
    H1=shadedErrorBar(f,mean(a,2),std(a,0,2),'r');
    H2=shadedErrorBar(f,mean(b,2),std(b,0,2),'b',0.1);
    legend([H1.mainLine,H2.mainLine],'LOB1-canPC1','LOB2-canPC2','location','northeast')
    set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
    title(odName)
    ylabel('Coherence')
    xlabel('Hz')
    ylim([0 1.2])
    
    subplot(2,1,2)
    hold on
    H1=shadedErrorBar(f,mean(c,2),std(c,0,2),'r');
    H2=shadedErrorBar(f,mean(d,2),std(d,0,2),'b',0.1);
    legend([H1.mainLine,H2.mainLine],'LOB1-LOB2','canPC1-canPC2','location','northeast')
    set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
    title(odName)
    ylabel('Coherence')
    xlabel('Hz')
    ylim([0 1.2])
    if sff == 1
        fname = [direc,'coherence plots\',odName,'_coherence_LOB-LPC_',recPhase];
        saveas(gcf,fname)
    end
end

%% plot coherence values LOB-ROB
for iOd=1:5
    odName=odorList{iOd};
    figure
    a = C.(odName).('LOB1_ROB1'); 
    b = C.(odName).('LOB2_ROB2'); 
    c = C.(odName).('LOB1_LOB2'); 
    d = C.(odName).('ROB1_ROB2');
    
    subplot(2,1,1)
    hold on
    H1=shadedErrorBar(f,mean(a,2),std(a,0,2),'r');
    H2=shadedErrorBar(f,mean(b,2),std(b,0,2),'b',0.1);
    legend([H1.mainLine,H2.mainLine],'LOB1-ROB1','LOB2-ROB2','location','northeast')
    set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
    title(odName)
    ylabel('Coherence')
    xlabel('Hz')
    ylim([0 1.2])
    
    subplot(2,1,2)
    hold on
    H1=shadedErrorBar(f,mean(c,2),std(c,0,2),'r');
    H2=shadedErrorBar(f,mean(d,2),std(d,0,2),'b',0.1);
    legend([H1.mainLine,H2.mainLine],'LOB1-LOB2','ROB1-ROB2','location','northeast')
    set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
    title(odName)
    ylabel('Coherence')
    xlabel('Hz')
    ylim([0 1.2])
    if sff == 1
        fname = [direc,'coherence plots\',odName,'_coherence_LOB-ROB_',recPhase];
        saveas(gcf,fname)
    end
end

%% %%%%Spectrograms%%%%%

params.tapers = [2 3];
params.pad = 0;
params.Fs = sf;
params.fpass = [15 100];
params.err = [2 0.05];
params.trialave = 1;

movingwin = [0.4 0.1];

ROI = 2:6*sf; 

for iCh = 1:6
    chName = chan_names{iCh};
    dmatS.(odorList{1}).(chName) = odor1trials.(chName)(ROI,:);
    dmatS.(odorList{2}).(chName) = odor2trials.(chName)(ROI,:);
    dmatS.(odorList{3}).(chName) = odor3trials.(chName)(ROI,:);
    dmatS.(odorList{4}).(chName) = odor4trials.(chName)(ROI,:);
    dmatS.(odorList{5}).(chName) = odor5trials.(chName)(ROI,:);
end

for iOd=1:5
    odName = odorList{iOd};
    for iCh=1:6
        chName = chan_names{iCh};
        [Spect.(odName).(chName),tSpec,fSpec,SerrSpec]=mtspecgramc(dmatS.(odName).(chName),movingwin,params);
    end
end

%%%%%The following lines appear to be used to analyze "late" data

save([direc,'spectrogram_variables_',recPhase,'.mat'], ...
    'Spect','fSpec', 'tSpec')

%% plotting spectrograms
if ~exist([direc,'spectrograms'],'dir')
    mkdir(direc, 'spectrograms');
end

for iOd=1:5
    odName = odorList{iOd};
    figure;
    for iCh=1:length(usechan)
        subplot(3,2,iCh);
        chName = chan_names{iCh};
        pcolor(tSpec,fSpec,log(Spect.(odName).(chName))');
        shading interp; 
        %ylim([30 90]);
        hold on;
        line([2,2],[0,100],'Color','w','LineWidth',2); %Ensure Pre window%%%
        title([odName,' ',chName,' ',recPhase])
        xlabel('Time (s)')
        ylabel('Frequency (Hz)');
        set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
    end
    if sff == 1
        fname = [direc,'spectrograms\',odName '_spectrograms_',recPhase,'.png'];
        saveas(gcf,fname)
    end
end


