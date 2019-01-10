%% 64Channel SiProbe Data Filtering %%%%%%%%%%%%%%%%%%%%%%%%
% The data from this probe (RL75) is noisy with 60hz background noise and
% other artefact. This program uses the power spectra of each trial to help
% sort good trials from bad. Human input is needed to: 
% 1) specify which trials have big, obvious artifact across all channels
%    (from unplugging, head movements, etc.) 
% 2) mark individual trials that are swamped by aggressive 60hz noise
%    (borderline cases only)
%
% The program automatically saves after each channel is combed for 60hz
% noise, and loads the saved file if it already exists in the folder. No
% need to sort through an entire recording in one sitting. 
%
% Author: Emily Tao, 2018

%% SECTION 1: load and initialize data

addpath('\\ssdfiles.uchicago.edu\kaylab\Code.Repository\MATLAB\chronux_2_11');

clearvars;

% open file
[file,path] = uigetfile('*.mat');
load(fullfile(path,file));
if exist(fullfile(path,'autosave.mat'),'file')
    load(fullfile(path,'autosave.mat'));
end

% decimate from 3k to 1.5k smpl
nTrials = size(allChannelData,2);
if size(allChannelData,1) > 4500
    allChannelDataDec(1:4500,1:nTrials,1:64) = nan;
    for i=1:64
        for j=1:nTrials
            allChannelDataDec(1:4500,j,i) = downsample(allChannelData(:,j,i),2);
        end
    end
    allChannelData = allChannelDataDec;
    goodTrialsTTL=goodTrials;
    save(fullfile(path,file),'allChannelData','goodTrialsTTL','-append');
end

% plotting variables consistent for all trials
triallength = size(allChannelData,1);
smpl = 1500; 
trialtime = -1:1/smpl:1.9995;
trialtime = trialtime';
trialtime = trialtime(1:triallength,1);

%% SECTION 2:
% look for bad trials wrt signal quality - eg cord disconnected, no signal
% if basic signal quality is bad, don't bother sorting for 60hz noise 
% because it will be all bad.

goodTrialsDt_ch = []; %matrix of good trials per channel
ch = str2num(input('Test channel:\n>>  ', 's')); %input any channel to test
ind=1;
while ~isempty(ch) %test as many channels as desired; quit w/empty input
    dmat = allChannelData(:,:,ch);
    goodTrialsDt_ch(ind,:) = ones(size(nspk)); %initialize all 'good' trials
    
    figure; %plotting one by one each trial that has TTL data
    for i=1:nTrials 
        if goodTrialsTTL(i)
            plot( trialtime, dmat(:,i), 'k' ); 
            xlim( [-1,2] );
            title(['Trial Number  ' num2str(i)]);
            xlabel('Time (s)');                
            saveTrial = input('Save Trial? 1=yes, 0=no\n>>  ', 's'); 
            pause(0.01);
            if isempty(saveTrial) 
                %do nothing; keep as is
            elseif strcmp(saveTrial,'0')
                goodTrialsDt_ch(ind,i) = 0; 
                %mark trial as 'bad'
            else
                %do nothing; keep as is
            end    
        end
    end
    ind = ind+1;
    ch = str2num(input('Test channel: 1-64 to test, ENTER to stop\n>>  ', 's'));
end

%if any trial tested has 'good' data on any tested channel, keep it
goodTrialsDt = logical(sum(goodTrialsDt_ch,1));
save(fullfile(path,file),'goodTrialsDt','-append');

%% SECTION 3: filter trials with 60 Hz noise
% computer sorts 'clear cut' cases based on power spectra. borderline cases
% shown to human eye for final decision

if ~exist('saveTrialsAll','var') 
    saveTrialsAll = repmat(goodTrialsDt',1,64);
end
ch = str2num(input('Start from channel:\n>>  ', 's'));

%set acceptable 'clear cut' cases
HIGHLIMIT = 2.4; % if diff in log pwr greater than this, auto-BAD
LOWLIMIT = 1.7; % if diff in log pwr lower than this, auto-OK
doNextChannel=true;
while doNextChannel 
    
    dmat = allChannelData(:,:,ch);  % LFP data matrix; be sure to put in the right channel
    goodTrialsNs = goodTrialsTTL .* goodTrialsDt;

    meanlog60pwr= zeros(size(nspk)); % mean power around 60hz
    meanlogContext= zeros(size(nspk)); % mean power of freq slightly above/below 60hz
    dPwr= zeros(size(nspk)); % difference in power of 60hz vs surrounding freq's

    % pwr spectra per trial
    params.tapers = [3 5];
    params.pad = 1;
    params.Fs = 1500;
    params.fpass = [3 110];
    params.err = [2 0.05];
    params.trialave = 0;

    for i=1:nTrials
        if goodTrialsNs(i)
            [S,f,Serr]=mtspectrumc(dmat(:,i),params);
            meanlog60pwr(i)=mean(log(S(f>59.5 & f<60.5)));
            meanlogContext(i) = mean(log(S((f>58 & f<58.5)|(f>61.5&f<62))));
            dPwr(i)= meanlog60pwr(i)-meanlogContext(i);
            maxPwr55 = max(log(S(f<55)));
            % if dPwr higher than HIGHLIMIT or mean power at 60hz is
            % greater than power at all lower frequencies, mark BAD
            if dPwr(i) > HIGHLIMIT || meanlog60pwr(i) > maxPwr55
                goodTrialsNs(i) = false; 
            end
        end
    end
    fprintf('Bad trials found: %i\n', length(find(dPwr>HIGHLIMIT)));
    fprintf('Good trials found: %i\n', length(find(dPwr<LOWLIMIT & dPwr>0)));
    fprintf('Remaining: %i\n', length(find(dPwr>LOWLIMIT & dPwr<=HIGHLIMIT)));

    % show border cases for human input
    for i=1:nTrials
        if dPwr(i)>LOWLIMIT && dPwr(i)<=HIGHLIMIT

            % plot the raw signal
            subplot(2,1,1);
            plot( trialtime, dmat(:,i), 'k' );
            xlim( [-1,2] );
            title(['Channel ' num2str(ch) '  Trial Number  ' num2str(i)]);
            xlabel('Time (s)');

            % plot the power spectra
            subplot(2,1,2);
            [S,f,Serr]=mtspectrumc(dmat(:,i),params);
            plot( f, log(S), 'k' );
            title(['dPwr = ' num2str(dPwr(i))]);

            % prompt user
            saveTrial = input('Save Trial? 1=yes, 0=no\n>>  ', 's'); 
            pause(0.01);
            if isempty(saveTrial) 
                %do nothing; keep as is
            elseif strcmp(saveTrial,'0')
                goodTrialsNs(i) = false; 
            else
                %do nothing; keep as is
            end    
        end
    end
    clf;
    saveTrialsAll(1:nTrials,ch) = logical(goodTrialsNs);
    
    % prompt user to continue
    doNextQuery = input('Do next channel? 1=yes, 0=no\n>>  ', 's');
    if ch==64
        disp('Channel 64 complete. You are finished!');
        doNextChannel = false;
    elseif strcmp(doNextQuery,'0')
        doNextChannel = false;
    else
        ch = ch+1;
    end
    
    % save trial decision to autosave file after each channel is completed
    % fileName= [recName '_DBT.mat'];
    save(fullfile(path,'autosave.mat'), 'saveTrialsAll', 'trialtime');
end

%% check the bad trials
% verify that the computer is sorting bad trials correctly
ch=41;
figure;
for i=1:nTrials
    if goodTrialsTTL(i) && ~saveTrialsAll(i,ch)
        dmat = allChannelData(:,:,ch);
        
        % plot the raw signal
        subplot(2,1,1);
        plot( trialtime, dmat(:,i), 'k' );
        xlim( [-1,2] );
        title(['Channel ' num2str(ch) '  Trial Number  ' num2str(i)]);
        xlabel('Time (s)');
        
        % plot the power spectra
        subplot(2,1,2);
        [S,f,Serr]=mtspectrumc(dmat(:,i),params);
        plot( f, log(S), 'k' );
        title(['dPwr = ' num2str(dPwr(i))]);
        
        % prompt user; user can recategorize if desired
        saveTrial = input('Save Trial? 1=no, 0=yes\n>>  ', 's'); 
        pause(0.01);
        if isempty(saveTrial) 
            %do nothing; keep as is
        elseif strcmp(saveTrial,'0')
            goodTrialsNs(i) = true; 
        else
            %do nothing; keep as is
        end    
    end
end

%% check the good trials
% verify that the computer is sorting good trials correctly
ch=41;
figure;
for i=1:nTrials
    if goodTrialsTTL(i) && saveTrialsAll(i,ch)
        dmat = allChannelData(:,:,ch);
        
        % plot the raw signal
        subplot(2,1,1);
        plot( trialtime, dmat(:,i), 'k' );
        xlim( [-1,2] );
        title(['Channel ' num2str(ch) '  Trial Number  ' num2str(i)]);
        xlabel('Time (s)');
        
        % plot the power spectra
        subplot(2,1,2);
        [S,f,Serr]=mtspectrumc(dmat(:,i),params);
        plot( f, log(S), 'k' );
        title(['dPwr = ' num2str(dPwr(i))]);
        
        % prompt user; user can recategorize if desired
        saveTrial = input('Save Trial? 1=no, 0=yes\n>>  ', 's'); 
        pause(0.01);
        if isempty(saveTrial) 
            %do nothing; keep as is
        elseif strcmp(saveTrial,'0')
            goodTrialsNs(i) = true; 
        else
            %do nothing; keep as is
        end    
    end
end