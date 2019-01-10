%% Odor presentation analysis
% plot power spectrum data from different recording days to compare across
% time / drug condition

clearvars;

%initialize channels
ch1 = 1;
ch2 = 2;
ch3 = 3;
ch4 = 4;
ch5 = 5;
ch6 = 6;

%save figure flag 1=save 0=no
sff=1;

%names of folders of recs to compare
rec1fold='111418_LowAg';
rec2fold='111218_Saline';
recPhase='Pre';

%load power data from both folders
load([rec1fold '/workspace_' recPhase '.mat'],...
    'f','YLp','chan_names','ylabp','XL','all_inds1','Rname',...
    'Nsess','odorList','usechan');
rec1_pwr_data = load([rec1fold '/workspace_' recPhase '.mat'], 'logpwr','direc');
rec2_pwr_data = load([rec2fold '/workspace_' recPhase '.mat'], 'logpwr','direc');

%save plots to rec1 folder labeled with both rec names
savePath = ['power compare plots ' rec1fold ' ' rec2fold];
if ~exist([rec1_pwr_data.direc,savePath],'dir')
    mkdir(rec1_pwr_data.direc,savePath);
    mkdir([rec1_pwr_data.direc,savePath],'fig');
    mkdir([rec1_pwr_data.direc,savePath],'png');
end

%% plot the rec1 and rec2 power plots on the same axis 
for iOd=1:5
    odName=odorList{iOd};
    for iChanPair=1:2:length(usechan)
        chanPair = chan_names(iChanPair:iChanPair+1);
        figure;
        hold on
        H3=shadedErrorBar(f,mean(rec2_pwr_data.logpwr.(odName).(chanPair{1}),2),std(rec2_pwr_data.logpwr.(odName).(chanPair{1}),0,2),'k',0.5);
        H4=shadedErrorBar(f,mean(rec2_pwr_data.logpwr.(odName).(chanPair{2}),2),std(rec2_pwr_data.logpwr.(odName).(chanPair{2}),0,2),'k',0.1);
        H1=shadedErrorBar(f,mean(rec1_pwr_data.logpwr.(odName).(chanPair{1}),2),std(rec1_pwr_data.logpwr.(odName).(chanPair{1}),0,2),'m',0.5);
        H2=shadedErrorBar(f,mean(rec1_pwr_data.logpwr.(odName).(chanPair{2}),2),std(rec1_pwr_data.logpwr.(odName).(chanPair{2}),0,2),'m',0.1);
        hold off
        ylim(YLp);
        legend([H1.mainLine H2.mainLine ...
            H3.mainLine H4.mainLine],...
            [chanPair{1},' ',rec1fold(8:end)],[chanPair{2},' ',rec1fold(8:end)],...
            [chanPair{1},' ',rec2fold(8:end)],[chanPair{2},' ',rec2fold(8:end)],...
            'location','northeast')
        %legend for single electrode comparison plots
%         legend([H2.mainLine H4.mainLine],...
%             [chanPair{2},' post'],[chanPair{2},' pre'],...
%             'location','northeast')
        set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
        title([rec1fold(8:end) ' vs ' rec2fold(8:end) ': ' odorList{iOd} ' ' recPhase]);
        ylabel(ylabp)
        xlim(XL)
        if sff == 1
            fname = [rec1_pwr_data.direc,'power compare plots ',rec1fold,' ',rec2fold,'/fig/',...
                chanPair{1},'&',chanPair{2},'_',odorList{iOd},'_power_' recPhase '.fig'];
            saveas(gcf,fname)
            fname = [rec1_pwr_data.direc,'power compare plots ',rec1fold,' ',rec2fold,'/png/',...
                chanPair{1},'&',chanPair{2},'_',odorList{iOd},'_power_' recPhase '.png'];
            saveas(gcf,fname)
        end
    end
end
