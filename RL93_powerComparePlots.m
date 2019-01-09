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

%load power data pre and post infusion
load('workspace_pre.mat','f','YLp','chan_names','ylabp','XL','all_inds1','Rname',...
            'Nsess','odorList','direc');
post_pwr_data = load('workspace_post.mat', 'logpwr');
pre_pwr_data = load('workspace_pre.mat', 'logpwr');

%make folder
if ~exist([direc,'power compare plots'],'dir')
        mkdir(direc, 'power compare plots');
end

% code for old format workspace
% post_pwr_data = load('workspace_post.mat', 'logpwr1', 'logpwr2', 'logpwr3', 'logpwr4');
% pre_pwr_data = load('workspace_pre.mat', 'logpwr1', 'logpwr2', 'logpwr3', 'logpwr4');


%% plot the pre- and post- infusion power plots on the same axis
for iOd=1:5
    odName=odorList{iOd};
    for iChanPair=1:2:6
        chanPair = chan_names(iChanPair:iChanPair+1);
        figure;
        hold on
        H3=shadedErrorBar(f,mean(pre_pwr_data.logpwr.(odName).(chanPair{1}),2),std(pre_pwr_data.logpwr.(odName).(chanPair{1}),0,2),'k',0.5);
        H4=shadedErrorBar(f,mean(pre_pwr_data.logpwr.(odName).(chanPair{2}),2),std(pre_pwr_data.logpwr.(odName).(chanPair{2}),0,2),'k',0.1);
        H1=shadedErrorBar(f,mean(post_pwr_data.logpwr.(odName).(chanPair{1}),2),std(post_pwr_data.logpwr.(odName).(chanPair{1}),0,2),'m',0.5);
        H2=shadedErrorBar(f,mean(post_pwr_data.logpwr.(odName).(chanPair{2}),2),std(post_pwr_data.logpwr.(odName).(chanPair{2}),0,2),'m',0.1);
        hold off
        ylim(YLp);
        legend([H1.mainLine H2.mainLine ...
            H3.mainLine H4.mainLine],...
            [chanPair{1},' post'],[chanPair{2},' post'],...
            [chanPair{1},' pre'],[chanPair{2},' pre'],...
            'location','northeast')
        %legend for single electrode comparison plots
%         legend([H2.mainLine H4.mainLine],...
%             [chanPair{2},' post'],[chanPair{2},' pre'],...
%             'location','northeast')
        set(gca,'box', 'off','FontWeight','bold', 'Fontsize', 14, 'Linewidth', 1.5)
        title([Nsess(1:6) ' ' Nsess(8:end) ': ' odorList{iOd} ' post vs. pre']);
        ylabel(ylabp)
        xlim(XL)
        if sff == 1
            fname = [direc,'power compare plots\',chanPair{1},'&',chanPair{2},'_',odorList{iOd},'_power_postvpre.png'];
            saveas(gcf,fname)
        end
    end
end
