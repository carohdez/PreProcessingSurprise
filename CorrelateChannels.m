function [] = CorrelateChannels(channel, session, recording)
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
cd /mnt/homes/home024/pmurphy/meg_data/surprise/  % specify path where dataset is located
savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
dir_behav = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Data/';
dir_ET = [dir_behav,'QNV','/S',num2str(session),'/Eyetracking/'];
%load([savepath 'QNV_3_1_ascdat.mat']); % Eyelink data in ascii
%load([savepath 'qnvMEG_data.mat']); % preprocessed MEG file

dataset = ['QNV-3_Surprise_20170519_',recording,'.ds'];
if session == 1, dataset = ['QNV-1_Surprise_20170517_',recording,'.ds'];end
if session == 2, dataset = ['QNV-2_Surprise_20170518_',recording,'.ds'];end
if session == 4, dataset = ['QNV-4_Surprise_20170523_',recording,'.ds'];end

MEG_fileinfo = ft_read_event(dataset);
cfg = [];
cfg.trl = [1 MEG_fileinfo(end).sample, 0];
cfg.channel = 'EEG057';
cfg.datafile = [dataset];
cfg.continuous = 'yes';
%cfg.hpfilter = 'yes';
cfg.bpfilter = 'yes'; % Apply band pass filter to channel
%cfg.hpfreq = 1;
cfg.bpfreq = [0.029 400];
MEG_data = ft_preprocessing(cfg);   

channel = 1; % because once we preprocessed just the EEG057, we get just that channel

% For the regular sessions, we obtain the triggers from MEG UPPT001 channel
triggersMEG = [];
for e = 1:length(MEG_fileinfo)
    if strcmp(MEG_fileinfo(e).type,'UPPT001')
        triggersMEG(end+1,:) = [MEG_fileinfo(e).value MEG_fileinfo(e).sample];        
    end
end

triggersET = [];
all_trMEG = [];
all_trET = [];

% subplot(211)
% plot(MEG_data.trial{1,1}(channel,:))
ETfile_ini = 1;
ETfile_end = 5;
if(session==1) && strcmp(recording,'02')
    ETfile_ini = 6;
    ETfile_end = 8;
end
if(session==2 || session ==3) && strcmp(recording,'02')
    ETfile_ini = 6;
    ETfile_end = 9;
end
if(session==4) && strcmp(recording,'01')
    ETfile_end = 4;
end
if(session==4) && strcmp(recording,'02')
    ETfile_ini = 5;
    ETfile_end = 7;
end
    
for ETfile = ETfile_ini:ETfile_end
    %load([savepath 'QNV_' num2str(session) '_' num2str(ETfile) '_ascdat.mat']); % Eyelink data in ascii
    ascii_file = dir([dir_ET,'QNV_',num2str(session),'_',num2str(ETfile),'.asc']);
    ascdat = read_eyelink_ascNK_AU([dir_ET, ascii_file(1).name]);
    %ascdat = read_eyelink_ascNK_AU([dir_ET, 'QNV_',num2str(session),'_',num2str(ETfile),'.asc']);
    ET_start_smp = []; ET_end_smp = [];
    ET_event_val = cell(1,length(ascdat.msg));
    ET_event_valnum = zeros(1,length(ascdat.msg)); ET_event_smp = ET_event_valnum;
    for i=1:length(ascdat.msg)
        strtok=regexp(ascdat.msg{i},'(\t|\s)','split');
        if ~isempty(find(ascdat.dat(1,:) == str2double(strtok{2})));% si el id msg aparece en dat
            if ~strcmp('Start',strtok{3}) && ~strcmp('TRIALID',strtok{3})
                ET_event_smp(i) = find(ascdat.dat(1,:) == str2double(strtok{2})); % vector with indexes (column) of file ascdat correspondent to triggers (0 if there is not trigger value)
                ET_event_valnum(i) = str2double(strtok{3}); % vector with id trigger

            end
        end
    end

    % add extra channel with triggers to ascdat.dat to be able to upsample the trigger samples in the right way.
    ascdat.dat(5,:) = 0;
    ascdat.dat(5,ET_event_smp(ET_event_smp>0)) = ET_event_valnum(ET_event_smp>0);

    ET_start_smp = [ET_start_smp find(ascdat.dat(5,:) == 1)]; % only samples corresponding to trigger 1 (begining of the block)
    ET_end_smp = [ET_end_smp find(ascdat.dat(5,:) == 2)]; % only samples corresponding to trigger 2 (end of the block)  

    % Define trial for whole run ET           
    ET_data = [];
    ET_data.label ={'h'; 'v'; 'p';'trg'};%{'h'; 'v'; 'p'};
    ET_data.fsample = ascdat.fsample;
    ET_data.trial{1,1} = ascdat.dat(2:end,ET_start_smp(1):ET_end_smp(end));
    ET_data.time{1,1} = 0:1/ET_data.fsample:length(ascdat.dat(1,ET_start_smp(1):ET_end_smp(end)))/ET_data.fsample-1/ET_data.fsample;
    ET_data.sampleinfo(1,:) = [1 length(ascdat.dat(1,ET_start_smp(1):ET_end_smp(end)))];
    ET_data.cfg.trl(1,:) = [ET_start_smp(1), ET_end_smp(end), 0];

    %upsample
    cfg = [];
    %cfg.resamplefs = 1200;
    cfg.resamplefs = (765127/765172)*1200;
    cfg.detrend = 'no';
    ET_data2 = ft_resampledata(cfg, ET_data);

    % Correlation to find optimal match: 2-step approach
    % rough and broad search
    hcorr_sample = []; % MEG index sample max correlation
    hcoeff_corr = 0; % max coefficient correlation
    for i = 1:2 % 1=xpos, 2=ypos, 3=pupil
        % Apply band pass filter
        % ET_data2.trial{1}(i,:) = ft_preproc_bandpassfilter(ET_data2.trial{1}(i,:), 1200, [1 15], [], [], 'onepass');
        neg_search = 1;
        indx = [];
        %nsteps = length(MEG_data.trial{1})-length(ET_data2.trial{1})-neg_search;
        nsteps = length(MEG_data.trial{1})-length(ET_data2.trial{1});
        stepsize = 250;
        steps = 1:stepsize:nsteps;
        r = zeros(size(steps));
        for j = 1:length(r)
            %r(j) = corr(squeeze(MEG_data.trial{1,1}(channel,steps(j):(steps(j)+length(ET_data2.trial{1})-neg_search)))', squeeze(ET_data2.trial{1}(i,neg_search:end))');
            r(j) = corr(squeeze(MEG_data.trial{1,1}(channel,steps(j):(steps(j)+length(ET_data2.trial{1})-1)))', squeeze(ET_data2.trial{1}(i,:))');
        end
        [~, indx] = max(abs(r));
        bs_hcorr_sample = steps(indx);
        save([savepath,'qnv',num2str(session),'_EOGcoef_',num2str(ETfile),num2str(i),'_broad.mat'],'r','-v7.3');
        
        % slow and detailed search
        nsteps = stepsize;
        r = zeros(1, nsteps*2);
        steps = steps(indx)-nsteps:steps(indx)+nsteps;
        if ~isempty(steps)
            for j = 1:length(steps)-1
                %r(j) = corr(squeeze(MEG_data.trial{1,1}(channel,steps(j):(steps(j)+length(ET_data2.trial{1})-neg_search)))', squeeze(ET_data2.trial{1}(i,neg_search:end))');
                r(j) = corr(squeeze(MEG_data.trial{1,1}(channel,steps(j):(steps(j)+length(ET_data2.trial{1})-1)))', squeeze(ET_data2.trial{1}(i,:))');
            end
    %         subplot(212)
    %         plot(r)
    %         title(sprintf('Detailed search corr  %s file %d session %d',measure_name, ETfile, session ))   
            [coeff, indx] = max(abs(r));
            if abs(coeff) > hcoeff_corr
                hcoeff_corr = coeff;
                hcorr_sample = steps(indx);
    %            hcorr_sample = bs_hcorr_sample + (indx - stepsize);
            end;

        else
            fprintf('detailed research not possible, steps==0');
        end

        save([savepath,'qnv',num2str(session),'_EOGcoef_',num2str(ETfile),num2str(i),'_detailed.mat'],'r','-v7.3');

        %save plots
        %print([savepath,'qnv',num2str(session),'_EOGcoef_',num2str(ETfile),'_',measure_name,'_plot'],'-dpng');

    end
    % title(sprintf('Eyelink %s %d over MEG after correlation peak detection',measure_name, ETfile )) 
    % legend('eog','x pos','y pos','pupil')
    %print([savepath,'qnv_EOG_ET_',ETfile,'_plot'],'-dpng');

    % Throwing out all samples pre-block start to align with ET data segment used for initial correlation
    ascdat.dat = ascdat.dat(:,ET_start_smp(1):end);
    ET_event_smp = ET_event_smp-ET_event_smp(ET_event_valnum==1)+1;
    
    % Upsampling triggers
    ET_smp_upsampled = round(ET_event_smp(ET_event_smp>0)*((765127/765172)*1.2));
    ET_event_smp_tr = ET_event_smp(find(ET_event_smp>0));
    for i=1:length(ET_smp_upsampled)
        triggersET(end+1, :) = [ascdat.dat(5,ET_event_smp_tr(i)) hcorr_sample + ET_smp_upsampled(i)];
    end
end
% startblMEG = find(triggersMEG(:,1)==1);
% startblET = find(triggersET(:,1)==1);
% endblMEG = find(triggersMEG(:,1)==2);
% endblET = find(triggersET(:,1)==2);
% 
% figure
% subplot(511)
% % Obtain the triggers for each block to compare
% for i=1:length(startblMEG)
%     subplot(5,1,i)
%     ctrigsET = triggersET(startblET(i):endblET(i),:);   % pull only triggers for this block
%     ctrigsMEG = triggersMEG(startblMEG(i):endblMEG(i),:);   % pull only triggers for this block
%     plot(ctrigsET(ctrigsET(:,1)==11,2)-ctrigsMEG(ctrigsMEG(:,1)==11,2))
%     xlabel('# of trial')
%     ylabel('diff in samples')
% %     subplot(5,2,i)
% %     plot(MEG_data.trial{1,1}(channel,startblMEG(i):endblMEG(i)))
% %     xlabel('EOG')
%     ctrigsET(:,3)=i;
%     ctrigsMEG(:,3)=i;
%     all_trMEG = [all_trMEG; ctrigsMEG];
%     all_trET = [all_trET; ctrigsET];
% end
% print([savepath,'qnv_discrep_MEG_ET_',num2str(session),'_',recording,'_plot'],'-dpng');

% save triggers of the recording
save([savepath,'qnv_triggersMEG_',num2str(session),'_',recording,'.mat'],'triggersMEG','-v7.3');
save([savepath,'qnv_triggersET_',num2str(session),'_',recording,'.mat'],'triggersET','-v7.3');
end