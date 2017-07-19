function [] = MEGPreprocessingSummary(code)

dir_prepr = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
subjfile = dir([dir_prepr,'*_Surprise_Prepr*','.mat']);
%subjfile = dir([dir_prepr,'PDP-2_Surprise_Prepr*','.mat']);
prepr_summary = [];
name = {};
tr_rejected = [];
for i=1:length(subjfile),
    load([dir_prepr, subjfile(i).name],'remaining_tr');
    load([dir_prepr, subjfile(i).name],'trl');
    name = subjfile(i).name;
    %remaining_tr= data.remaining_tr;
    nr_trials_beg=remaining_tr(1);
    %nr_tr_rejected = remaining_tr-nr_trials_beg;
    subject_info = {name remaining_tr any(trl(:,14)~=trl(:,15))};
    prepr_summary = [prepr_summary;subject_info];
end
%save data file
savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
save([savepath,'summary_preprocessed.mat'],'prepr_summary','-v7.3');
%
%load('/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/GSB-1_Surprise_Preprocessed_02.mat','trl')



%%% comparacion de metodos saccades %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_prepr = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
subjfile = dir([dir_prepr,'*_Surprise_Preprocessed_0*','.mat']);
prepr_summary = [];
name = {};
tr_rejected = [];
tr_end_sacc = [];
tr_end_diff = [];
subject = {};
session = {};
recording = {};
for i=1:length(subjfile),
    load([dir_prepr, subjfile(i).name],'remaining_tr');
    load([dir_prepr, subjfile(i).name],'trl');
    name = subjfile(i).name;
    subject = name(1:3);
    session =  name(5);
    recording = name(29:30);
    nr_trials_beg = remaining_tr(1);
    nr_trials_end = remaining_tr(length(remaining_tr));
    remaining_tr_no_sacc = remaining_tr;
    try
        load([dir_prepr, subject,'-',session,'_Surprise_Preprocessed_sacc_',recording,'.mat'],'remaining_tr');
        nr_trials_end_sac = remaining_tr(length(remaining_tr));
    catch
        nr_trials_end_sac = -1;
        remaining_tr = [];
    end
    
    
    %subject_info = {subject session recording nr_trials_beg nr_trials_end nr_trials_end_sac nr_trials_end-nr_trials_end_sac remaining_tr_no_sacc remaining_tr};
    subject_info = {subject session recording nr_trials_end-nr_trials_end_sac remaining_tr_no_sacc remaining_tr};
    prepr_summary = [prepr_summary;subject_info];
end
%save data file
savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
save([savepath,'summary_preprocessed_sacc.mat'],'prepr_summary','-v7.3');


for i=1:6
    [x, y, z] = eye_voltage2gaze(data.trial{tr_sacc(i)}, ranges, screen_x, screen_y, ch_mapping); figure, hold on, plot(x,'b'), plot(y,'r')
end
    
%[x, y, z] = eye_voltage2gaze(data.trial{tr_sacc(7)}, ranges, screen_x, screen_y, ch_mapping); figure, hold on, plot(x,'b'), plot(y,'r')

