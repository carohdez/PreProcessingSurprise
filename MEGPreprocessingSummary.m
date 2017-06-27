function [] = MEGPreprocessingSummary(code)

dir_prepr = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
subjfile = dir([dir_prepr,'*_Surprise_Prepr*','.mat']);
prepr_summary = [];
name = {};
tr_rejected = [];
for i=1:length(subjfile),
    load([dir_prepr, subjfile(i).name]);
    name = subjfile(i).name;
    remaining_tr= data.remaining_tr;
    nr_trials_beg=remaining_tr(1);
    %nr_tr_rejected = remaining_tr-nr_trials_beg;
    subject_info = {name remaining_tr data.trialinfo};
    prepr_summary = [prepr_summary;subject_info];
end
%save data file
savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';
save([savepath,'summary_preprocessed.mat'],'prepr_summary','-v7.3');
    