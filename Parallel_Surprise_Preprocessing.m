% Top-level script for running multiple 2-accumulator LCA fits in parallel
% Uses qsubcellfun.m to submit jobs to TORQUE


% Path stuff
addpath(genpath('/mnt/homes/home024/chernandez'))
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221/qsub'

% Subject/model variant stuff
allsubj = {'DHB' 'EXF' 'TFD'};
sessions = {'1' '2' '3'}; 
recordings = {'01' '02'}; 


% Construct cell arrays for calling jobs
subj_in={}; sess_in={}; rec_in={};
for i = 1:length(allsubj);
    for j = 1:length(sessions);
        for k = 1:length(recordings);
            subj_in{end+1} = allsubj{i};
            sess_in{end+1} = sessions{j};
            rec_in{end+1} = recordings{k};
        end
    end
end

% Submit jobs to TORQUE
setenv('TORQUEHOME', 'yes')   % not sure what this does..........   
cd('~/Data/')

memreq = 4*1024^3;  % memory required (1GB = 1024^3)
timreq = 5*60*60;   % time required (hours*minutes*seconds)

% qsubcellfun(@MEGpreprocessing, subj_in, sess_in, rec_in, 'compile', 'no',...
%         'memreq', memreq, 'timreq', timreq, 'stack', 1, 'StopOnError', false, 'backend', 'torque');
    
subj_in={1}
qsubcellfun(@MEGPreprocessingSummary,subj_in, 'compile', 'no',...
        'memreq', memreq, 'timreq', timreq, 'stack', 1, 'StopOnError', false, 'backend', 'torque');

