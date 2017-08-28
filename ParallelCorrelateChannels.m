% Top-level script for running multiple 2-accumulator LCA fits in parallel
% Uses qsubcellfun.m to submit jobs to TORQUE


% Path stuff
addpath(genpath('/mnt/homes/home024/chernandez'))
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221/qsub'

% Subject/model variant stuff
all_chans = {370,371,372,372,373,374,375};
all_chans = {362};
all_sessions = {1,2,4};
all_sessions = {3};
all_recordings = {'01','02'};
%all_recordings = {'01'};

% Construct cell arrays for calling jobs
chann_in={};
session_in={};
recording_in={};
for i = 1:length(all_chans);
    for j = 1:length(all_sessions);
        for k = 1:length(all_recordings)
            chann_in{end+1} = all_chans{i};
            session_in{end+1} = all_sessions{j};
            recording_in{end+1} = all_recordings{k};            
        end
    end
end

memreq = 4*1024^3;  % memory required (1GB = 1024^3)
timreq = 3*60*60;   % time required (hours*minutes*seconds)

qsubcellfun(@CorrelateChannels, chann_in, session_in, recording_in, 'compile', 'no',...
        'memreq', memreq, 'timreq', timreq, 'stack', 1, 'StopOnError', false, 'backend', 'torque');

