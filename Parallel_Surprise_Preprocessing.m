% Top-level script for running multiple 2-accumulator LCA fits in parallel
% Uses qsubcellfun.m to submit jobs to TORQUE


% Path stuff
addpath(genpath('/mnt/homes/home024/chernandez'))
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221/qsub'

% Subject/model variant stuff
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ'};
%allsubj = {'EMB' 'HBC' 'JTB' 'OMF' 'PDP' 'QNV' 'TNB' 'TSJ'};%Subjects with exceptions
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

% QNV
subj_in{end+1} = 'QNV';
sess_in{end+1} = '4';
rec_in{end+1} = '01';

subj_in{end+1} = 'QNV';
sess_in{end+1} = '4';
rec_in{end+1} = '02';

%Subjects 3 reccordings
% EMB, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-8
subj_in{end+1} = 'EMB';
sess_in{end+1} = '2';
rec_in{end+1} = '03';

% JTB, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
subj_in{end+1} = 'JTB';
sess_in{end+1} = '2';
rec_in{end+1} = '03';

% PDP, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
subj_in{end+1} = 'PDP';
sess_in{end+1} = '2';
rec_in{end+1} = '03';

% TNB, Session 3
% rec 01: block 1, rec 02: blocks 2-5, rec 03: blocks 6-9
subj_in{end+1} = 'TNB';
sess_in{end+1} = '3';
rec_in{end+1} = '03';


% TSJ, Session 2
% rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
subj_in{end+1} = 'TSJ';
sess_in{end+1} = '4';
rec_in{end+1} = '03';

% Submit jobs to TORQUE
setenv('TORQUEHOME', 'yes')   % not sure what this does..........   
cd('~/Data/')

memreq = 4*1024^3;  % memory required (1GB = 1024^3)
timreq = 2*60*60;   % time required (hours*minutes*seconds)

qsubcellfun(@MEGpreprocessing, subj_in, sess_in, rec_in, 'compile', 'no',...
        'memreq', memreq, 'timreq', timreq, 'stack', 1, 'StopOnError', false, 'backend', 'torque');

%To run the summary
% subj_in={1}
% qsubcellfun(@MEGPreprocessingSummary,subj_in, 'compile', 'no',...
%         'memreq', memreq, 'timreq', timreq, 'stack', 1, 'StopOnError', false, 'backend', 'torque');

