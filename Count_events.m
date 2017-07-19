% % Script for loading CTF MEG data file for Surprise study and counting
% % number of events per block
clear, clc

addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults

cd /mnt/homes/home024/pmurphy/meg_data/surprise/  % specify path where dataset is located
%cd /mnt/homes/home024/chernandez/
cfg = [];
cfg.dataset = 'QNV-4_Surprise_20170523_01.ds';   % specify dataset to be loaded

%hdr = ft_read_header(cfg.dataset);  % read header
%dat = ft_read_data(cfg.dataset);   % read data
event = ft_read_event(cfg.dataset);  % read events

% %data = ft_preprocessing(cfg);

% Extract event markers
real_events = [];
for e = 1:length(event)
        if strcmp(event(e).type,'UPPT001')
                    real_events(end+1,1) = event(e).value;
        end
end

%%%%%%%%%%%%%%%%%%%
%%% TASK BLOCKS %%%
%%%%%%%%%%%%%%%%%%%
% Find start/end points of blocks
s_nums = find(real_events==1);
e_nums = [];
for s = 1:length(s_nums)
        e = s_nums(s)+1;
            while real_events(e)~=1 && real_events(e)~=2 && e<length(real_events)
                        e = e+1;
            end
                e_nums(end+1) = e;
end

% Loop through blocks and count events per block
num_blocks = length(s_nums);
for s = 1:length(s_nums)
        onsets(s) = length(find(real_events(s_nums(s):e_nums(s))==11));
            samples(s) = length(find(real_events(s_nums(s):e_nums(s))==21));
                postmasks(s) = length(find(real_events(s_nums(s):e_nums(s))==31));
                    
                    respcues(s) = length(find(real_events(s_nums(s):e_nums(s))==41));
                        resps(s) = length(find(real_events(s_nums(s):e_nums(s))==42 | real_events(s_nums(s):e_nums(s))==43 | real_events(s_nums(s):e_nums(s))==44));
                            feedbacks(s) = length(find(real_events(s_nums(s):e_nums(s))==51 | real_events(s_nums(s):e_nums(s))==52 | real_events(s_nums(s):e_nums(s))==53));
end

% Display event counts per block
for s = 1:length(s_nums)
        fprintf('\nTASK block %d:',s)
            fprintf('\nNumber of trial starts = %d',onsets(s))
                fprintf('\nNumber of samples = %d',samples(s))
                    fprintf('\nNumber of postmasks = %d',postmasks(s))
                        fprintf('\nNumber of response cues = %d',respcues(s))
                            fprintf('\nNumber of responses = %d',resps(s))
                                fprintf('\nNumber of feedback tones = %d\n',feedbacks(s))
end

%%%%%%%%%%%%%%%%%%%%%%
%%% SENSORY BLOCKS %%%
%%%%%%%%%%%%%%%%%%%%%%
% Find start/end points of blocks
s_nums = find(real_events==3);
e_nums = [];
for s = 1:length(s_nums)
        e = s_nums(s)+1;
            while real_events(e)~=1 && real_events(e)~=4 && e<length(real_events)
                        e = e+1;
            end
                e_nums(end+1) = e;
end

% Loop through blocks and count events per block
num_blocks = length(s_nums);
for s = 1:length(s_nums)
        Sonsets(s) = length(find(real_events(s_nums(s):e_nums(s))==13));
            Soffsets(s) = length(find(real_events(s_nums(s):e_nums(s))==14));
                Stargs(s) = length(find(real_events(s_nums(s):e_nums(s))==23));
                    Sresps(s) = length(find(real_events(s_nums(s):e_nums(s))==24));
end

% Display event counts per block
for s = 1:length(s_nums)
        fprintf('\nSENSORY block %d:',s)
            fprintf('\nNumber of samples = %d',Sonsets(s))
                fprintf('\nNumber of sample offsets = %d',Soffsets(s))
                    fprintf('\nNumber of targets = %d',Stargs(s))
                        fprintf('\nNumber of responses = %d\n',Sresps(s))
end

%%%%%%%%%%%%%%%%%%%%%%
%%% MOTOR BLOCKS %%%
%%%%%%%%%%%%%%%%%%%%%%
% Find start/end points of blocks
s_nums = find(real_events==5);
e_nums = [];
for s = 1:length(s_nums)
        e = s_nums(s)+1;
            while real_events(e)~=1 && real_events(e)~=6 && e<length(real_events)
                        e = e+1;
            end
                e_nums(end+1) = e;
end

% Loop through blocks and count events per block
num_blocks = length(s_nums);
for s = 1:length(s_nums)
        Mrests(s) = length(find(real_events(s_nums(s):e_nums(s))==16));
            Mfixs(s) = length(find(real_events(s_nums(s):e_nums(s))==26));
                Mcues(s) = length(find(real_events(s_nums(s):e_nums(s))==36));
                    Mgos(s) = length(find(real_events(s_nums(s):e_nums(s))==46));
                        Mresps(s) = length(find(real_events(s_nums(s):e_nums(s))==56 | real_events(s_nums(s):e_nums(s))==57));
end

% Display event counts per block
for s = 1:length(s_nums)
        fprintf('\nMOTOR block %d:',s)
            fprintf('\nNumber of rest cue onsets = %d',Mrests(s))
                fprintf('\nNumber of active cue onsets = %d',Mfixs(s))
                    fprintf('\nNumber of response cue onsets = %d',Mcues(s))
                        fprintf('\nNumber of go cue onsets = %d',Mgos(s))
                            fprintf('\nNumber of responses = %d\n',Mresps(s))
end