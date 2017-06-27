function [trl, event] = ft_trialfun_surprise(cfg)

trials_behav = [];
dir_behav = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Data/';

%triggers of interest
init_tr = [11];             %onset og pre-sequence mask
gocue_tr = [41];            %cue for response
response_tr = [42,43,44];   %left, right, bad response
feedback_tr = [51,52,53];   %correct, miss, bad response
end_tr = [61];              %offset of feedback/onset of break period
sampleon_tr = [21];         %onset individual sample
startblock_tr = [1];        %start block
endblock_tr = [2];          %end block

trl = [];
% read the header, contains the sampling frequency
hdr = ft_read_header(cfg.dataset);

%read info from ds name
%example: cfg.dataset = 'TFD-1_Surprise_20170315_01.ds';
% if isfield(cfg, 'dataset')
%     subject = cfg.dataset(1:3); %TODO convert to 1x1 
%     session = cfg.dataset(5);
%     session_date = cfg.dataset(16:23);; %TODO convert to 1x1
%     session_part = cfg.dataset(26);
% end

subject = cfg.subj;
session = cfg.sess;
recording = cfg.rec;

%obtain info from behavioral
dir_behaviour = [dir_behav,subject,'/S',session,'/Behaviour/'];
dir_samples = [dir_behav,subject,'/S',session,'/Sample_seqs/'];
blockbehav = 0;
samples_x_trial = 0;

subjfile = dir([dir_behaviour,'*','.mat']);
sessfile = dir([dir_samples,'*','.mat']);
    
for i=1:length(subjfile),               %block
    load([dir_behaviour, subjfile(i).name]);
    load([dir_samples, sessfile(i).name]);
    nrtrials_block = length(Behav);
    blockbehav = blockbehav + 1;
    for j=1:length(Behav),              %trial
        samples_x_trial = sum(~isnan(stimIn(j,:)),2); 
        new_trial = [blockbehav j samples_x_trial 0];
        if (strcmp(recording,'01') && i<6)  || (strcmp(recording,'02') && i>5)
            if (i == 1 || i == 6) && j==1
              trials_behav = new_trial;
            else
              trials_behav = [trials_behav; new_trial];
            end
        end
    end   
end

% read the events    subject = cfg.dataset(1:3);
if isfield(cfg, 'event')
   event = cfg.event;
else
   event = ft_read_event(cfg.dataset);
end

sel = true(1, length(event));  
selinit = true(1, length(event));
selgocue =  true(1, length(event));
selresponse = true(1, length(event));
selfeedback = true(1, length(event));
selend = true(1, length(event)); 
selsampleon = true(1, length(event)); 
selstartblock = true(1, length(event)); 
selendblock = true(1, length(event));

% select all events of the trigger channel
if isfield(cfg.trialdef, 'eventtype') && ~isempty(cfg.trialdef.eventtype)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).type, cfg.trialdef.eventtype);
  end
end

% select all events for the triggers of interest
if isfield(cfg.trialdef, 'eventvalue') && ~isempty(cfg.trialdef.eventvalue)
  for i=1:numel(event)
    selinit(i) = sel(i) && ismatch(event(i).value, init_tr);
    selgocue(i) = sel(i) && ismatch(event(i).value, gocue_tr);
    selresponse(i) = sel(i) && ismatch(event(i).value, response_tr);
    selfeedback(i) = sel(i) && ismatch(event(i).value, feedback_tr);
    selend(i) = sel(i) && ismatch(event(i).value, end_tr); 
    selsampleon(i) = sel(i) && ismatch(event(i).value, sampleon_tr);
    selstartblock(i) = sel(i) && ismatch(event(i).value, startblock_tr);
    selendblock(i) = sel(i) && ismatch(event(i).value, endblock_tr);
    sel(i) = sel(i) && ismatch(event(i).value, cfg.trialdef.eventvalue);
  end
end

% convert from boolean vector into a list of indices of the events of interest
sel = find(sel);
selinit = find(selinit);
selgocue = find(selgocue);
selresponse = find(selresponse);
selfeedback = find(selfeedback);
selend = find(selend);
selsampleon = find(selsampleon);
selstartblock = find(selstartblock);
selendblock = find(selendblock);

nr_trial=0;         %to count trial within dataset
nr_trial_wblock=0;  %to count trial within block
nr_block=0;

%End interrupted
if(sel(length(sel)))>selfeedback(length(selfeedback))
    sel(length(sel))=[]; %remove last element
end
%Late begining
if(sel(length(1)))>selfeedback(length(1))
    selfeedback(length(selfeedback))=[]; %remove last element
end

%calculate the real number of trial within block based on behavioral,
%starting from the end of sel list (it assumes that MEG data always contain
%the end of recording, but sometimes the beggining is missing)
trbehav_ind = length(trials_behav); 
for i=length(sel):-1:1
    trials_behav(trbehav_ind,4) = sel(i);
    trbehav_ind = trbehav_ind - 1;
end
nr_trial_meg=0;
%create trials based on start: mask, trigger 11
for i=sel 
    nr_trial_meg = nr_trial_meg + 1;

    %obtain the trial and block info from behavioral
    trbehav_ind=find(trials_behav(:,4)==i,4);
    nr_trial = trbehav_ind;
    nr_trial_wblock = trials_behav(trbehav_ind,2);
    nr_block = trials_behav(trbehav_ind,1);
    samples_trial = trials_behav(trbehav_ind,3);
    
    % determine where the trial starts with respect to the event
    if ~isfield(cfg.trialdef, 'prestim')
        trloff = event(i).offset;
        trlbeg = event(i).sample;
    else
    % shift the begin sample with the specified amount
        trlsammask = event(i).sample;
        trlbeg = event(i).sample + round(-cfg.trialdef.prestim*hdr.Fs);% ch here calculation of sample based on sg and freq sample (1200). verify if this is the correct way to calcule sample from time!!
    end
    
    %obtain the index for the triggers of interest 
    k=selgocue(nr_trial_meg);
    trlsamgocue = event(k).sample;
    k=selresponse(nr_trial_meg);
    trlsamresponse = event(k).sample;
    k=selfeedback(nr_trial_meg);
    trlsamfeedback = event(k).sample;  
    %k=selend(i);
    %trlsamend = event(k).sample;
    
    trloff=0; %TODO
    
    %get the number of samples per trial
    samples_trial_meg = 0;
    %i es el indice , zB 369
    %necesito la siguiente posicion de eventos de inicio
    %y asi obtendre el indice del siguiente inicio
    
    if nr_trial_meg+1 <= length(sel)  
        nextini = sel(nr_trial_meg+1);
    else
        nextini = selsampleon(length(selsampleon))+1;%last trial
    end
    for k = selsampleon
      if k<nextini && k>i 
          samples_trial_meg = samples_trial_meg +1;
      end
    end
    
    % determine the number of samples that has to be read (excluding the begin sample)
    %TODO: verify, they do this in the general_fun, maybe we don't
    %need that here

    %calculate duration and end (here we use the go cue)
    trlend = trlsamgocue;
    trldur = trlsamgocue - trlbeg; %in samples

    %TODO: obtain from behavioral data
    %responseTime, relative to the go on cue
    %accuracy
    trldur_frommask = trlsamgocue - trlsammask;
    trldur_sg = trldur / hdr.Fs;
    trldur_frommask_sg = trldur_frommask / hdr.Fs;
    trloff = round(-cfg.trialdef.prestim*hdr.Fs);
    
    new_trial = [trlbeg,trlend,trloff,nr_block,nr_trial_wblock, trldur, trldur_frommask, trldur_sg, trldur_frommask_sg, trlsammask, trlsamgocue,trlsamresponse,trlsamfeedback,samples_trial,samples_trial_meg];

    %transform the rest of variables so everything is consistent
    %subject,session,session_date,session_part,block,nr_trial,
    %trlsammask, trlsamgocue,trlsamresponse,trlsamfeedback,samples_trial
    
    if length(trl) == 0
      trl = new_trial;
    else
      trl = [trl; new_trial];
    end

    cfg.trl = trl;               
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION returns true if x is a member of array y, regardless of the class of x and y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = ismatch(x, y)
if isempty(x) || isempty(y)
  s = false;
elseif ischar(x) && ischar(y)
  s = strcmp(x, y);
elseif isnumeric(x) && isnumeric(y)
  s = ismember(x, y);
elseif ischar(x) && iscell(y)
  y = y(strcmp(class(x), cellfun(@class, y, 'UniformOutput', false)));
  s = ismember(x, y);
elseif isnumeric(x) && iscell(y) && all(cellfun(@isnumeric, y))
  s = false;
  for i=1:numel(y)
    s = s || ismember(x, y{i});
  end
else
  s = false;
end
