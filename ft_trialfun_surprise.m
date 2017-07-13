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

% obtain info from behavioral:
% id block (1-9), id trial within block (1-76), samples per trial (1-12)
% col nr. 4 remains empty, is completed with the number of sample of the
% first element of the block
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
        if length(trials_behav) == 0
            trials_behav = new_trial;
        else
            trials_behav = [trials_behav; new_trial];
        end
    end   
end

% read the events    subject = cfg.dataset(1:3);
if isfield(cfg, 'event')
   event = cfg.event;
else
   event = ft_read_event(cfg.dataset);
end




% create vectors that will contain the samples of triggers of interest
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
if sel(end)>selgocue(end)
    sel(end)=[]; %remove last element
end
if length(selstartblock)>length(selendblock)
    selendblock = [selendblock selfeedback(end)]; % consider the last selfeedback as the end of the block
end


%Late begining
if sel(1)>selgocue(1)
    selgocue(end)=[]; % remove last element of gocue
end
if length(selstartblock)<length(selendblock)
    selstartblock = [ 1 selstartblock]; % add a beggining, so the first block could be processed
end

end_block=0; % sample of current block
first_meg_block = true;
% loop meg file blocks
for start_block = selstartblock
    nr_block = nr_block+1;
    end_block = selendblock(find(selendblock>start_block,1, 'first')); 
    
    if first_meg_block
        % calculate the id block relative to the whole session
        % regular case
        if strcmp(recording,'02') && ...
            ~((strcmp(subject,'EMB') && strcmp(session,'2')) ||...
            (strcmp(subject,'HBC') && strcmp(session,'3')) ||...
            (strcmp(subject,'JTB') && strcmp(session,'2')) ||...
            (strcmp(subject,'OMF') && strcmp(session,'1')) ||...
            (strcmp(subject,'PDP') && strcmp(session,'2')) ||...
            (strcmp(subject,'QNV') && strcmp(session,'4')) ||...
            (strcmp(subject,'TNB') && strcmp(session,'3')) ||...
            (strcmp(subject,'TSJ') && strcmp(session,'2')) )
           nr_block  = nr_block+5;
        end

        % exceptional cases

        % EMB, Session 2
        % rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-8
        if strcmp(subject,'EMB') && strcmp(session,'2') 
            if strcmp(recording,'02') 
                nr_block = 0;
            else if strcmp(recording,'03')
                nr_block = nr_block+5;
                end
            end
        end        
        % HBC, Session 3
        % rec 01: blocks 1-3, rec 02: blocks 4-8
        if strcmp(subject,'HBC') && strcmp(session,'3') 
            if strcmp(recording,'03')
                nr_block = nr_block+3;
            end
        end
        % JTB, Session 2
        % rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
        if strcmp(subject,'JTB') && strcmp(session,'2') 
            if strcmp(recording,'02') 
                nr_block = 0;
            else if strcmp(recording,'03')
                nr_block = nr_block+5;
                end
            end
        end      
        % OMF, Session 1
        % rec 01: blocks 1-4, rec 02: blocks 5-8
        if strcmp(subject,'OMF') && strcmp(session,'1') 
            if strcmp(recording,'02')
                nr_block = nr_block+4;
            end
        end       
        % PDP, Session 2
        % rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
        if strcmp(subject,'PDP') && strcmp(session,'2') 
            if strcmp(recording,'02') 
                nr_block = nr_block+1;
            else if strcmp(recording,'03')
                nr_block = nr_block+4;
                end
            end
        end
        % QNV, Session 1
        % rec 01: blocks 1-4, rec 02: blocks 5-9
        if strcmp(subject,'QNV') && strcmp(session,'4') 
            if strcmp(recording,'02')
                nr_block = nr_block+4;
            end
        end
        % TNB, Session 3
        % rec 01: block 1, rec 02: blocks 2-5, rec 03: blocks 6-9
        if strcmp(subject,'TNB') && strcmp(session,'3') 
            if strcmp(recording,'02') 
                nr_block = nr_block+1;
            else if strcmp(recording,'03')
                nr_block = nr_block+5;
                end
            end
        end
        % TSJ, Session 2
        % rec 01: blocks 1-5, rec 02: blocks 0, rec 03: blocks 6-7
        if strcmp(subject,'TSJ') && strcmp(session,'2') 
            if strcmp(recording,'02') 
                nr_block = nr_block+2;
            else if strcmp(recording,'03')
                nr_block = nr_block+5;
                end
            end
        end
    end
    % calculate the real number of trial within block based on behavioral,
    % starting from the end of sel list (it assumes that MEG data always contain
    % the end of recording, but sometimes the beggining is missing)
    
    % we need to extract from behavioral matrix only those blocks correspondent to the current recording
    trials_behav_block = trials_behav(trials_behav(:,1)==nr_block,:);   % behav trials for this block
    % create vectors with only the samples from meg data for this block
    sel_block = sel(:,sel(1,:)<=end_block & sel(1,:)>=start_block);  
    selgocue_block = selgocue(:,selgocue(1,:)<=end_block & selgocue(1,:)>=start_block);  
    selresponse_block = selresponse(:,selresponse(1,:)<=end_block & selresponse(1,:)>=start_block);  
    selfeedback_block = selfeedback(:,selfeedback(1,:)<=end_block & selfeedback(1,:)>=start_block);  
    selsampleon_block = selsampleon(:,selsampleon(1,:)<=end_block & selsampleon(1,:)>=start_block);  
    
    % if behav contains less trials than MEG (zB KSV session 1 rec 02)
    if length(sel_block) > length(trials_behav_block)
        sel_block(end)=[]; %remove last element
    end
    % if behav contains more trials than MEG (zB GSB session 1 rec 02)
    if length(sel_block) < length(trials_behav_block)
        trials_behav_block(end,:)=[]; %remove last element
    end
    trbehav_ind = length(trials_behav_block);    
    for i=length(sel_block):-1:1
        trials_behav_block(trbehav_ind,4) = sel_block(i);               % setting the sample from meg data in behav trials
        trbehav_ind = trbehav_ind - 1;
    end
    nr_trial_meg=0;
    
    % create trials based on start: mask, trigger 11
    for i=sel_block 
        nr_trial_meg = nr_trial_meg + 1;

        % obtain the trial and block info from behavioral
        trbehav_ind=find(trials_behav_block(:,4)==i,4);
        nr_trial_wblock = trials_behav_block(trbehav_ind,2);
        %nr_block = trials_behav_block(trbehav_ind,1);
        samples_trial = trials_behav_block(trbehav_ind,3);

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
        k=selgocue_block(nr_trial_meg);
        trlsamgocue = event(k).sample;
        k=selresponse_block(nr_trial_meg);
        trlsamresponse = event(k).sample;
        k=selfeedback_block(nr_trial_meg);
        trlsamfeedback = event(k).sample;  
        
        trloff=0; %TODO

        %get the number of samples per trial
        samples_trial_meg = 0;
        %i es el indice , zB 369
        %necesito la siguiente posicion de eventos de inicio
        %y asi obtendre el indice del siguiente inicio

        if nr_trial_meg+1 <= length(sel_block)  
            nextini = sel_block(nr_trial_meg+1);
        else
            nextini = selsampleon_block(length(selsampleon_block))+1;%last trial
        end
        for k = selsampleon_block
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
    first_meg_block = false;
end % end of meg file blocks loop



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
