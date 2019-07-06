% script to preprocess the data from all SFC experiments for ERP analysis - starting point is re-referenced eye-artefact-free data
% including:
%       epoching
%       artefact removal
%       low-pass filtering
%       averaging (ERP)
%
% REQUIREMENT: import_raw_data.m and prepare_EEG_data.m must have been executed before to load raw data and create eye-blink config file

clear all;

% include SPM12
addpath(genpath('/home/jan/spm12'))

%% Some organizational stuff
% set specifics for each experiment

% SFC
subj_SFC = 2:25;
subj_SFC(ismember(subj_SFC, [3 11 23])) = []; 

project_dir_SFC = ''; % set project directoy

res_dir_SFC = fullfile(project_dir_SFC, 'processed_data', filesep);

file_mask_SFC = 'fMMMSFC_f2_locked';

% SFCwE
subj_SFCwE = 2:25;

project_dir_SFCwE = ''; % set project directoy

res_dir_SFCwE = fullfile(project_dir_SFCwE, 'processed_data');

file_mask_SFCwE = 'eog_SFCwE_SVD_corrected';

% SFC_DD
subj_SFC_DD = 6:31;  % list of subjects
subj_SFC_DD(ismember(subj_SFC_DD, [16 17])) = [];

project_dir_SFC_DD = ''; % set project directoy

res_dir_SFC_DD = fullfile(project_dir_SFC_DD, 'processed_data');

file_mask_SFC_DD = 'SFC_DD_ERP_f2_locked';

TOE_model_SFC_DD = 'SFC_DD_posteriors_simple_sigma_v_prior_bias_std_erf_direct24';

% SFC_DD_BP
subj_SFC_DD_BP = 18; % 1:20;  % list of subjects

project_dir_SFC_DD_BP = ''; % set project directoy

res_dir_SFC_DD_BP = fullfile(project_dir_SFC_DD_BP, 'processed_data');

file_mask_SFC_DD_BP = 'TfMdSFC_DD_BP_f2_locked';

TOE_model_SFC_DD_BP = 'SFC_DD_BP_posteriors_simple_sigma_v_prior_bias_std_erf_direct18';


% SFC_DD_BP_NM
subj_SFC_DD_BP_NM = 1:18;  % list of subjects

project_dir_SFC_DD_BP_NM = ''; % set project directoy

res_dir_SFC_DD_BP_NM = fullfile(project_dir_SFC_DD_BP_NM, 'processed_data');

file_mask_SFC_DD_BP_NM = 'TfMdSFC_DD_BP_NM_f2_locked';

TOE_model_SFC_DD_BP_NM = 'SFC_DD_BP_NM_posteriors_simple_sigma_v_prior_bias_std_erf_direct18';

% SFC_DD_NM
subj_SFC_DD_NM = 1:19;  % list of subjects

project_dir_SFC_DD_NM = ''; % set project directoy

res_dir_SFC_DD_NM = fullfile(project_dir_SFC_DD_NM, 'processed_data');

file_mask_SFC_DD_NM = 'TfMdSFC_DD_NM_f2_locked';

TOE_model_SFC_DD_NM = 'SFC_DD_NM_posteriors_simple_sigma_v_prior_bias_std_erf_direct17';


%% collect all the specifics from individual experiments

exp_label    = {'SFC_DD', 'SFC_DD_BP', 'SFC_DD_BP_NM', 'SFC_DD_NM'};
subjs        = {subj_SFC_DD, subj_SFC_DD_BP, subj_SFC_DD_BP_NM, subj_SFC_DD_NM};
project_dirs = {project_dir_SFC_DD, project_dir_SFC_DD_BP, project_dir_SFC_DD_BP_NM, project_dir_SFC_DD_NM};
res_dirs     = {res_dir_SFC_DD, res_dir_SFC_DD_BP, res_dir_SFC_DD_BP_NM, res_dir_SFC_DD_NM};
file_masks   = {file_mask_SFC_DD, file_mask_SFC_DD_BP, file_mask_SFC_DD_BP_NM, file_mask_SFC_DD_NM};
TOE_models   = {TOE_model_SFC_DD, TOE_model_SFC_DD_BP, TOE_model_SFC_DD_BP_NM, TOE_model_SFC_DD_NM};


% exp_label    = {'SFC', 'SFCwE'};
% subjs        = {subj_SFC, subj_SFCwE};
% project_dirs = {project_dir_SFC, project_dir_SFCwE};
% res_dirs     = {res_dir_SFC, res_dir_SFCwE};
% file_masks   = {file_mask_SFC, file_mask_SFCwE};

num_exp = length(exp_label);

%% define work steps - make sure that the file_mask matches the right conditions!!!
file_ext  = '_ERP_f2_locked';
% new_file_ext = '_ERP_f2_locked_f1_bc';

ebfname = 'ebf_conf.mat';
bt_file = 'visual_bad_trials'; 

log_filename = ['log_meta' file_ext];

% ------- processing steps -------
EPOCHS           = 0;
ARTEFACTS        = 1;
LP_FILTER        = 1;
ERP              = 1;

%% parameter settings
lpf_cutoff      = 30;       % low-pass or notch in Hz [default: 30 Hz, notch: [49 51]]
lpf_type        = 'low';    % 'stop' (requires 1x2 vector for lpf_cutoff representing stop-band) or 'low' [default: 'low']

epoching        = 2;            % where to center epoch: 1= f1_locked, 2 = f2_locked, 3 =resp_locked
epoching_labels = {'f1_locked', 'f2_locked', 'resp_locked'};

epoch_interval  = [-2250 3500]; % stim-locked: [-500 1500] resp-locked: [-1000 500] % ms

artefact_threshold = 80; % in µV <<<<<<<<<<<<<<<<<<<<<<<<<<< WATCH OUT!!!!!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><
artefact_twin      = [-1 0];

baseline       = [-100 0];     % baseline for ERP calculation in ms


% print out the specified processing steps
fprintf('\nSPECIFICATIONS:\n-------\nfilter cut-off: %dHz (%s)\nEpoching: %s\nEpoch Interval: %d to %dms\nArtefact Threshold: %dmV\nBaseline (ERP): %d to %d msec\n\n', ...
    lpf_cutoff(1), lpf_type, epoching_labels{epoching}, epoch_interval(1), epoch_interval(2), artefact_threshold, baseline)
fprintf('SELECTED OPERATIONS:\n------\n%d - Low-pass/Notch Filter\n%d - Cut Data (epoching)\n%d - Remove Artefacts\n%d - compute ERP\\n\n', ...
    LP_FILTER,EPOCHS,ARTEFACTS,ERP);
reply = input('CONTINUE??? y/n [y]: ', 's');
if isempty(reply)
    reply = 'y';
end

if strcmp(reply, 'n'), error('Aborted!'), end

% write the settings into a text file
fid = fopen([log_filename '_summary.txt'], 'w');
fprintf(fid,'\nSPECIFICATIONS:\n-------\nfilter cut-off: %dHz (%s)\nEpoching: %s\nEpoch Interval: %d to %dms\nArtefact Threshold: %dmV\nBaseline (ERP): %d to %d msec\n\n', ...
    lpf_cutoff(1), lpf_type, epoching_labels{epoching}, epoch_interval(1), epoch_interval(2), artefact_threshold, baseline);
fprintf(fid, 'SELECTED OPERATIONS:\n------\n%d - Low-pass/Notch Filter\n%d - Cut Data (epoching)\n%d - Remove Artefacts\n%d - compute ERP\\n\n', ...
    LP_FILTER,EPOCHS,ARTEFACTS,ERP);
fclose(fid);

disp('>>>>>>>>>>>>>>> Start!!!! <<<<<<<<<<<<<<<<<<<<<')

% loop over all experiments
for e = 1:num_exp
    
    exp_name = exp_label{e};
    
    % get the specifics for each experiment
    res_dir = res_dirs{e};
    cd(res_dir);            % work in result directory
    
    subj = subjs{e};
    
    file_mask = file_masks{e};
    
    % loop over all subjects
    for n=subj

        subj_ID = sprintf('%02d',n); % string of subject ID

        disp(['========================== subject ' subj_ID ' (' exp_name ') ========================='])

        subj_dir = ['subject',subj_ID];

        % change into subject dir
        cd(fullfile(res_dir, subj_dir))

        try
            load(log_filename)
        catch e_catch
            disp('Cannot load log file. Will create one.')
        end  


        %% define conditions & cut data into epochs around the trigger
        if EPOCHS
            target_file = [file_mask subj_ID];
            try
                D = spm_eeg_load(target_file);
            catch
                disp(['Epoching: No matching file found for subject ' subj_ID])
                continue
            end

            % copy to a new dataset which codes the epoching (e.g. add a f2_locked or resp_locked to the filename)
            S         = [];
            S.D       = D;
            S.outfile = [exp_name file_ext subj_ID];
            D         = spm_eeg_copy(S);

            tmp=D.events;
            for j=1:length(tmp)
                if isempty(tmp(j).value)
                    tmp(j).value=999;
                end
            end
            evt=[tmp.value];
            
            preprocess_log.strange_trigger = sum(evt == 255);

            types = {tmp.type};
            time_stamps = [tmp.time];
            block_onsets = time_stamps(strcmp(types,'Epoch'));
            preprocess_log.block_onsets = block_onsets/60.0; % block onsets in minutes

            
            % trigger values used in experiment
            % 1         - chose left
            % 2         - chose right
            % 3         - saccade onset
            % 7         - blink
            % 10        - incorrect trial
            % 11        - correct trial
            % 44        - broken fixation
            % 55        - broken fixation
            % 66        - QuaeroSys Error
            % 99        - no response recorded
            % 131-134   - f1 trigger
            % 201-220   - f2 trigger
            %
            % trigger sequences
            % ideal case: 13x 2xx 3 1/2 10/11
            % no resp:    13x 2xx 99  10

            bad_trigger = [55 66 99];
            correct_trg = [51 11];
            error_trg   = [50 10];

            if sum(ismember(evt, [301:320 401:420])) > 0
                evtlog = unique(evt(ismember(evt, [301:320 401:420])));
            else
            
                evtlog=[];
                evt = [evt repmat(99, [1 40])]; % add some dummy values in the end to make sure not to exceed the max num of evts
                for j=1:length(tmp)
                    if ismember(evt(j), [101:120 201:220])
                        % find idx of feedback trigger (10 or 11) relative to j
                        new_f1_idx = find(ismember(evt(j+1:j+40), [131:134]), 1, 'first');
                        % and extract all trigger values of current trial
                        trial_trg = evt(j+1:j+new_f1_idx);

                        if epoching == 1 % on f1

                            if any(ismember(correct_trg, trial_trg)) && ~any(ismember(bad_trigger, trial_trg))
                                tmp(j-1).value = evt(j-1)+200; % correct f1-locked
                                evtlog = [evtlog; tmp(j-1).value];

                            elseif any(ismember(error_trg, trial_trg)) && ~any(ismember(bad_trigger, trial_trg))
                                tmp(j-1).value = evt(j-1)+300; % incorrect f1-locked
                                evtlog = [evtlog; tmp(j-1).value];
                            end

                        elseif epoching == 2 % on f2;

                            if any(ismember(correct_trg, trial_trg)) && ~any(ismember(bad_trigger, trial_trg))
                                tmp(j).value = evt(j)+100; % on f2 correct trials
                                evtlog = [evtlog; tmp(j).value];

                            elseif any(ismember(error_trg, trial_trg)) && ~any(ismember(bad_trigger, trial_trg))
                                
                                tmp(j).value = evt(j)+200; % on f2 incorrect
                                evtlog = [evtlog; tmp(j).value];

                            else
                                fprintf('%d\t%d\t%d\t%d\t%d <<<<<<<<<< NO MATCH \n', evt(j:j+4))
                            end

                        elseif epoching == 3 % on saccade
                    
                            % find last saccade before evaluation
                            sacc_idx = find(trial_trg == 2, 1, 'last');
                            resp_idx = find(trial_trg == 1, 1, 'first');

                            if isempty(sacc_idx)
                                fprintf([repmat('%d\t', 1, new_f1_idx+1) ' <<<<<<<<<< NO MATCH\n'], evt(j:j+new_f1_idx))
                                
                                if ismember(51, trial_trg) && ~any(ismember(bad_trigger, trial_trg))
                                    tmp(j+resp_idx).value = evt(j)+100; % on response correct trials
                                    tmp(j+resp_idx).time = tmp(j+resp_idx).time - 0.2; % subtract the fixation period

                                elseif ismember(50, trial_trg) && ~any(ismember(bad_trigger, trial_trg))
                                    tmp(j+resp_idx).value = evt(j)+200; % on response incorrect
                                    tmp(j+resp_idx).time = tmp(j+resp_idx).time - 0.2; % subtract the fixation period
                                end
                        
                                evtlog = [evtlog; tmp(j+resp_idx).value];

                            elseif any(ismember(correct_trg, trial_trg)) && ~any(ismember(bad_trigger, trial_trg))
                                tmp(j+sacc_idx).value = evt(j)+100; % on response correct trials
                                evtlog = [evtlog; tmp(j+sacc_idx).value];

                            elseif any(ismember(error_trg, trial_trg)) && ~any(ismember(bad_trigger, trial_trg))
                                tmp(j+sacc_idx).value = evt(j)+200; % on response incorrect
                                evtlog = [evtlog; tmp(j+sacc_idx).value];
                            end
                        end
                    end
                end

                evtlog = unique(evtlog);
                D      = events(D, 1, tmp); % store recoded trigger in events
                save(D);
            end
                
            % define epochs
            S    = [];
            S.D  = D;
            S.bc = 0; % baseline correction: off

            % define trials
            S.timewin = epoch_interval;

            % correct trials
            for j=1:length(evtlog)
                S.trialdef(j).conditionlabel = num2str(evtlog(j));
                S.trialdef(j).eventtype      = 'STATUS';
                S.trialdef(j).eventvalue     = evtlog(j);
            end
            S.reviewtrials = 0;
            S.save         = 0;
            D              = spm_eeg_epochs(S);

            save(D); % save bad trials from bad block

            preprocess_log.trials_in_bad_block = sum(D.reject);
            preprocess_log.epoching.label      = epoching_labels{epoching};
            preprocess_log.epoching.interval   = epoch_interval;
        end

        %% remove bad trials
        if ARTEFACTS
            
%             target_file = ['e' exp_name file_ext subj_ID];
            target_file = ['e' file_mask subj_ID];
     
            try
                D = spm_eeg_load(target_file);
                load(fullfile(['subject' subj_ID], [bt_file subj_ID])) % load manually detected bad trial nr
            catch
                disp(['Artefact Correction: No matching file found for subject ' subj_ID])
                continue
            end

            % set manually the previously detected bad trials
%             old_name = D.fname;
%             new_name = strrep(old_name, 'eTfMd', 'aeq');
            new_name = ['aeq' exp_name file_ext subj_ID]; 

            S         = [];
            S.D       = D;
            S.outfile = new_name;
            D         = spm_eeg_copy(S);

            D = D.badtrials(bad_trials,1);

            nbad = numel(bad_trials);
            Dnbad = numel(D.badtrials);

            save(D);
            % end of manual bad trial setting

%             S                            = [];
%             S.D                          = D;
%             S.badchanthresh              = 0.2;
%             S.methods.channels           = 'EEG';
%             S.methods.fun                = 'jump_twin';
%             S.methods.settings.threshold = artefact_threshold; % µV
%             S.methods.settings.twin      = artefact_twin;
%             D                            = spm_eeg_artefact(S);
% 
            preprocess_log.total_rejected_trials = sum(D.reject);
            preprocess_log.bad_channels          = D.badchannels;

            S   = [];
            S.D = D;
            D   = spm_eeg_remove_bad_trials(S);

            if epoching == 1
                preprocess_log.trials_per_stim_pair = count_trials(D, [331:334 431:434], 1:8);
            else
                preprocess_log.trials_per_stim_pair = count_trials(D, [301:320 401:420], 1:40);
            end
        end

        %% low-pass filter or notch filter
        if LP_FILTER

            target_file = ['raeq' exp_name file_ext subj_ID];

            try
                D = spm_eeg_load(target_file);
            catch
                disp(['LP filter: No matching file found for subject ' subj_ID])
                continue
            end
            
%             % copy to a new dataset which codes the epoching (e.g. add a f2_locked or resp_locked to the filename)
%             S         = [];
%             S.D       = D;
%             S.outfile = ['rae' exp_name new_file_ext subj_ID];
%             D         = spm_eeg_copy(S);

            S      = [];
            S.D    = D;
            S.band = lpf_type;
            S.freq = lpf_cutoff;
            D      = spm_eeg_filter(S);

            preprocess_log.filter.low_pass_cutoff = lpf_cutoff;
        end

        %% Baseline correction & averaging
        if ERP
            target_file = ['fraeq' exp_name file_ext subj_ID];

            try
                D = spm_eeg_load(target_file);
            catch
                disp(['ERP computation: No matching file found for subject ' subj_ID])
                continue
            end
            
            % copy to a new dataset which codes the epoching (e.g. add a f2_locked or resp_locked to the filename)
%             old_name  = D.fname;
%             new_name  = strrep(old_name, file_ext, new_file_ext);
%             
%             S         = [];
%             S.D       = D;
%             S.outfile = new_name;
%             D         = spm_eeg_copy(S);
            
            S         = [];
            S.D       = D;
            S.timewin = baseline;
            D         = spm_eeg_bc(S);

%             S        = [];
%             S.robust = 0;
%             S.D      = D;
%             S.review = 0;
%             spm_eeg_average(S);

            preprocess_log.ERP_baseline = baseline;
        end

        save(fullfile(res_dir, subj_dir, log_filename), 'preprocess_log');

    end
end
