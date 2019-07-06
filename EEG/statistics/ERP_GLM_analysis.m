% ERP analysis script. Compute the linear trend in each electrode for each participant
clear all

project_dir = '/set/your/path/to/data';
cd(project_dir);

TOE_dir  = '/set/your/path/to/model';

file_ID = 'bfraeSFC_DD_ERP_f2_locked';
log_fname = '_log.mat';

% load model params of Bayesian model
load(fullfile(TOE_dir, 'SFC_DD_posteriors_simple_sigma_v_prior_bias_std_erf_direct24'))
subj_ID_model = subj_ID;

subj = {'06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'}; % SFC_DD
% subj = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '12', '14', '15', '16', '17', '18', '19', '20'}; % SFC_DD_BP
% subj = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18'}; % SFC_DD_BP_NM
% subj = {'01', '02', '03', '04', '05', '07', '08', '09', '10', '11', '12', '13', '15', '16', '17', '18', '19'}; % SFC_DD_NM

N = length(subj);

%           01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
stim_set = [16 16 16 16 16 20 20 20 20 20 24 24 24 24 24 28 28 28 28 28;...
            12 14 16 18 20 16 18 20 22 24 20 22 24 26 28 24 26 28 30 32];
        
ortho_conds = [303 304 305 307 308 309 310 311 312 313 314 315 317 318];
ortho_conds = [ortho_conds ortho_conds+100];
ortho_trial_idx = [3 4 5 7 8 9 10 11 12 13 14 15 17 18];

non_zero_diff_conds = sort([1:5:20 2:5:20 4:5:20 5:5:20]);  % for studies with zero diff trials

neg_diff_conds = sort([1:5:20 2:5:20]);
pos_diff_conds = sort([4:5:20 5:5:20]);

stim_diff = diff(stim_set);

f1        = stim_set(1,:);
f2        = stim_set(2,:);
f_mean    = mean(stim_set(:));       


ORTHO_SUBSET = false;

effect_IDX = 1; % 1 - subjective stim diff, 2 - f1 mod, 3 - f2 & subjective diff, 4 - subjective diff & f2, 5 - objective diff of log freqs, 
                % 6 - obj diff & abs(obj diff), 7 - 1-stimulus probability, 8 - seperate model for pos and neg diffs

if effect_IDX == 1 || effect_IDX == 3 || effect_IDX == 4 || effect_IDX == 5 || effect_IDX == 6 || effect_IDX == 9
    ncovariates = 6;
elseif effect_IDX == 2 || effect_IDX == 7
    ncovariates = 4;
elseif effect_IDX == 8 
    ncovariates = 8;    
end



for n = 1:N
    % for each subject
    subj_ID = subj{n};
    
    fprintf('>>>>>>>> Subject %02d from %d <<<<<<<<<<<<<<<<<<<\n', n, N);
    
    cd(fullfile(['subject' subj{n}]))
    
    D            = spm_eeg_load([file_ID subj{n}]);
    chans        = D.meegchannels;
    trial_conds  = str2double(D.conditions);
    
    % take only non_zero diff trials
    non_zero_diff_trials = find(ismember(mod(trial_conds,100), non_zero_diff_conds));
    trial_conds = trial_conds(non_zero_diff_trials);
    
    ntrials      = numel(trial_conds);
    
    % initialize beta
    if n==1
        samples_ind = D.indsample(-2):D.indsample(2);
        betas = zeros(numel(chans), ncovariates, numel(samples_ind), N);
    end
    
    
    
    trial_conds_correct   = trial_conds(trial_conds < 400); % get only the condition from 1 - 16 for correcht trials
    trial_conds_incorrect = trial_conds(trial_conds > 400); % get only the condition from 1 - 16 for incorrect trials
       
    if ORTHO_SUBSET
         % get index of all ortho trials
        ortho_trials = find(ismember(trial_conds, ortho_conds));
        ntrials = numel(ortho_trials);
        
        % remove ortho trials
        for idx = 301:320
            if ~ismember(idx, ortho_conds)
                trial_conds_correct(trial_conds_correct==idx) = []; 
            end
        end

        for idx = 401:420
            if ~ismember(idx, ortho_conds)
                trial_conds_incorrect(trial_conds_incorrect==idx) = [];
            end
        end
    end
    
    ntrials_correct   = numel(trial_conds_correct);
    ntrials_incorrect = numel(trial_conds_incorrect);
    
    % Bayesian model
    subj_model_idx = find(ismember(subj_ID_model, str2double(subj_ID)));
    
    subj_log_f1             = sort(repmat(unique(subj_posterior(subj_model_idx).muX(1,:)), 1,5)); % get the possible posteriors of f1 for this subject
    subj_log_f2             = log(stim_set(2,:));
    subj_diffs_bm           = subj_log_f2 - subj_log_f1;
    subj_diff_correct       = subj_diffs_bm(trial_conds_correct-300);
    abs_subj_diff_correct   = abs(subj_diff_correct);
    subj_diff_incorrect     = subj_diffs_bm(trial_conds_incorrect-400);
    abs_subj_diff_incorrect = abs(subj_diff_incorrect);
    
    % simple model
%     w = subject_weights(n);
%     abs_subj_diff_correct   = abs(stim_set(2,trial_conds_correct-300) - (w.*stim_set(1,trial_conds_correct-300) + (1-w).*f_mean));
%     abs_subj_diff_incorrect = abs(stim_set(2,trial_conds_incorrect-400) - (w.*stim_set(1,trial_conds_incorrect-400) + (1-w).*f_mean));
%     subj_diff_correct       = (stim_set(2,trial_conds_correct-300) - (w.*stim_set(1,trial_conds_correct-300) + (1-w).*f_mean));
%     subj_diff_incorrect     = (stim_set(2,trial_conds_incorrect-400) - (w.*stim_set(1,trial_conds_incorrect-400) + (1-w).*f_mean));
    
    f1_correct              = stim_set(1,trial_conds_correct-300);
    f1_incorrect            = stim_set(1,trial_conds_incorrect-400);
    
    f2_correct              = stim_set(2,trial_conds_correct-300);
    f2_incorrect            = stim_set(2,trial_conds_incorrect-400);

    obj_diff_log_correct   = log(f2_correct) - log(f1_correct);
    obj_diff_log_incorrect = log(f2_incorrect) - log(f1_incorrect);
    
    
    obj_diff_correct   = f2_correct - f1_correct;
    obj_diff_incorrect = f2_incorrect - f1_incorrect;
    
     % load log files to compute stimulus occurences
    load(fullfile(project_dir, '..', 'log_files', [subj_ID log_fname]));
    
    f2_vals = 12:2:32;
    counts = histc(mylog.flutter(:), f2_vals);

    f2_prob          = counts./sum(counts);
    f2_prob_per_cond = zeros(1, size(stim_set,2));
    
    for f2 = f2_vals
        f2_prob_per_cond(stim_set(2,:) == f2) = 1 - f2_prob(f2==f2_vals); % Take 1 - p to get a positive correlation between amplitude and probability
    end

    f2_probs_correct   = f2_prob_per_cond(trial_conds_correct-300);
    f2_probs_incorrect = f2_prob_per_cond(trial_conds_incorrect-400);

    
    if effect_IDX == 1 % subjective stim diff
        % set up Design Matrix      
        X                          = zeros(ntrials, ncovariates);
        X(1:ntrials_correct,1)     = subj_diff_correct-mean(subj_diff_correct);
        X(ntrials_correct+1:end,2) = subj_diff_incorrect-mean(subj_diff_incorrect);
        X(1:ntrials_correct,3)     = abs_subj_diff_correct-mean(abs_subj_diff_correct);
        X(ntrials_correct+1:end,4) = abs_subj_diff_incorrect-mean(abs_subj_diff_incorrect);
        X(1:ntrials_correct,5)     = 1;
        X(ntrials_correct+1:end,6) = 1;
        
%         corr_orth_X   = spm_orth(X(:,[1 3]));
%         incorr_orth_X = spm_orth(X(:,[2 4]));
%         
%         X(:,[1 3])    = corr_orth_X;
%         X(:,[2 4])    = incorr_orth_X;
        
    elseif effect_IDX == 2 % f1 mod
        % set up design matrix        
        X                          = zeros(ntrials, ncovariates);
        X(1:ntrials_correct,1)     = f1_correct-mean(f1_correct);
        X(ntrials_correct+1:end,2) = f1_incorrect-mean(f1_incorrect);
        X(1:ntrials_correct,3)     = 1;
        X(ntrials_correct+1:end,4) = 1;
        
    elseif effect_IDX == 3 % f2 mod & subj diff mod
        
        X                          = zeros(ntrials, ncovariates);
        X(1:ntrials_correct,1)     = f2_correct-mean(f2_correct);
        X(ntrials_correct+1:end,2) = f2_incorrect-mean(f2_incorrect);
        X(1:ntrials_correct,3)     = subj_diff_correct-mean(subj_diff_correct);
        X(ntrials_correct+1:end,4) = subj_diff_incorrect-mean(subj_diff_incorrect);
        X(1:ntrials_correct,5)     = 1;
        X(ntrials_correct+1:end,6) = 1;
        
        corr_orth_X   = spm_orth(X(:,[1 3]));
        incorr_orth_X = spm_orth(X(:,[2 4]));
        
        X(:,[1 3])    = corr_orth_X;
        X(:,[2 4])    = incorr_orth_X;
        
    elseif effect_IDX == 4 % subj diff mod & f2 mod
        
        X                          = zeros(ntrials, ncovariates);
        X(1:ntrials_correct,1)     = subj_diff_correct-mean(subj_diff_correct);
        X(ntrials_correct+1:end,2) = subj_diff_incorrect-mean(subj_diff_incorrect);
        X(1:ntrials_correct,3)     = f2_correct-mean(f2_correct);
        X(ntrials_correct+1:end,4) = f2_incorrect-mean(f2_incorrect);
        X(1:ntrials_correct,5)     = 1;
        X(ntrials_correct+1:end,6) = 1;
        
        corr_orth_X   = spm_orth(X(:,[1 3]));
        incorr_orth_X = spm_orth(X(:,[2 4]));
        
        X(:,[1 3])    = corr_orth_X;
        X(:,[2 4])    = incorr_orth_X;

    elseif effect_IDX == 5 % obj diff (of log freqs) mod
        
        X                          = zeros(ntrials, ncovariates);
        X(1:ntrials_correct,1)     = obj_diff_log_correct-mean(obj_diff_log_correct);
        X(ntrials_correct+1:end,2) = obj_diff_log_incorrect-mean(obj_diff_log_incorrect);
        X(1:ntrials_correct,3)     = abs(obj_diff_log_correct)-mean(abs(obj_diff_log_correct));
        X(ntrials_correct+1:end,4) = abs(obj_diff_log_incorrect)-mean(abs(obj_diff_log_incorrect));
        X(1:ntrials_correct,5)     = 1;
        X(ntrials_correct+1:end,6) = 1;
        
    elseif effect_IDX == 6 % obj diff mod
        
        X                          = zeros(ntrials, ncovariates);
        X(1:ntrials_correct,1)     = obj_diff_correct-mean(obj_diff_correct);
        X(ntrials_correct+1:end,2) = obj_diff_incorrect-mean(obj_diff_incorrect);
        X(1:ntrials_correct,3)     = abs(obj_diff_correct)-mean(abs(obj_diff_correct));
        X(ntrials_correct+1:end,4) = abs(obj_diff_incorrect)-mean(abs(obj_diff_incorrect));
        X(1:ntrials_correct,5)     = 1;
        X(ntrials_correct+1:end,6) = 1;
        
    elseif effect_IDX == 7 % f2 probabilities
        % set up Design Matrix      
        X                          = zeros(ntrials, ncovariates);
        X(1:ntrials_correct,1)     = f2_probs_correct-mean(f2_probs_correct);
        X(ntrials_correct+1:end,2) = f2_probs_incorrect-mean(f2_probs_incorrect);
        X(1:ntrials_correct,3)     = 1;
        X(ntrials_correct+1:end,4) = 1;
        
    elseif effect_IDX == 8 % seperate distance effects for positive and negative differences
        
        trials_correct_neg = find(ismember(trial_conds_correct-300, neg_diff_conds));
        trials_correct_pos = find(ismember(trial_conds_correct-300, pos_diff_conds));
        trials_incorrect_neg = find(ismember(trial_conds_incorrect-400, neg_diff_conds));
        trials_incorrect_pos = find(ismember(trial_conds_incorrect-400, pos_diff_conds));
        
        ntrials_correct_neg    = numel(trials_correct_neg);
        ntrials_correct_pos    = numel(trials_correct_pos);
        ntrials_incorrect_neg  = numel(trials_incorrect_neg);
        ntrials_incorrect_pos  = numel(trials_incorrect_pos);
        
        % set up design matrix        
        X                                                                 = zeros(ntrials_correct_neg+ntrials_correct_pos+ntrials_incorrect_neg+ntrials_incorrect_pos, ncovariates); % not all incorrect trials cause of zero diff trials
        X(1:ntrials_correct_neg,1)                                        = subj_diff_correct(trials_correct_neg) - mean(subj_diff_correct(trials_correct_neg));          % negative difference, correct
        X(1+ntrials_correct_neg:ntrials_correct,2)                        = subj_diff_correct(trials_correct_pos) - mean(subj_diff_correct(trials_correct_pos));          % positive difference, correct
        X(1+ntrials_correct:ntrials_correct+ntrials_incorrect_neg,3)      = subj_diff_incorrect(trials_incorrect_neg) - mean(subj_diff_incorrect(trials_incorrect_neg));  % negative difference, incorrect
        X(1+ntrials_correct+ntrials_incorrect_neg:ntrials_correct+ntrials_incorrect_neg+ntrials_incorrect_pos,4)...
                                                                          = subj_diff_incorrect(trials_incorrect_pos) - mean(subj_diff_incorrect(trials_incorrect_pos));  % positive difference, incorrect
        X(1:ntrials_correct_neg,5)                                        = 1;
        X(1+ntrials_correct_neg:ntrials_correct,6)                        = 1;
        X(1+ntrials_correct:ntrials_correct+ntrials_incorrect_neg,7)      = 1;
        X(1+ntrials_correct+ntrials_incorrect_neg:ntrials_correct+ntrials_incorrect_neg+ntrials_incorrect_pos,8)...
                                                                          = 1;     
        
        trial_order = [trials_correct_neg'; trials_correct_pos'; trials_incorrect_neg'; trials_incorrect_pos'];     
        
    elseif effect_IDX == 9
        % exclude all trials with large positive subjective difference (large = upper border of 4th class: 0.09)
        trials_correct_subset = find(subj_diff_correct <= 0.09);
        trials_incorrect_subset = find(subj_diff_incorrect <= 0.09);
        
        ntrials_correct_subset   = numel(trials_correct_subset); 
        ntrials_incorrect_subset = numel(trials_incorrect_subset);
        ntrials_subset = ntrials_correct_subset + ntrials_incorrect_subset;
        
        % set up Design Matrix      
        X                          = zeros(ntrials_subset, ncovariates);
        X(1:ntrials_correct_subset,1)     = subj_diff_correct(trials_correct_subset)-mean(subj_diff_correct(trials_correct_subset));
        X(ntrials_correct_subset+1:end,2) = subj_diff_incorrect(trials_incorrect_subset)-mean(subj_diff_incorrect(trials_incorrect_subset));
        X(1:ntrials_correct_subset,3)     = abs_subj_diff_correct(trials_correct_subset)-mean(abs_subj_diff_correct(trials_correct_subset));
        X(ntrials_correct_subset+1:end,4) = abs_subj_diff_incorrect(trials_incorrect_subset)-mean(abs_subj_diff_incorrect(trials_incorrect_subset));
        X(1:ntrials_correct_subset,5)     = 1;
        X(ntrials_correct_subset+1:end,6) = 1;
        
        % index the correct trials with correct trial numbers of subset and the incorrect trials by offsetting the trial numbers of subset for incorrect trials by correct trial number of full set
        trial_order = [trials_correct_subset ntrials_correct+trials_incorrect_subset]; 
            
    end
        
    
    for c=chans
        % estimate beta for each channel
        if ORTHO_SUBSET
            y = squeeze(D(c,samples_ind,ortho_trials))';
        else
            y = squeeze(D(c,samples_ind,non_zero_diff_trials))';
        end
        
        if effect_IDX == 8 || effect_IDX == 9
            y = y(trial_order,:);
        end
        
        % estimate beta 
%         beta = (inv(X'*X)*X')*y;
        betas(c,:,:,n) = X\y;
        
    end
    
    cd('..')
    
end

save('ERP_GLM_subj_abs_diff_coh_preproc', 'betas')

disp('Done')
