% correlate proportion of correct responses (PCRs) with subjectively perceived absolute differences (as proxy for confidence)

project_dir = '/path/to/your/data';
TOE_dir     = '/path/to/your/model';

subj_SFC          = {'02', '04', '05', '06', '07', '08', '10', '13', '14', '15', '16', '17', '19', '20', '21', '22', '24', '25'}; % 18 subjects
subj_SFCwE        = {'02', '03', '04', '06', '07', '08', '09', '10', '12', '13', '14', '15', '16', '18', '19', '20', '21', '22', '23', '24', '25'};
subj_SFC_DD       = {'06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'};
subj_SFC_DD_BP    = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '12', '14', '15', '16', '17', '18', '19', '20'}; % SFC_DD_BP
subj_SFC_DD_BP_NM = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18'}; % SFC_DD_BP_NM
subj_SFC_DD_NM    = {'01', '02', '03', '04', '05', '07', '08', '09', '10', '11', '12', '13', '15', '16', '17', '18', '19'}; % SFC_DD_NM
subj = subj_SFC_DD_NM;

file_ID = 'bfraeqSFC_DD_NM_ERP_f2_locked';

% load(TOE_dir, 'SFC_posteriors_simple_sigma_v_prior_bias_std_erf_direct18');
% load(TOE_dir, 'SFCwE_posteriors_simple_sigma_v_prior_bias_std_erf_direct23');
% load(fullfile(TOE_dir, 'SFC_DD_posteriors_simple_sigma_v_prior_bias_std_erf_direct24'))
% load(fullfile(TOE_dir, 'SFC_DD_BP_posteriors_simple_sigma_v_prior_bias_std_erf_direct18'))
% load(fullfile(TOE_dir, 'SFC_DD_BP_NM_posteriors_simple_sigma_v_prior_bias_std_erf_direct18'))
load(fullfile(TOE_dir, 'SFC_DD_NM_posteriors_simple_sigma_v_prior_bias_std_erf_direct17'))
subj_ID_model = subj_ID;

stim_set = [16 16 16 16 16 20 20 20 20 20 24 24 24 24 24 28 28 28 28 28;...
            12 14 16 18 20 16 18 20 22 24 20 22 24 26 28 24 26 28 30 32];
         
non_zero_diffs  = sort([1:5:20 2:5:20 4:5:20 5:5:20]);

ortho_trial_idx = [3 4 5 7 8 9 10 11 12 13 14 15 17 18];      

% stimulus set for SFC and SFCwE
% stim_set = [16 16 16 16 20 20 20 20 24 24 24 24 28 28 28 28;...
%             12 14 18 20 16 18 22 24 20 22 26 28 24 26 30 32];        
% non_zero_diffs = 1:16;

stim_diffs = diff(stim_set);       

subj_sample_diffs = -0.5:0.05:0.5;

blue = [0 0 1];
lblue = [0.7 .78 1];
lred  = [1 0.69 0.39];
red = [1 0 0];

rho_abs_diff_pcr  = nan(numel(subj),1);
rho_obj_diffi_pcr = nan(numel(subj),1);
p   = nan(numel(subj),1);
pcr = nan(numel(subj),size(stim_set,2));
prob_chose_f1 = nan(numel(subj),size(stim_set,2));
prob_chose_f2 = nan(numel(subj),size(stim_set,2));
prob_f1_stim = nan(numel(subj), size(stim_set,2));
prob_f1_model = nan(numel(subj), numel(subj_sample_diffs));
subj_diffs    = nan(numel(subj), size(stim_set,2));

for n=1:length(subj)

    subj_ID = subj{n};
    
    try
        D = spm_eeg_load(fullfile(project_dir, 'processed_data', ['subject' subj_ID], [file_ID subj_ID '.mat']));
        disp(['=========== subject ' subj_ID ' ==========='])
    catch
        fprintf('cannot load file from subject %s\n', subj_ID)
        continue
    end
    
    % compute PCRs per condition
    conds = str2double(D.conditions);
    
    correct_trials = find(round(conds/100) == 3);
    incorrect_trials = find(round(conds/100) == 4);
    chose_f1_trials_eq = find(round(conds/100) == 5);
    chose_f2_trials_eq = find(round(conds/100) == 6);
    
    % for each condition compute the ratio of correct responses
    for i=1:size(stim_set,2)
        
        % compute PCRs
        cond_trials = find(mod(conds, 100) == i);
        nr_total    = numel(cond_trials);                               % total number of trials for condition i
        nr_correct  = numel(intersect(cond_trials, correct_trials));    % number of correct trials for condition i

        pcr(n,i) = nr_correct/nr_total;
        
        % compute prob to chose f1 or f2 
        if stim_diffs(i) > 0
            nr_chose_f1 = nr_total - nr_correct;
            nr_chose_f2 = nr_correct;
        elseif stim_diffs(i) < 0
            nr_chose_f1 = nr_correct;
            nr_chose_f2 = nr_total - nr_correct;
        elseif stim_diffs(i) == 0
            nr_chose_f1 = numel(intersect(cond_trials, chose_f1_trials_eq));   % number of trials with choice "f1" for condition i
            nr_chose_f2 = numel(intersect(cond_trials, chose_f2_trials_eq));   % number of trials with choice "f1" for condition i         
        end
        
        prob_chose_f1(n,i) = nr_chose_f1/nr_total;
        prob_chose_f2(n,i) = nr_chose_f2/nr_total;
        
    end
    
    % get subjective differences
    subj_model_idx = find(ismember(subj_ID_model, str2double(subj_ID)));
    
    if isempty(subj_model_idx), continue, end
    
    subj_log_f1         = sort(repmat(unique(subj_posterior(subj_model_idx).muX(1,:)), 1,5)); % get the possible posteriors of f1 for this subject
    subj_log_f2         = log(stim_set(2,:));
    subj_diffs_bm       = subj_log_f2 - subj_log_f1;
    abs_subj_diffs_bm   = abs(subj_diffs_bm);

    
    % correlate the two
    [rho_abs_diff_pcr(n), p(n)] = corr(pcr(n,non_zero_diffs)', abs_subj_diffs_bm(non_zero_diffs)');
    
    [rho_obj_diffi_pcr(n), p(n)] = corr(pcr(n,non_zero_diffs)', abs(stim_diffs(non_zero_diffs))');
    
    % compute model responses
    est_params  = cell2mat({subj_posterior(subj_model_idx).muTheta});
    est_sigma   = est_params(1);
    est_v_prior = est_params(4);
    est_bias    = cell2mat({subj_posterior(subj_model_idx).muPhi});
    
    prob_f1_model(n,:) = psychometric_resp_func_model(subj_sample_diffs, est_v_prior, est_sigma, est_bias);
    
    prob_f1_stim(n,:) = psychometric_resp_func_model(subj_diffs_bm, est_v_prior, est_sigma, est_bias);
    
    subj_diffs(n,:) = subj_diffs_bm;
end

% pcr = pcr(:, non_zero_diffs);

save('behavioral_data_SFC_DD_NM', 'rho_obj_diffi_pcr', 'rho_abs_diff_pcr', 'pcr', 'subj_diffs', 'prob_f1_model', 'prob_f1_stim', 'subj_sample_diffs', 'prob_chose_f1', 'prob_chose_f2')


edges = -1:0.05:1;
count_subj_diff_rho = histc(rho_abs_diff_pcr, edges);
count_obj_diff_rho  = histc(rho_obj_diffi_pcr, edges);
figure
hold on
bar(edges, count_obj_diff_rho, 'FaceColor', [.8 .8 .8], 'EdgeColor', [.8 .8 .8])
bar(edges, count_subj_diff_rho, 'FaceColor', [1 0 0], 'EdgeColor', [1 0 0])

% 
figure
hold on
plot(stim_set(1,stim_diffs==-4), prob_f1_stim(:,stim_diffs==-4), 'b', 'linewidth', 2)
plot(stim_set(1,stim_diffs==-4), pcr(:,stim_diffs==-4), 'bs')
plot(stim_set(1,stim_diffs==-2), prob_f1_stim(:,stim_diffs==-2), 'color', lblue, 'linewidth', 2)
plot(stim_set(1,stim_diffs==-2), pcr(:,stim_diffs==-2), 's', 'color', lblue)
plot(stim_set(1,stim_diffs==2), 1-prob_f1_stim(:,stim_diffs==2), 'color', lred, 'linewidth', 2)
plot(stim_set(1,stim_diffs==2), pcr(:,stim_diffs==2), 's', 'color', lred)
plot(stim_set(1,stim_diffs==4), 1-prob_f1_stim(:,stim_diffs==4), 'r', 'linewidth', 2)
plot(stim_set(1,stim_diffs==4), pcr(:,stim_diffs==4), 'rs')

figure
hold on
plot(subj_diffs(:,stim_diffs==-4), prob_chose_f1(:,stim_diffs==-4), 'b.', 'MarkerSize', marker_size, 'MarkerFaceColor', 'b')
plot(subj_diffs(:,stim_diffs==-2), prob_chose_f1(:,stim_diffs==-2), '.', 'color', lblue, 'MarkerSize', marker_size, 'MarkerFaceColor', lblue)
plot(subj_diffs(:,stim_diffs==2), prob_chose_f1(:,stim_diffs==2), '.', 'color', lred, 'MarkerSize', marker_size, 'MarkerFaceColor', lred)
plot(subj_diffs(:,stim_diffs==4), prob_chose_f1(:,stim_diffs==4), 'r.', 'MarkerSize', marker_size, 'MarkerFaceColor', 'r')
plot(subj_sample_diffs, nanmean(prob_f1_model), 'k', 'linewidth', 2)
