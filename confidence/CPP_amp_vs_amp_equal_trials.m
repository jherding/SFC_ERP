% split data of each subject into trials with high and low subjectively perceived evidence (according to behavioral model) and plot CPP amplitude for correct and incorrect trials
% -> check prediction of confidence
clear all


% include SPM12
addpath(genpath('/home/jan/spm12'))

stats_dir = '/home/jan/projects/2016-04-05_SFC_ERP/stats';
data_dir = {'/home/jan/data1/EEG/SFC_DD/processed_data', ...
           '/home/jan/data1/EEG/SFC_DD_BP/processed_data', '/home/jan/data1/EEG/SFC_DD_BP_NM/processed_data', '/home/jan/data1/EEG/SFC_DD_NM/processed_data'};

TOE_dir  = '/home/jan/projects/2015-02-13_TOE_Bayes';
      
       
exp_label = {'SFC_DD', 'SFC_DD_BP', 'SFC_DD_BP_NM', 'SFC_DD_NM'};

subj_SFC_DD       = {'06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'};
subj_SFC_DD_BP    = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '12', '14', '15', '16', '17', '18', '19', '20'}; % SFC_DD_BP
subj_SFC_DD_BP_NM = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18'}; % SFC_DD_BP_NM
subj_SFC_DD_NM    = {'01', '02', '03', '04', '05', '07', '08', '09', '10', '11', '12', '13', '15', '16', '17', '18', '19'}; % SFC_DD_NM

subj = {subj_SFC_DD, subj_SFC_DD_BP, subj_SFC_DD_BP_NM, subj_SFC_DD_NM};

N_total = numel(subj_SFC_DD) + numel(subj_SFC_DD_BP) + numel(subj_SFC_DD_BP_NM) + numel(subj_SFC_DD_NM); 
alpha_crit = 0.025;

% channel of interest
% coi = {'Pz'};
coi = {'CPz'};


%           01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
stim_set = [16 16 16 16 16 20 20 20 20 20 24 24 24 24 24 28 28 28 28 28;...
            12 14 16 18 20 16 18 20 22 24 20 22 24 26 28 24 26 28 30 32];

stim_diffs = diff(stim_set);      

zero_diff_trials = 403:5:420;

% for SFC_DD
% diff_class_borders = [-0.605 -0.1; ...
%                       -0.1    0.1
%                        0.1   0.38;];
% diff_class_borders = [0 0.09; ...
%                       0.09 0.605];

diff_class_borders = [0 0.09; ...
                    0.09 0.17
                    0.17 0.605];                  
             
                   
b2r_map = b2r(-1, 1);              

color_scale = b2r_map([48 88 108 148 188 228], :);

% index of subject (considering all - manual counter)
subj_idx = 1;

class_idx     = zeros(N_total,size(diff_class_borders,1));
group_subj_diffs = cell(N_total, 1);
trial_nr_correct = nan(N_total, size(diff_class_borders,1));
trial_nr_incorrect = nan(N_total, size(diff_class_borders,1));

p3_peaks = nan(N_total, size(diff_class_borders,1)*2);

for e=1:numel(exp_label)
    
    disp(['>>>>>>>>>>>>>>>> ' exp_label{e} ' <<<<<<<<<<<<<<<<'])
    
    

    file_mask = ['bfraeq' exp_label{e} '_ERP_f2_locked'];

    exp_subj = subj{e};
    
    load(fullfile(TOE_dir, [exp_label{e} '_posteriors_simple_sigma_v_prior_bias_std_erf_direct' num2str(numel(exp_subj))]))
    subj_ID_model = subj_ID;
    
    for n = 1:length(exp_subj)
        subj_ID = exp_subj{n}; % string of subject 

        subj_file_name = fullfile(data_dir{e}, ['subject' subj_ID], [file_mask subj_ID]);

        try
            D = spm_eeg_load(subj_file_name);
            disp(['=================== Subject' subj_ID ' ===================='])
        catch
            disp(['cannot load ' subj_file_name])
            continue
        end

        
        % compute subjective differences for individual subject
        trial_conds = str2double(D.conditions);
        
        % take only non_zero diff trials
        zero_diff_trials = find(ismember(mod(trial_conds,100), 3:5:20));
%         trial_conds = trial_conds(zero_diff_trials);
        
        cond_idx   = mod(trial_conds,100);   % 1 - 20 according to stimulus pair
        choice_idx = round(trial_conds/100); % 5 - chose f2 < f1; 6 - chose f2 > f1

        ntrials      = numel(trial_conds);    


        % Bayesian model
        subj_model_idx = find(ismember(subj_ID_model, str2double(subj_ID)));

        subj_log_f1      = sort(repmat(unique(subj_posterior(subj_model_idx).muX(1,:)), 1,5)); % get the possible posteriors of f1 for this subject
        subj_log_f2      = log(stim_set(2,:));
        subj_diffs_bm    = subj_log_f2 - subj_log_f1;
        subj_diffs       = subj_diffs_bm(cond_idx);

        % binary index for correct and incorrect trials
        correct_eq_trials   = find((subj_diffs < 0 & choice_idx == 5) | (subj_diffs > 0 & choice_idx == 6));
        incorrect_eq_trials = find((subj_diffs < 0 & choice_idx == 6) | (subj_diffs > 0 & choice_idx == 5));
        
        % group the differences into classes
        subj_class_idx = zeros(size(subj_diffs));
        for d=1:length(subj_diffs)

            within_class = zeros(1,size(diff_class_borders,1));

            for b=1:size(diff_class_borders,1)
                within_class(b) = abs(subj_diffs(d)) > diff_class_borders(b,1) && abs(subj_diffs(d)) < diff_class_borders(b,2);
            end
            subj_class_idx(d) = find(within_class);
        end
                
        % compute p3 peak
        p3_data = squeeze(D.selectdata(coi, [0.5 0.8], []));
        
        for i=unique(subj_class_idx)
        
            trial_idx_correct = intersect(find(subj_class_idx == i), correct_eq_trials);            
            trial_idx_error = intersect(find(subj_class_idx == i), incorrect_eq_trials);
        
            
            p3_peaks(subj_idx, i)                            = mean(mean(p3_data(:,trial_idx_correct),2));
            p3_peaks(subj_idx, i+size(diff_class_borders,1)) = mean(mean(p3_data(:,trial_idx_error),2));            
            
        end 
        
%         str_idx = strread(num2str(zero_diff_trials),'%s');
%         
%         group_data(subj_idx, :, :) = mean(D.selectdata(coi,[-2 2],str_idx),1);
        
        subj_idx = subj_idx + 1;
    
    end
    
    
end


%% save peaks

save('/your/path/', 'p3_peaks')

%% plot the peaks of the CPP for this data set
N_per_class = sum(~isnan(p3_peaks));

figure('name', 'CPP amp vs. evidence (5 classes)')
hold on
errorbar(mean(diff_class_borders,2), squeeze(nanmean(p3_peaks(:,1:2))), squeeze(nanstd(p3_peaks(:,1:2)))./sqrt(N_per_class(1:2)), 'o-', 'linewidth', 3, 'Color', 'g', 'MarkerFaceColor', 'g')
errorbar(mean(diff_class_borders,2), squeeze(nanmean(p3_peaks(:,3:4))), squeeze(nanstd(p3_peaks(:,3:4)))./sqrt(N_per_class(3:4)), 'o-', 'linewidth', 3, 'Color', 'r', 'MarkerFaceColor', 'r')
% set(gca, 'Xtick', [1 2 3], 'XTicklabel', {'low', 'medium', 'high'})
% xlim([0.9 3.1])
xlabel('abs. SPFD')
ylabel('CPP amplitude [\muV]')
legend('correct', 'incorrect')
legend boxoff


figure('name', 'CPP peak vs signed diff')
hold on
% plot(1:size(diff_class_borders,1), median_rt_per_subj_diff, 'color', [.8 .8 .8])
% hold on
% plot(1:size(diff_class_borders,1), nanmean(median_rt_per_subj_diff), 'k', 'linewidth', 3)
errorbar(mean(diff_class_borders,2), nanmean(p3_peaks(:, 1:6)), nanstd(p3_peaks(:, 1:6))/sqrt(N_total), 'color', 'g', 'linewidth', 3 )
errorbar(mean(diff_class_borders,2), nanmean(p3_peaks(:, 7:12)), nanstd(p3_peaks(:, 7:12))/sqrt(N_total), 'color', 'r', 'linewidth', 3 )
xlabel('signed SPFD')
ylabel('CPP peak [\muV]')
legend('correct', 'incorrect')
legend boxoff

