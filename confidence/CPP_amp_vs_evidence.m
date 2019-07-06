% split data of each subject into trials with high, medium and low subjectively perceived evidence (according to behavioral model) and plot CPP amplitude for correct and incorrect trials
% -> check prediction of confidence
clear all


% include SPM12
addpath(genpath('/home/jan/spm12'))

stats_dir = '/home/jan/projects/2016-04-05_SFC_ERP/stats';
data_dir = {'/home/jan/data1/EEG/DMTS_dec/processed_data', '/home/jan/data2/EEG_data/2014-10-29_SFCwE/processed_data', '/home/jan/data1/EEG/SFC_DD/processed_data', ...
           '/home/jan/data1/EEG/SFC_DD_BP/processed_data', '/home/jan/data1/EEG/SFC_DD_BP_NM/processed_data', '/home/jan/data1/EEG/SFC_DD_NM/processed_data'};

TOE_dir  = '/home/jan/projects/2015-02-13_TOE_Bayes';
      
       
exp_label = {'SFC', 'SFCwE', 'SFC_DD', 'SFC_DD_BP', 'SFC_DD_BP_NM', 'SFC_DD_NM'};

subj_SFC          = {'02', '04', '05', '06', '07', '08', '10', '13', '14', '15', '16', '17', '19', '20', '21', '22', '24', '25'}; % 18 subjects
subj_SFCwE        = {'02', '03', '04', '06', '07', '08', '09', '10', '12', '13', '14', '15', '16', '18', '19', '20', '21', '22', '23', '24', '25'}; % 21 subjects
subj_SFC_DD       = {'06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'};
subj_SFC_DD_BP    = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '12', '14', '15', '16', '17', '18', '19', '20'}; % SFC_DD_BP
subj_SFC_DD_BP_NM = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18'}; % SFC_DD_BP_NM
subj_SFC_DD_NM    = {'01', '02', '03', '04', '05', '07', '08', '09', '10', '11', '12', '13', '15', '16', '17', '18', '19'}; % SFC_DD_NM

subj = {subj_SFC, subj_SFCwE, subj_SFC_DD, subj_SFC_DD_BP, subj_SFC_DD_BP_NM, subj_SFC_DD_NM};

N_total = numel(subj_SFC) + numel(subj_SFCwE) + numel(subj_SFC_DD) + numel(subj_SFC_DD_BP) + numel(subj_SFC_DD_BP_NM) + numel(subj_SFC_DD_NM); 

n_correct_conds = 6;
n_error_conds = 6;

correct_conds = 501:506;
error_conds   = 507:512;

% channel of interest
coi = {'CPz'};



% for SFC and SFCwE
% diff_class_borders = [-0.6 -0.33; ...
%                       -0.33 -0.18; ...
%                       -0.18  -0.09; ...
%                       -0.09   0; ...
%                        0      0.09; ...
%                        0.09   0.17; ...
%                        0.17   0.33;];
                   
% for SFC_DD
diff_class_borders = [-0.605 -0.17; ...
                      -0.17  -0.09; ...
                      -0.09   0; ...
                       0      0.09; ...
                       0.09   0.17; ...
                       0.17   0.38;];     
% abs_diff_borders = [0 0.09; ...
%                     0.09 0.17
%                     0.17 0.605];
                           
plot_idx_correct = 1:size(diff_class_borders,1);
plot_idx_incorrect = size(diff_class_borders,1)+1:size(diff_class_borders,1)*2;

b2r_map = b2r(-1, 1);
color_scale = b2r_map([48 88 108 148 188 228], :);

p3_amp = nan(N_total, size(diff_class_borders,1)*2);

% index of subject (considering all - manual counter)
subj_idx = 1;

class_idx     = zeros(N_total,size(diff_class_borders,1));
group_subj_diffs = cell(N_total, 1);
trial_nr_correct = nan(N_total, size(diff_class_borders,1));
trial_nr_incorrect = nan(N_total, size(diff_class_borders,1));

p3_peaks = nan(N_total, size(diff_class_borders,1)*2); % N x diff classes * 2 because one for correct and one for incorrect

p3_amp_correct = cell(N_total,1);
p3_amp_incorrect = cell(N_total,1);

subj_diffs_correct   = cell(N_total,1);
subj_diffs_incorrect = cell(N_total,1);

for e=1:numel(exp_label)
    
    disp(['>>>>>>>>>>>>>>>> ' exp_label{e} ' <<<<<<<<<<<<<<<<'])
    
    file_mask = ['mbfrae' exp_label{e} '_ERP_f2_locked'];

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


        if subj_idx == 1
            samples_ind      = D.indsample(-2):D.indsample(2);
            group_data       = nan(N_total, numel(samples_ind), size(diff_class_borders,1));
            group_data_error = nan(size(group_data));
            t_ERP            = D.time(samples_ind);
        end
        
        % compute subjective differences for individual subject
        trial_conds = str2double(D.conditions);
        
        
        cond_idx   = mod(trial_conds,100);   % 1 - 20 according to stimulus pair

        ntrials      = numel(trial_conds);    


        % Bayesian model
        subj_model_idx = find(ismember(subj_ID_model, str2double(subj_ID)));
        
        if e <=2
            %           01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 
            stim_set = [16 16 16 16 20 20 20 20 24 24 24 24 28 28 28 28;...
                        12 14 18 20 16 18 22 24 20 22 26 28 24 26 30 32];
            subj_log_f1 = sort(repmat(unique(subj_posterior(subj_model_idx).muX(1,:)), 1,4)); % get the possible posteriors of f1 for this subject
            
        else
            
            %           01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
            stim_set = [16 16 16 16 16 20 20 20 20 20 24 24 24 24 24 28 28 28 28 28;...
                        12 14 16 18 20 16 18 20 22 24 20 22 24 26 28 24 26 28 30 32];
            subj_log_f1 = sort(repmat(unique(subj_posterior(subj_model_idx).muX(1,:)), 1,5)); % get the possible posteriors of f1 for this subject
        end
        
        subj_log_f2      = log(stim_set(2,:));
        subj_diffs_bm    = subj_log_f2 - subj_log_f1;
        subj_diffs       = subj_diffs_bm(cond_idx);

        % correct and incorrect trials
        correct_trials   = find(trial_conds < 400);
        incorrect_trials = find(trial_conds > 400);
        
        % get the p3 amplitude for each subj diff level
        subj_diffs_correct{subj_idx} = subj_diffs(correct_trials);
        subj_diffs_incorrect{subj_idx} = subj_diffs(incorrect_trials);
        
        % compute p3 peak
%         p3_data = squeeze(D.selectdata(coi, [0.5 0.8], strread(num2str(cond_idx), '%s')));
        p3_data                  = squeeze(D(D.indchannel(coi),D.indsample(0.5):D.indsample(0.8),:)); 
        p3_amp_correct{subj_idx} = mean(p3_data(:,correct_trials),1);
        p3_amp_incorrect{subj_idx} = mean(p3_data(:,incorrect_trials),1);
        
        % group the differences into classes
        subj_class_idx = zeros(size(subj_diffs));
        for d=1:length(subj_diffs)

            within_class = zeros(1,size(diff_class_borders,1));

            for b=1:size(diff_class_borders,1)
                within_class(b) = subj_diffs(d) > diff_class_borders(b,1) && subj_diffs(d) < diff_class_borders(b,2);
            end
            subj_class_idx(d) = find(within_class);
        end
                
        
        for i=unique(subj_class_idx)
        
            trial_idx_correct = intersect(find(subj_class_idx == i), correct_trials);   
            trial_idx_error = intersect(find(subj_class_idx == i), incorrect_trials);
            
            p3_peaks(subj_idx, i)                            = mean(mean(p3_data(:,trial_idx_correct),2));
            p3_peaks(subj_idx, i+size(diff_class_borders,1)) = mean(mean(p3_data(:,trial_idx_error),2));            
            
        end 
        
        subj_idx = subj_idx + 1;
    
    end
    
    
end

%% save peaks

save('/your/path/', 'p3_peaks')

%% plot CPP amplitude vs. subjectively perceived absolute evidence for this study

N_per_class = sum(~isnan(p3_peaks));

figure('name', 'CPP amp vs. evidence (5 classes)')
hold on
errorbar(mean(diff_class_borders,2), squeeze(nanmean(p3_peaks(:,plot_idx_correct))), squeeze(nanstd(p3_peaks(:,plot_idx_correct)))./sqrt(N_per_class(plot_idx_correct)), 'o-', 'linewidth', 3, 'Color', 'g', 'MarkerFaceColor', 'g')
errorbar(mean(diff_class_borders,2), squeeze(nanmean(p3_peaks(:,plot_idx_incorrect))), squeeze(nanstd(p3_peaks(:,plot_idx_incorrect)))./sqrt(N_per_class(plot_idx_incorrect)), 'o-', 'linewidth', 3, 'Color', 'r', 'MarkerFaceColor', 'r')
% set(gca, 'Xtick', [1 2 3], 'XTicklabel', {'low', 'medium', 'high'})
% xlim([0.9 3.1])
xlabel('abs. SPFD')
ylabel('CPP amplitude [\muV]')
legend('correct', 'incorrect')
legend boxoff


figure('name', 'CPP amp vs. evidence (detailed)')
hold on
plot([subj_diffs_correct{:}], [p3_amp_correct{:}], '.g')
plot([subj_diffs_incorrect{:}], [p3_amp_incorrect{:}], '.r')
set(gca, 'Xtick', [1 2 3], 'XTicklabel', {'low', 'medium', 'high'})
% xlim([0.9 3.1])
xlabel('abs. SPFD')
ylabel('PCR')
legend('high', 'medium', 'low')
legend boxoff

