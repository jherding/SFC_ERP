% predict confidence based on model - for individual studies
clear all

TOE_dir     = ''; % set path to model directory

subj_SFC          = {'02', '04', '05', '06', '07', '08', '10', '13', '14', '15', '16', '17', '19', '20', '21', '22', '24', '25'}; % 18 subjects
subj_SFCwE        = {'02', '03', '04', '06', '07', '08', '09', '10', '12', '13', '14', '15', '16', '18', '19', '20', '21', '22', '23', '24', '25'};
subj_SFC_DD       = {'06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'};
subj_SFC_DD_BP    = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '12', '14', '15', '16', '17', '18', '19', '20'}; % SFC_DD_BP
subj_SFC_DD_BP_NM = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18'}; % SFC_DD_BP_NM
subj_SFC_DD_NM    = {'01', '02', '03', '04', '05', '07', '08', '09', '10', '11', '12', '13', '15', '16', '17', '18', '19'}; % SFC_DD_NM
subj = subj_SFC_DD_NM;

% load the model parameters for individual studies
% load('/home/jan/projects/2015-02-13_TOE_Bayes/SFC_posteriors_simple_sigma_v_prior_bias_std_erf_direct18');
% load('/home/jan/projects/2015-02-13_TOE_Bayes/SFCwE_posteriors_simple_sigma_v_prior_bias_std_erf_direct23');
% load(fullfile(TOE_dir, 'SFC_DD_posteriors_simple_sigma_v_prior_bias_std_erf_direct24'))
% load(fullfile(TOE_dir, 'SFC_DD_BP_posteriors_simple_sigma_v_prior_bias_std_erf_direct18'))
% load(fullfile(TOE_dir, 'SFC_DD_BP_NM_posteriors_simple_sigma_v_prior_bias_std_erf_direct18'))
load(fullfile(TOE_dir, 'SFC_DD_NM_posteriors_simple_sigma_v_prior_bias_std_erf_direct17'))
subj_ID_model = subj_ID;


subj_sample_diffs = -0.4:0.001:0.4;
x_sample = -log(5):0.01:log(5);

Nsamples = 100000; % sample 100000 trials for each stimulus pair

confidence_incorrect = nan(numel(subj),size(subj_sample_diffs,2));
confidence_correct   = nan(numel(subj),size(subj_sample_diffs,2));

avg_perceived_evidence_correct = nan(numel(subj),size(subj_sample_diffs,2));
avg_perceived_evidence_incorrect = nan(numel(subj),size(subj_sample_diffs,2));

high_confidence_PCR = nan(numel(subj),size(subj_sample_diffs,2));
low_confidence_PCR = nan(numel(subj),size(subj_sample_diffs,2));

all_samples = zeros(numel(subj),size(subj_sample_diffs,2), Nsamples);

for n=1:length(subj)

    subj_ID = subj{n};
    
    % get subjective differences
    subj_model_idx = find(ismember(subj_ID_model, str2double(subj_ID)));
    
    if isempty(subj_model_idx)
        continue
    else
        disp(['========== subject ' subj_ID '============'])
    end

    posterior_var       = subj_posterior(subj_model_idx).muX(2,end);
    global_prior        = subj_posterior(subj_model_idx).muTheta(4);
    
    if global_prior < 0 || posterior_var < 0
        disp(['negative variances!'])
        continue
    end
    
    bias = subj_posterior(subj_model_idx).muPhi;
    
    
    
%     figure
    for d = subj_sample_diffs
        
        prob_chose_f2 = 0.5*(1 + erf(d/(sqrt(2)*sqrt(posterior_var))));  % Prob F1>F2
        prob_chose_f1 = 0.5*(1 + erf(-d/(sqrt(2)*sqrt(posterior_var))));  % Prob F1>F2

        % translate the bias (in change of probability to chose f1) into a
        % criterion shift
        criterion_shift = norminv(prob_chose_f1 + bias, d, sqrt(posterior_var)) - norminv(prob_chose_f1, d, sqrt(posterior_var));
        if isnan(criterion_shift)
            criterion_shift = 0;
        end
        
        % draw samples from percept distribution to compute confidence as probability of correct choices given the percept sample
        percept_samples = d + sqrt(posterior_var) .* randn(Nsamples,1);
        all_samples(n,d==subj_sample_diffs, :) = percept_samples;
        neg_percepts = percept_samples(percept_samples < 0);
        pos_percepts = percept_samples(percept_samples > 0);
        
        confidence_chose_f1 = (0.5*(1 + erf(abs(neg_percepts - criterion_shift) ./ (sqrt(2)*sqrt(posterior_var)) )));
        confidence_chose_f2 = (0.5*(1 + erf(abs(pos_percepts - criterion_shift) ./ (sqrt(2)*sqrt(posterior_var)) )));
        
        % compute PCRs for high and low confidence
        high_confidence_PCR(n, d == subj_sample_diffs) = mean([sign(neg_percepts(confidence_chose_f1 > 0.8)) == sign(d); sign(pos_percepts(confidence_chose_f2 > 0.8)) == sign(d)]);
        low_confidence_PCR(n, d == subj_sample_diffs) = mean([sign(neg_percepts(confidence_chose_f1 <= 0.8)) == sign(d); sign(pos_percepts(confidence_chose_f2 <= 0.8)) == sign(d)]);

        
        avg_confidence_chose_f1 = mean(confidence_chose_f1);
        avg_confidence_chose_f2 = mean(confidence_chose_f2);
        
        posterior_diff = normpdf(x_sample,d,sqrt(posterior_var));
        posterior_diff = posterior_diff/sum(posterior_diff);
        
        neg_half = posterior_diff(x_sample<criterion_shift)/sum(posterior_diff(x_sample<criterion_shift));
        pos_half = posterior_diff(x_sample>criterion_shift)/sum(posterior_diff(x_sample>criterion_shift));
        
        % compute expected value aka weighted average or center of mass on each side of decision boundary
        avg_perceived_evidence_f1 = abs(x_sample(x_sample<criterion_shift)*neg_half');
        avg_perceived_evidence_f2 = abs(x_sample(x_sample>criterion_shift)*pos_half');
        
        
        if d > 0
            confidence_correct(n,d==subj_sample_diffs)   = avg_confidence_chose_f2;
            confidence_incorrect(n,d==subj_sample_diffs) = avg_confidence_chose_f1;
            
            avg_perceived_evidence_correct(n,d==subj_sample_diffs)     = avg_perceived_evidence_f2;
            avg_perceived_evidence_incorrect(n,d==subj_sample_diffs)   = avg_perceived_evidence_f1;
        else
            confidence_correct(n,d==subj_sample_diffs)   = avg_confidence_chose_f1;
            confidence_incorrect(n,d==subj_sample_diffs) = avg_confidence_chose_f2;
            
            avg_perceived_evidence_correct(n,d==subj_sample_diffs)     = avg_perceived_evidence_f1;
            avg_perceived_evidence_incorrect(n,d==subj_sample_diffs)   = avg_perceived_evidence_f2;
        end
        
%         plot(x_sample, posterior_diff, 'r')
%         hold on
%         plot([0 0], ylim, 'k', 'linewidth', 3)
%         plot([conf_chose_f1 conf_chose_f1], ylim, 'g')
%         plot([conf_chose_f2 conf_chose_f2], ylim, 'b')
%         plot([bias bias], ylim, 'k--')
%         text(d, 0.01, sprintf('d = %1.2f', d))
%         hold off
%         waitforbuttonpress
        
    end
        
        
        
end

save('model_confidence_SFC_DD_NM', 'confidence_correct', 'confidence_incorrect', 'avg_perceived_evidence_correct', 'avg_perceived_evidence_incorrect', 'high_confidence_PCR', 'low_confidence_PCR', 'subj_sample_diffs')

% confidence
figure
hold on
plot(subj_sample_diffs(subj_sample_diffs>0), nanmean([confidence_correct(:,subj_sample_diffs>0); fliplr(confidence_correct(:,subj_sample_diffs<0))]), 'g', 'linewidth', 3)
plot(subj_sample_diffs(subj_sample_diffs>0), nanmean([confidence_incorrect(:,subj_sample_diffs>0); fliplr(confidence_incorrect(:,subj_sample_diffs<0))]), 'r', 'linewidth', 3)

figure
hold on
plot(subj_sample_diffs, nanmean(confidence_correct), 'g', 'linewidth', 3)
plot(subj_sample_diffs, nanmean(confidence_incorrect), 'r', 'linewidth', 3)

% avg perceived evidence
figure
hold on
plot(subj_sample_diffs(subj_sample_diffs>0), nanmean([avg_perceived_evidence_correct(:,subj_sample_diffs>0); fliplr(avg_perceived_evidence_correct(:,subj_sample_diffs<0))]), 'g', 'linewidth', 3)
plot(subj_sample_diffs(subj_sample_diffs>0), nanmean([avg_perceived_evidence_incorrect(:,subj_sample_diffs>0); fliplr(avg_perceived_evidence_incorrect(:,subj_sample_diffs<0))]), 'r', 'linewidth', 3)

figure
hold on
plot(subj_sample_diffs, nanmean(avg_perceived_evidence_correct), 'g', 'linewidth', 3)
plot(subj_sample_diffs, nanmean(avg_perceived_evidence_incorrect), 'r', 'linewidth', 3)

% plot PCRs for high and low confidence
figure
hold on
plot(subj_sample_diffs(subj_sample_diffs>0), nanmean([high_confidence_PCR(:,subj_sample_diffs>0); fliplr(high_confidence_PCR(:,subj_sample_diffs<0))]), 'k', 'linewidth', 3)
plot(subj_sample_diffs(subj_sample_diffs>0), nanmean([low_confidence_PCR(:,subj_sample_diffs>0); fliplr(low_confidence_PCR(:,subj_sample_diffs<0))]), 'color', [.5 .5 .5], 'linewidth', 3)


figure
hold on
plot(subj_sample_diffs, nanmean(high_confidence_PCR), 'k', 'linewidth', 3)
plot(subj_sample_diffs, nanmean(low_confidence_PCR), 'color', [.5 .5 .5], 'linewidth', 3)
