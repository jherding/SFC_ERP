% demo script for model to compute confidence as in Donner et al., 2017;
% Lak et al., 2014 and as roughly described elsewhere Sanders et al. 2016
% etc.

f1 = 24;
f2 = 32; 

% stimulus vals
f1 = log(f1);          % f1 (freq on log scale)
f2 = log(f2);          % f2 (freq on log scale)

% model parameters
%--------------------
sigma     = 0.05;      % precision of encoding (variance of stimulus likelihood for f1 and f2)
gamma     = 2;      % decay parameter (multiplied with variance of f1 likelihood or posterior to model increase in memory uncertainty due to retention)
mu_global = log(22);      % mean of global prior (on log scale)

fwhm         = log(32) - log(12);
v_global     = (fwhm/(2*sqrt(2*log(2))))^2; % variance of global prior

bias = 0;

% model states
% --------------
m1  = mu_global;
v1  = v_global;

in.Gamma = 'Likelihood';

%% Up-dating the hidden states pertaining to the first stimulation

% To account for TOE, either multiply the likelihood variance by gamma....
switch in.Gamma
    case 'Likelihood'
        v_lhood1 = sigma*gamma;
    case 'Posterior'
        v_lhood1 = sigma;
end

v1p = 1/(1/v1 + 1/v_lhood1);         % variance of posterior for f1
m1p = v1p * (m1/v1 + f1/v_lhood1);   % mean of posterior for f1 (freq on log scale)

% ... or the posterior variance.
switch in.Gamma
    case 'Posterior'
        v1p = gamma*v1p;
end

%% Up-dating the hidden states pertaining to the second stimulation

v_lhood2 = sigma;   % variance of stimulus likelihood (same for f1 and f2)

%full posterior model
m2_full  = m1p;
v2_full  = v1p;
v2p_full = 1/(1/v2_full + 1/v_lhood2);         % variance of posterior for f2
m2p_full = v2p_full * (m2_full/v2_full + f2/v_lhood2);   % mean of posterior for f2 (freq on log scale)
        
%diff model (f2 - f1')
m2  = m1;    % set to arbitrary value
v2  = v1;  % set to infinty for display purposes only... prior not used
m2p  = f2;                           % mean of likelihood of f2 serves as mean of posterior for f2
v2p  = v_lhood2;                     % variance likelihood of f2 serves as variance of posterior for f2

% difference distribution
Dm = m2p - m1p;    % difference of posterior means
Sv = v1p + v2p;    % sum of variance

x_sample = -log(100):0.001:log(100);

% set difference
Dm = 0.1;
Ypost_diff  = normpdf(x_sample,Dm,sqrt(Sv));  % posterior of f1
Ypost_diff  = Ypost_diff/sum(Ypost_diff);     % ... normalized

% split distri in correct and incorrect choices and normalize them aka make
% them pdfs
neg_part = Ypost_diff(x_sample<0);
neg_part = neg_part/sum(neg_part);

pos_part = Ypost_diff(x_sample>0);
pos_part = pos_part/sum(pos_part);

avg_error_evidence   = x_sample(x_sample<0)*neg_part';
avg_correct_evidence = x_sample(x_sample>0)*pos_part';


figure
hold on
area(x_sample(x_sample<0), Ypost_diff(x_sample<0), 'Facecolor', 'r', 'EdgeColor', 'none')
area(x_sample(x_sample>0), Ypost_diff(x_sample>0), 'Facecolor', 'g', 'EdgeColor', 'none')
plot([0 0], [0 max(Ypost_diff)], 'k-', 'linewidth', 2)
plot([avg_error_evidence avg_error_evidence], [0 0.25*max(Ypost_diff)], 'color', [.5 .5 .5], 'linewidth', 5)
plot([avg_correct_evidence avg_correct_evidence], [0 0.25*max(Ypost_diff)], 'color', [.5 .5 .5], 'linewidth', 5)
xlim([-2 2])
axis off
set(gcf, 'renderer', 'painters')

%% simulate confidence for range of differences
diff_levels = 0:0.01:0.4;

% confidence for negative differences: here = incorrect choices since only
% positive difference levels tested
confidence_neg_1 = zeros(1,numel(diff_levels));
confidence_neg_2 = zeros(1,numel(diff_levels));

% confidence for positive differences: here = correct choices since only
% positive difference levels tested
confidence_pos_1 = zeros(1,numel(diff_levels));
confidence_pos_2 = zeros(1,numel(diff_levels));


avg_confidence_correct   = nan(1, numel(diff_levels));
avg_confidence_incorrect = nan(1, numel(diff_levels));

for i = 1:numel(diff_levels)
    
    Dm = diff_levels(i);
    prob_chose_f2 = 0.5*(1 + erf(Dm/(sqrt(2)*sqrt(Sv))));  % Prob F1>F2
    prob_chose_f1 = 0.5*(1 + erf(-Dm/(sqrt(2)*sqrt(Sv))));  % Prob F1>F2

    % translate the bias (in change of probability to chose f1) into a
    % criterion shift
    criterion_shift = norminv(prob_chose_f1 + bias, Dm, sqrt(Sv)) - norminv(prob_chose_f1, Dm, sqrt(Sv));

    %samples
    percepts = Dm + sqrt(Sv).*randn(1,100000);
    
    avg_confidence_correct(i) = mean(0.5*(1 + erf(abs(percepts(percepts>0)-criterion_shift) ./ (sqrt(2)*sqrt(Sv)) )));
    avg_confidence_incorrect(i) = mean(0.5*(1 + erf(abs(percepts(percepts<0)-criterion_shift) ./ (sqrt(2)*sqrt(Sv)) )));
    
    x_sample = -log(1000):0.001:log(1000);
    neg_sample = x_sample(x_sample <= 0);
    pos_sample = x_sample(x_sample >= 0);

    Ypost_diff  = normpdf(x_sample,Dm,sqrt(Sv));  % posterior of f1
    Ypost_diff  = Ypost_diff/sum(Ypost_diff);     % ... normalized

    % split distri in correct and incorrect choices and normalize them aka make
    % them pdfs
    neg_part = Ypost_diff(x_sample<=0);
    neg_part = neg_part/sum(neg_part);

    pos_part = Ypost_diff(x_sample>=0);
    pos_part = pos_part/sum(pos_part);


    % compute average perceived differences for positive and negative
    % differences AND according probabilities

    avg_neg_evidence1 = x_sample(x_sample<0)*neg_part';
    neg_evidence1 = (x_sample(x_sample<0).*neg_part).*numel(neg_part);
    avg_neg_evidence2 = x_sample(x_sample<0)*Ypost_diff(x_sample<0)';

    confidence_neg_1(i) = (0.5*(1 + erf(abs(avg_neg_evidence1-criterion_shift) ./ (sqrt(2)*sqrt(Sv)) )));
    confidence_neg_2(i) = 0.5*(1 + erf(abs(avg_neg_evidence2-criterion_shift) ./ (sqrt(2)*sqrt(Sv)) ));

    avg_pos_evidence1 = x_sample(x_sample>0)*pos_part';
    pos_evidence1 = (x_sample(x_sample>0).*neg_part).*numel(neg_part);
    avg_pos_evidence2 = x_sample(x_sample>0)*Ypost_diff(x_sample>0)';

    confidence_pos_1(i) = (0.5*(1 + erf(abs(avg_pos_evidence1-criterion_shift) ./ (sqrt(2)*sqrt(Sv)) )));
    confidence_pos_2(i) = 0.5*(1 + erf(abs(avg_pos_evidence2-criterion_shift) ./ (sqrt(2)*sqrt(Sv)) ));
end


figure('name', 'confidence')
hold on
plot(diff_levels, avg_confidence_correct, 'g', 'linewidth', 3)
plot(diff_levels, confidence_pos_1, 'k', 'linewidth', 2)
plot(diff_levels, avg_confidence_incorrect, 'r', 'linewidth', 3)
plot(diff_levels, confidence_neg_1, 'k', 'linewidth', 2)
plot([0.1 0.3], [avg_confidence_correct(diff_levels==0.1) avg_confidence_correct(diff_levels==0.2)], 'og', 'Markerfacecolor', 'g', 'MarkerSize', 10)
plot([0.1 0.3], [avg_confidence_incorrect(diff_levels==0.1) avg_confidence_incorrect(diff_levels==0.2)], 'or', 'Markerfacecolor', 'r', 'MarkerSize', 10)
ylim([0.5 1])
set(gca, 'Ytick', [.5 1], 'Xtick', [0 .2 .4], 'FontSize', 15)
set(gcf, 'renderer', 'painters')
% 




