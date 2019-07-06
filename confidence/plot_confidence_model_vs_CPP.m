% plot model prediction for confidence and CPP amplitude behavior

data_dir = '/srv/projects/2015-02-13_TOE_Bayes/confidence'; %'/home/jan/projects/2016-04-05_SFC_ERP';


% CPP features
load(fullfile(data_dir, 'p3_peaks_equal_trials_three_classes'));

N_per_class = sum(~isnan(p3_peaks));
plot_idx_correct = 1:3;
plot_idx_incorrect = 4:6;


% plot model predictions and CPP features
% for SFC_DD
diff_class_borders = [ 0      0.09; ...
                       0.09   0.17; ...
                       0.09   0.605;];

data_subj_diffs = mean(diff_class_borders, 2);
                   
N_total = 116;

%%
figure



axes('FontSize', 15, 'position', [0.6 0.6 0.3 0.3])
hold on
errorbar(data_subj_diffs, squeeze(nanmean(p3_peaks(:,plot_idx_correct))), squeeze(nanstd(p3_peaks(:,plot_idx_correct)))./sqrt(N_per_class(plot_idx_correct)), 'o-', 'linewidth', 3, 'Color', 'g', 'MarkerFaceColor', 'g')
errorbar(data_subj_diffs, squeeze(nanmean(p3_peaks(:,plot_idx_incorrect))), squeeze(nanstd(p3_peaks(:,plot_idx_incorrect)))./sqrt(N_per_class(plot_idx_incorrect)), 'o-', 'linewidth', 3, 'Color', 'r', 'MarkerFaceColor', 'r')
set(gca, 'YAxisLocation', 'right')
xlim([0 0.4])
ylim([2.5 4.5])
% xlabel('abs. SPFD')
ylabel('CPP amplitude [\muV]')
legend('correct', 'incorrect')
% legend boxoff
set(gcf, 'renderer', 'painters')



