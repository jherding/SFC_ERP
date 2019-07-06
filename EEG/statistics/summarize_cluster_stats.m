% script to read all the stats result and summarize significant clusters in
% a table

clear all

stats_dir   = '/path/to/your/stats';

out_file = 'cluster_stats_summary.txt';
fid = fopen(fullfile(stats_dir, out_file), 'w');

if ~fid
    disp('ERROR: Cannot open file')
    break
end
fprintf(fid, 'SUMMARY OF CLUSTER STATS\n\n');
fprintf(fid, '#\tTYPE\t\tP-VALUE\t\tTIME WINDOW\t\t\tELECTRODES\n');


%% get the stats for all experiments, f1 ~= f2, EVIDENCE

print_name = 'all studies, f1~=f2, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);


fprintf(fid, 'CORRECT trials\n');

% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');



% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------

fname = 'SFC_ERP_all_studies_incorrect';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_ERP_all_studies_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


%% get the stats for all experiments, f1 ~= f2, DIFFICULTY

print_name = 'all studies, f1~=f2, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');

% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------

fname = 'SFC_ERP_all_studies_incorrect';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_ERP_all_studies_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

%% get the stats for all experiments, f1 == f2, EVIDENCE

print_name = 'all studies, f1==f2 (equal trials), EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_zero_diffs_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');
% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------

fname = 'SFC_ERP_all_studies_zero_diffs_incorrect';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_ERP_all_studies_zero_diffs_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


%% get the stats for all experiments, f1 == f2, DIFFICULTY

print_name = 'all studies, f1==f2 (equal trials), DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_zero_diffs_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');

% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------

fname = 'SFC_ERP_all_studies_zero_diffs_incorrect';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_ERP_all_studies_zero_diffs_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

%% get the stats for all experiments in orthogonal subset, f1 ~= f2, EVIDENCE
fprintf(fid, '\n\n\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid, '%%%%%%%% CONTROL ANALYSIS %%%%%%%%\n');
fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');


print_name = 'all studies, f1~=f2, ORTHOGONAL subset, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_correct_ORTHO';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');

% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------

fname = 'SFC_ERP_all_studies_incorrect_ORTHO';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_ERP_all_studies_correct_vs_incorrect_evidence_ORTHO';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


%% get the stats for all experiments in orthogonal subset, f1 ~= f2, DIFFICULTY

print_name = 'all studies, f1~=f2, ORTHOGONAL subset, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_correct_ORTHO';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');
% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------

fname = 'SFC_ERP_all_studies_incorrect_ORTHO';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_ERP_all_studies_correct_vs_incorrect_difficulty_ORTHO';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

%% get the stats for all experiments with orthogonalizd regressors, f1 ~= f2, EVIDENCE


print_name = 'all studies, f1~=f2, ORTHOGONALIZED regressors, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_correct_uncorr_reg';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');

% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

%% get the stats for all experiments with orthogonalizd regressors, f1 ~= f2, DIFFICULTY


print_name = 'all studies, f1~=f2, ORTHOGONALIZED regressors, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_correct_uncorr_reg';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');

% print peak electrodes
fprintf(fid, ' overall positive peak electrode (best 10%%): %s (', cluster_stat.pos.peak_elec{1});
for e=1:length(cluster_stat.pos.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.pos.top_elecs{e});
end
fprintf(fid, ')\n');

fprintf(fid, ' overall negative peak electrode (best 10%%): %s (', cluster_stat.neg.peak_elec{1});
for e=1:length(cluster_stat.neg.top_elecs)
    fprintf(fid, '%s, ', cluster_stat.neg.top_elecs{e});
end
fprintf(fid, ')\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end



%% get the stats for all experiments, checking for a modulation by f1 during retention

print_name = 'all studies, f1~=f2, MODULATION BY f1';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_f1_mod';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_ERP_all_studies_f1_mod';
stats_var = 'difficulty_stat'; % MIND THE CONFUSING NAME!

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

fname = 'SFC_ERP_all_studies_f1_mod_correct_vs_incorrect';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


%% get the stats for individual experiments, checking for a modulation by subjective evidence

fprintf(fid, '\n\n\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid, '%%%%%%%% INDIVIDUAL EXPERIMENTS %%%%%%%%\n');
fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

% =========== SFC =====================
print_name = 'SFC, f1~=f2, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_incorrect';
stats_var = 'content_stat'; % MIND THE CONFUSING NAME!

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

print_name = 'SFC, f1~=f2, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_incorrect';
stats_var = 'difficulty_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


% =========== SFCwE =====================
print_name = 'SFCwE, f1~=f2, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFCwE_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFCwE_incorrect';
stats_var = 'content_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFCwE_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


print_name = 'SFCwE, f1~=f2, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFCwE_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFCwE_incorrect';
stats_var = 'difficulty_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFCwE_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end



% =========== SFC_DD_BP_NM =====================
print_name = 'SFC_DD_BP_NM, f1~=f2, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_NM_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_NM_incorrect';
stats_var = 'content_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_BP_NM_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end



print_name = 'SFC_DD_BP_NM, f1~=f2, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_NM_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_NM_incorrect';
stats_var = 'difficulty_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_BP_NM_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end



% =========== SFC_DD_BP =====================
print_name = 'SFC_DD_BP, f1~=f2, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_incorrect';
stats_var = 'content_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_BP_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


print_name = 'SFC_DD_BP, f1~=f2, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_BP_incorrect';
stats_var = 'difficulty_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_BP_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end



% =========== SFC_DD_NM =====================
print_name = 'SFC_DD_NM, f1~=f2, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_NM_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_NM_incorrect';
stats_var = 'content_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_NM_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


print_name = 'SFC_DD_NM, f1~=f2, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_NM_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_NM_incorrect';
stats_var = 'difficulty_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_NM_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end



% =========== SFC_DD =====================
print_name = 'SFC_DD, f1~=f2, EVIDENCE';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_correct';
stats_var = 'content_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_incorrect';
stats_var = 'content_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_correct_vs_incorrect_evidence';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


print_name = 'SFC_DD, f1~=f2, DIFFICULTY';
fprintf(fid, '\n\n======= %s =======\n', print_name);


% -------------------- CORRECT TRIALS -----------------------------
fname = 'SFC_DD_correct';
stats_var = 'difficulty_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'CORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- INCORRECT TRIALS -----------------------------
fname = 'SFC_DD_incorrect';
stats_var = 'difficulty_stat'; 

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INCORRECT trials\n');


% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% -------------------- CORRECT VS INCORRECT (INTERACTION) -----------------------------

fname = 'SFC_DD_correct_vs_incorrect_difficulty';
stats_var = 'interaction_stat';

cluster_stat = get_statistics(fname, stats_var, stats_dir);

fprintf(fid, 'INTERACTION (correct vs. incorrect trials)\n');

% positive clusters
for i=1:length(cluster_stat.pos.p)
    fprintf(fid, '%2d\tpositive \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.pos.p(i), round(cluster_stat.pos.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.pos.elecs{i})
        fprintf(fid, '%s ', cluster_stat.pos.elecs{i}{e});
    end
    fprintf(fid, '\n');
end

% negative clusters
for i=1:length(cluster_stat.neg.p)
    fprintf(fid, '%2d\tnegative \t%1.3f\t\t%4d to %4d ms\t\t', i, cluster_stat.neg.p(i), round(cluster_stat.neg.twin(i,:).*1000));
    
    for e=1:length(cluster_stat.neg.elecs{i})
        fprintf(fid, '%s ', cluster_stat.neg.elecs{i}{e});
    end
    fprintf(fid, '\n');
end


%% close the file
fclose(fid);
