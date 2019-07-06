function sig_clusters = get_statistics(stats_fname, stats_var, stats_dir)

load(fullfile(stats_dir, stats_fname))

stat = eval(stats_var);

if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
    pos_cluster_pvals = [stat.posclusters(:).prob];
    pos_signif_clust = find(pos_cluster_pvals < 1);
    
    pos_time = zeros(numel(pos_signif_clust), 2);
    pos_elecs = cell(numel(pos_signif_clust),1);

    for pos_i = pos_signif_clust
        tmp = stat.time(any(stat.posclusterslabelmat == pos_i,1));
        pos_time(pos_i == pos_signif_clust,:) = [tmp(1) tmp(end)];

        pos_elecs{pos_i == pos_signif_clust} = stat.label(any(stat.posclusterslabelmat == pos_i,2));
    end
else
    pos_cluster_pvals = [];
    pos_signif_clust = [];
    pos_time = [];
    pos_elecs = [];
end
    
if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
    neg_cluster_pvals = [stat.negclusters(:).prob];
    neg_signif_clust = find(neg_cluster_pvals < 1);

    neg_time = zeros(numel(neg_signif_clust), 2);
    neg_elecs = cell(numel(neg_signif_clust),1);

    for neg_i = neg_signif_clust
        tmp = stat.time(any(stat.negclusterslabelmat == neg_i,1));  
        neg_time(neg_i == neg_signif_clust,:) = [tmp(1) tmp(end)];

        neg_elecs{neg_i == neg_signif_clust} = stat.label(any(stat.negclusterslabelmat == neg_i,2));
    end
else
    neg_cluster_pvals = [];
    neg_signif_clust = [];
    neg_time = [];
    neg_elecs = [];
end

% get peak electrodes (absolute peak + 10 % highest t-values)
[max_val, ~] = max(stat.stat,[],2);
top_elecs_pos = stat.label(find(max_val > quantile(max_val, .8)));

[~, peak_idx] = max(max_val);
peak_elec_pos = stat.label(peak_idx);

[min_val, ~] = min(stat.stat,[],2);
top_elecs_neg = stat.label(find(min_val < quantile(min_val, .2)));

[~, peak_idx] = min(min_val);
peak_elec_neg = stat.label(peak_idx);


% gather peak electrodes
sig_clusters.pos.peak_elec = peak_elec_pos;
sig_clusters.pos.top_elecs = top_elecs_pos;

sig_clusters.neg.peak_elec = peak_elec_neg;
sig_clusters.neg.top_elecs = top_elecs_neg;


% gather p-vals
sig_clusters.pos.p = pos_cluster_pvals(pos_signif_clust);
sig_clusters.neg.p = neg_cluster_pvals(neg_signif_clust);

% gather time intervals
sig_clusters.pos.twin = pos_time;
sig_clusters.neg.twin = neg_time;

% gather electrodes
sig_clusters.pos.elecs = pos_elecs;
sig_clusters.neg.elecs = neg_elecs;