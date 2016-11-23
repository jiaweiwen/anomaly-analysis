% This script is used to do time series clustering
% 1. Apply DTW (dynamic time warping) to measure the dissimalrity between
%    different time serier.
% 2. Use Hierachical Tree to do clustering based on DTW
% 3. Apply silhuette to determine the optimal number of clusters.

close all; clear; clc

load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat

load ../New_Data/vm_cpu_vm_cpu_mean_median_all
load ../New_Data/vm_mem_vm_mem_mean_median_all
load ../New_Data/vm_cpu_vm_mem_mean_median_all
load ../New_Data/vm_cpu_vm_mem_pair_mean_median_all

size_box_vm = size(box_vm_time_series_summary);

% Determine the maximum time length
max_time = 96;
grat_small = 900;
max_lag = 2;

mkdir('../New_Data_time_series_cluster_xcf_figure');
path = '../New_Data_time_series_cluster_xcf_figure/';

mkdir('../New_Data_XCF_for_clustering_figure');
path1 = '../New_Data_XCF_for_clustering_figure/';

size_bucket = 12;
cpu_usage_bucket = 10;
mem_usage_bucket = 10;
cpu_cap_bucket = 20;
mem_cap_bucket = 75;

possible_related_features = [];
for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary{1, box_id});
    
    % If we don't have time series
    if size_box <= 2
        continue;
    end

%     
%     if box_id > 1000
%         break;
%     end
    
    if numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1)) <= 70
        continue;
    end
   
    % First, extract all of the time series
                   
    total_vm_cpu_capacity = 0; used_mean_vm_cpu_capacity = 0;
    mean_vm_cpu_util_set = []; 
    
    total_vm_mem_capacity = 0; used_mean_vm_mem_capacity = 0;
    mean_vm_mem_util_set = []; 
    
    total_box_cpu_capacity = nanmean(box_vm_time_series_summary{1, box_id}{1, 1}(:, 3) ./ ...
                                     box_vm_time_series_summary{1, box_id}{1, 1}(:, 4) * 100);
    used_mean_box_cpu_capacity = nanmean(box_vm_time_series_summary{1, box_id}{1, 1}(:, 3));
    
    total_box_mem_capacity = nanmean(box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 3));
    used_mean_box_mem_capacity = nanmean(box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 3) .* ...
                                     box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 4) / 100);
    
    for vm_id = 2 : size_box
        total_vm_cpu_capacity = total_vm_cpu_capacity + nanmean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4) ./ ...
                                                                box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5) * 100);
        used_mean_vm_cpu_capacity = used_mean_vm_cpu_capacity + nanmean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4));
        
        mean_vm_cpu_util_set(end+1) = nanmean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5));
        
        total_vm_mem_capacity = total_vm_mem_capacity + nanmean(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 4));
        used_mean_vm_mem_capacity = used_mean_vm_mem_capacity + nanmean(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 4) .* ...
                                                                box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5) / 100);
        
        mean_vm_mem_util_set(end+1) = nanmean(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5));
    end
    mean_vm_cpu_util = nanmean(mean_vm_cpu_util_set);
    mean_vm_mem_util = nanmean(mean_vm_mem_util_set);
    
    vcap_bcap_cpu_ratio = total_vm_cpu_capacity / total_box_cpu_capacity;
    vusage_vcap_cpu_ratio = min(1, used_mean_vm_cpu_capacity / total_vm_cpu_capacity);
    vusage_bcap_cpu_ratio = min(1, used_mean_vm_cpu_capacity / total_box_cpu_capacity);
    vusage_busage_cpu_ratio = min(1, used_mean_vm_cpu_capacity / used_mean_box_cpu_capacity);
    
    vcap_bcap_mem_ratio = total_vm_mem_capacity / total_box_mem_capacity;
    vusage_vcap_mem_ratio = min(1, used_mean_vm_mem_capacity / total_vm_mem_capacity);
    vusage_bcap_mem_ratio = min(1, used_mean_vm_mem_capacity / total_box_mem_capacity);
    vusage_busage_mem_ratio = min(1, used_mean_vm_mem_capacity / used_mean_box_mem_capacity);
    
    possible_related_features(end+1, :) = [total_vm_cpu_capacity / 1024, used_mean_vm_cpu_capacity / 1024, mean_vm_cpu_util, ...
                                           total_box_cpu_capacity  / 1024, used_mean_box_cpu_capacity  / 1024, ...
                                           vcap_bcap_cpu_ratio, vusage_vcap_cpu_ratio, vusage_bcap_cpu_ratio, vusage_busage_cpu_ratio, ...
                                           total_vm_mem_capacity  / 1024, used_mean_vm_mem_capacity  / 1024, mean_vm_mem_util, ...
                                           total_box_mem_capacity  / 1024, used_mean_box_mem_capacity  / 1024, ...
                                           vcap_bcap_mem_ratio, vusage_vcap_mem_ratio, vusage_bcap_mem_ratio, vusage_busage_mem_ratio, ...
                                           size_box - 1];
end
labels_name = {'CPU Vcap', 'CPU Vusage', 'CPU Vutil', ...
               'CPU Bcap', 'CPU Busage', 'CPU Vcap / Bcap', ...
               'CPU Vusage / Vcap', 'CPU Vuage / Bcap', 'CPU Vusage / Busage', ...
               'MEM Vcap', 'MEM Vusage', 'MEM Vutil', ...
               'MEM Bcap', 'MEM Busage', 'MEM Vcap / Bcap', ...
               'MEM Vusage / Vcap', 'MEM Vuage / Bcap', 'MEM Vusage / Busage', 'Number of VMs'};
all_case_name = {'XCF: VM CPU v.s. VM CPU', 'XCF: VM RAM v.s. VM RAM', ...
                 'XCF: VM CPU v.s. VM RAM', 'XCF: VM CPU v.s. VM RAM (Pairwise)'};
% Check the correlation among the correaltion value and possible related
% features
xcf_cand = [];
all_case = {vm_cpu_vm_cpu_all, vm_mem_vm_mem_all, vm_cpu_vm_mem_all, vm_cpu_vm_mem_pair_all};
for case_id = 1 : numel(all_case)
    time_series = all_case{case_id}(:,1);
    nan_time_series = find(isnan(time_series));
    time_series(nan_time_series) = 0;
    for feature_id = 1 : numel(possible_related_features(1,:))
        temp = crosscorr(time_series, possible_related_features(:, feature_id), 1);
        xcf_cand(case_id, feature_id) = temp(2);
    end
    
    % boxplot the top related features based box plot
    for cand_id = 1 : numel(possible_related_features(1,:))
        fig = figure;
        set(fig, 'Position', [200, 200, 1200, 400]);
        % [~, cand_id] = max(xcf_cand(case_id, :));
        % First try: remove the 10%ile and 90%ile
        lower_bound = prctile(possible_related_features(:, cand_id), 5);
        upper_bound = prctile(possible_related_features(:, cand_id), 95);
        
        used_idx1 = find(possible_related_features(:, cand_id) < upper_bound);
        used_idx2 = find(possible_related_features(:, cand_id) > lower_bound);
        
        used_idx = intersect(used_idx1, used_idx2);
    
        cand_related_feature = possible_related_features(used_idx, cand_id);
        cand_time_series = time_series(used_idx);
        
        max_feature = max(cand_related_feature);
        scale = 2;
        if max_feature > 2
            cand_related_feature = ceil(log2(cand_related_feature));
        else
            interval = max_feature / 10;
            cand_related_feature = ceil(cand_related_feature / interval);
            scale = 1;
        end
        
        bucket_range = min(cand_related_feature) : max(cand_related_feature);

        xticklabel_last = {};
        if scale == 2
            for bucket_id = 1 : numel(bucket_range)
                xticklabel_last{end+1} = cellstr(num2str(bucket_range(bucket_id), '2^%d'));
            end
        else
            for bucket_id = 1 : numel(bucket_range)
                xticklabel_last{end+1} = cellstr(num2str(bucket_range(bucket_id) * interval, '%.2f'));
            end
        end
        
        xticklabel_name = [xticklabel_last{1}];
        for name_id = 2 : numel(xticklabel_last)               
            xticklabel_name(end+1) = strcat(xticklabel_last{name_id-1},'-',xticklabel_last{name_id});
           
        end
        % xticklabel_name = {'1','2','3','4','5'};
        subplot(1,2,1)
        boxplot(cand_time_series, cand_related_feature);
        xlabel(labels_name{cand_id}); ylabel(all_case_name{case_id});
        set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 13);
        set(gca, 'xlim', [0 numel(bucket_range)+1]); 
        if scale ~= 2
            xticklabel_rotate([1: numel(xticklabel_name)],45,xticklabel_name);
        else
            set(gca, 'xticklabel', xticklabel_name);
        end
        
        subplot(1,2,2)
        [X, N] = hist(cand_related_feature, bucket_range);
        X = 1 : numel(X);
        bar(X, N/sum(N) * 100);
        xlabel(labels_name{cand_id}); ylabel('Percentage (%) of Boxes');
        set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0 : 10 : 100], 'fontsize', 13);
        set(gca, 'xlim', [0 numel(X)+1]); 
        if scale ~= 2
            xticklabel_rotate([1: numel(xticklabel_name)],45,xticklabel_name);
        else
            set(gca, 'xticklabel', xticklabel_name);
        end
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'xcf_related_features_boxplot_', mat2str(case_id), '_', mat2str(cand_id)));
    end
    close all
end

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
plot(xcf_cand(1,:), 'r*-', 'linewidth', 2);
hold on
plot(xcf_cand(2,:), 'ko-', 'linewidth', 2);
hold on
plot(xcf_cand(3,:), 'm>-', 'linewidth', 2);
hold on
plot(xcf_cand(4,:), 'b<-', 'linewidth', 2);
hold on
plot(1:18, zeros(18,1), 'g--', 'linewidth', 2)
h = legend('vm cpu : vm cpu', 'vm mem : vm mem', 'vm cpu : vm mem', 'vm cpu : vm mem (pair)');
set(h, 'box', 'on', 'location', 'northeast');
set(gca, 'ylim', [-0.5 0.5]); set(gca, 'ytick', [-0.5 : 0.1 : 0.5], 'fontsize', 18);
set(gca, 'xlim', [0 18]);
xticklabel_rotate([1: feature_id], 45, labels_name);
ylabel('XCF'); xlabel('Possible Related Feature ID');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'mean_xcf_related_features'));



% all_case = {vm_cpu_vm_cpu_all, vm_mem_vm_mem_all, vm_cpu_vm_mem_all, vm_cpu_vm_mem_pair_all};
% all_case_name = {'vm_cpu_vm_cpu', 'vm_mem_vm_mem', 'vm_cpu_vm_mem', 'vm_cpu_vm_mem_pair'};
% % Box plot
% box_size_all = ceil(box_size_all / size_bucket);
% box_cpu_usage_all = ceil(box_cpu_usage_all / cpu_usage_bucket);
% box_mem_usage_all = ceil(box_mem_usage_all / mem_usage_bucket);
% box_cpu_capacity_all = ceil(box_cpu_capacity_all / 1024 / cpu_cap_bucket);
% box_ram_capacity_all = ceil(box_ram_capacity_all / 1024 / mem_cap_bucket);
% 
% box_size_all = box_size_all * size_bucket;
% box_cpu_usage_all = box_cpu_usage_all * cpu_usage_bucket;
% box_mem_usage_all = box_mem_usage_all * mem_usage_bucket;
% box_cpu_capacity_all = box_cpu_capacity_all * cpu_cap_bucket;
% box_ram_capacity_all = box_ram_capacity_all * mem_cap_bucket;
% for case_id = 1 : numel(all_case)
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,1), box_size_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Mean XCF'); xlabel('Box Consolidation');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_mean_xcf_box_size'));
%   
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,1), box_cpu_usage_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Mean XCF'); xlabel('Overall Box CPU Usage (%)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_mean_xcf_cpu_usage'));
% 
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,1), box_mem_usage_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Mean XCF'); xlabel('Overall Box RAM Usage (%)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_mean_xcf_mem_usage'));
%     
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,1), box_cpu_capacity_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Mean XCF'); xlabel('Box CPU Capacity (GHZ)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_mean_xcf_cpu_capacity'));
% 
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,1), box_ram_capacity_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Mean XCF'); xlabel('Box RAM Capacity (GB)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_mean_xcf_ram_capacity'));
% end
% 
% close all
% 
% for case_id = 1 : numel(all_case)
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,2), box_size_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Median XCF'); xlabel('Box Consolidation');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_median_xcf_box_size'));
%   
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,2), box_cpu_usage_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Median XCF'); xlabel('Overall Box CPU Usage (%)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_median_xcf_cpu_usage'));
% 
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,2), box_mem_usage_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Median XCF'); xlabel('Overall Box RAM Usage (%)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_median_xcf_mem_usage'));
%     
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,2), box_cpu_capacity_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Median XCF'); xlabel('Box CPU Capacity (GHZ)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_median_xcf_cpu_capacity'));
% 
%     fig = figure;
%     set(fig,'Position',[200, 200, 600, 400]);
%     boxplot(all_case{case_id}(:,2), box_ram_capacity_all');
%     set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
%     ylabel('Median XCF'); xlabel('Box RAM Capacity (GB)');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, all_case_name{case_id}, '_median_xcf_ram_capacity'));
% end
% 
