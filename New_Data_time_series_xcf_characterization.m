% This script is used to do time series clustering
% 1. Apply DTW (dynamic time warping) to measure the dissimalrity between
%    different time serier.
% 2. Use Hierachical Tree to do clustering based on DTW
% 3. Apply silhuette to determine the optimal number of clusters.

close all; clear; clc

load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat

size_box_vm = size(box_vm_time_series_summary);

% Determine the maximum time length
max_time = 96;
grat_small = 900;
max_lag = 2;

mkdir('../New_Data_time_series_cluster_xcf_figure');
path = '../New_Data_time_series_cluster_xcf_figure/';

mkdir('../New_Data_XCF_for_clustering_figure');
path1 = '../New_Data_XCF_for_clustering_figure/';

box_cpu_vm_cpu_all = []; box_mem_vm_mem_all = [];
box_cpu_vm_mem_all = []; box_mem_vm_cpu_all = [];
vm_cpu_vm_mem_all = []; 
vm_cpu_vm_mem_pair_all = [];
vm_cpu_vm_cpu_all = [];
vm_mem_vm_mem_all = [];

box_size_all = [];

box_cpu_usage_all = [];
box_mem_usage_all = [];

box_cpu_capacity_all = [];
box_ram_capacity_all = [];

size_bucket = 12;
cpu_usage_bucket = 10;
mem_usage_bucket = 10;
cpu_cap_bucket = 20;
mem_cap_bucket = 75;

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
    time = box_vm_time_series_summary{1, box_id}{1, 1}(:,1) - box_vm_time_series_summary{1, box_id}{1, 1}(1,1);
    pm_id = box_vm_time_series_summary{1, box_id}{1, 1}(1,2);
    %time_series_cpu = box_vm_time_series_summary{1, box_id}{1, 1}(:,4);
    all_time_series = [box_vm_time_series_summary{1, box_id}{1, 1}(:,4), ...
                       box_vm_time_series_summary_mem{1, box_id}{1, 1}(:,4)];
    for vm_id = 2 : size_box
%         time_series_cpu = [time_series_cpu, ...
%                   box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5)];
        all_time_series = [all_time_series, ...
                  box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5), ...
                  box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5)];
    end
    
    % Second, calculate the dissimilarities among these time series using
    % cross-correlation. Need to notice that, here we change the xcf into
    % the range of [0,1], which makes it more sense to measure the
    % dissimilarity
    box_cpu_vm_cpu = []; box_mem_vm_mem = [];
    box_cpu_vm_mem = []; box_mem_vm_cpu = [];
    vm_cpu_vm_mem = []; 
    vm_cpu_vm_cpu = [];
    vm_mem_vm_mem = [];
    vm_cpu_vm_mem_pair = [];
    total_metric_num = size(all_time_series, 2);
    d = zeros(total_metric_num); d_half = zeros(total_metric_num);
    cross_all = ones(total_metric_num);
%     for row = 2 : total_metric_num
%         for col = 1 : row -1
%             % assume there is no shift among the time series
%             temp = crosscorr(all_time_series(:,row), all_time_series(:, col), max_lag);
% %             d(row, col) = 1 - max(abs(temp));
% %             cross_all(row, col) = max(temp);
% %             cross_all(col, row) = max(temp);
% %             d(col, row) = d(row, col);
% %             d_half(row, col) = d(row, col);
%             
%             if col == 1 && rem(row,2) == 1
%                 box_cpu_vm_cpu(end+1) = max(abs(temp));
%             elseif col == 2 && rem(row,2) == 0
%                 box_mem_vm_mem(end+1) = max(abs(temp));
%             elseif col == 1 && rem(row, 2) == 0 && row ~= 2
%                 box_cpu_vm_mem(end+1) = max(abs(temp));
%             elseif col == 2 && rem(row,2) == 1 
%                 box_mem_vm_cpu(end+1) = max(abs(temp));           
%             elseif rem(col,2) == 1 && col ~= 1 && rem(row,2) == 0
%                 vm_cpu_vm_mem(end+1) = max(abs(temp));
%                 % Record the pair-wise correlation-coefficient
%                 if row == col + 1
%                     vm_cpu_vm_mem_pair(end+1) = max(abs(temp));
%                 end
%             elseif rem(col,2) == 1 && col ~= 1 && rem(row,2) == 1
%                 vm_cpu_vm_cpu(end+1) = max(abs(temp));
%             elseif rem(col,2) == 0 && col ~= 2 && rem(row,2) == 0
%                 vm_mem_vm_mem(end+1) = max(abs(temp));
%             end   
%         end
%     end    
%     
%     if numel(box_cpu_vm_cpu) ~= 0
%         box_cpu_vm_cpu_all(end+1, 1:2) = [nanmean(box_cpu_vm_cpu), nanmedian(box_cpu_vm_cpu)];
% 
%         box_mem_vm_mem_all(end+1, 1:2) = [nanmean(box_mem_vm_mem), nanmedian(box_mem_vm_mem)];
%     
%         box_cpu_vm_mem_all(end+1, 1:2) = [nanmean(box_cpu_vm_mem), nanmedian(box_cpu_vm_mem)];
%    
%         box_mem_vm_cpu_all(end+1, 1:2) = [nanmean(box_mem_vm_cpu), nanmedian(box_mem_vm_cpu)];
%     
%         vm_cpu_vm_mem_all(end+1, 1:2) = [nanmean(vm_cpu_vm_mem), nanmedian(vm_cpu_vm_mem)];
% 
%         vm_cpu_vm_mem_pair_all(end+1, 1:2) = [nanmean(vm_cpu_vm_mem_pair), nanmedian(vm_cpu_vm_mem_pair)];
%     
%         vm_cpu_vm_cpu_all(end+1, 1:2) = [nanmean(vm_cpu_vm_cpu), nanmedian(vm_cpu_vm_cpu)];
% 
%         vm_mem_vm_mem_all(end+1, 1:2) = [nanmean(vm_mem_vm_mem), nanmedian(vm_mem_vm_mem)];
% 
%         box_size_all(end+1) = size_box - 1;
% 
%         box_cpu_usage_all(end+1) = nanmean(box_vm_time_series_summary{1, box_id}{1,1}(:,end));
%         box_mem_usage_all(end+1) = nanmean(box_vm_time_series_summary_mem{1, box_id}{1,1}(:, end));
% 
%         box_util_non_zero_idx = find(box_vm_time_series_summary{1, box_id}{1,1}(:,end) > 0);
%         box_cpu_capacity_all(end+1) = nanmean(box_vm_time_series_summary{1, box_id}{1,1}(box_util_non_zero_idx,end -1) ./ ...
%                                               box_vm_time_series_summary{1, box_id}{1,1}(box_util_non_zero_idx,end)) * 100;
%         box_ram_capacity_all(end+1) = nanmean(box_vm_time_series_summary_mem{1, box_id}{1,1}(:,end -1));
%     end
        
end
% save('../New_Data/box_cpu_vm_cpu_mean_median_all', 'box_cpu_vm_cpu_all');
% save('../New_Data/box_mem_vm_mem_mean_median_all', 'box_mem_vm_mem_all');
% save('../New_Data/box_cpu_vm_mem_mean_median_all', 'box_cpu_vm_mem_all');
% save('../New_Data/box_mem_vm_cpu_mean_median_all', 'box_mem_vm_cpu_all');
% save('../New_Data/vm_cpu_vm_cpu_mean_median_all', 'vm_cpu_vm_cpu_all');
% save('../New_Data/vm_mem_vm_mem_mean_median_all', 'vm_mem_vm_mem_all');
% save('../New_Data/vm_cpu_vm_mem_mean_median_all', 'vm_cpu_vm_mem_all');
% save('../New_Data/vm_cpu_vm_mem_pair_mean_median_all', 'vm_cpu_vm_mem_pair_all');
% save('../New_Data/box_cpu_usage_all','box_cpu_usage_all');
% save('../New_Data/box_mem_usage_all','box_mem_usage_all');
% save('../New_Data/box_size_all','box_size_all');
% save('../New_Data/box_cpu_capacity_all', 'box_cpu_capacity_all');
% save('../New_Data/box_ram_capacity_all','box_ram_capacity_all');

% % Plot the cdf of XCF
% path = '../New_Data_VM_BOX_Figure/';
% fig = figure;
% set(fig,'Position',[200, 200, 600, 400]);
% set(gca, 'fontsize', 18);
% [f_1, x_1] = ecdf(box_cpu_vm_cpu_all);
% [f_2, x_2] = ecdf(box_mem_vm_mem_all);
% 
% [f_3, x_3] = ecdf(box_cpu_vm_mem_all);
% [f_4, x_4] = ecdf(box_mem_vm_cpu_all);
% 
% plot(x_1, f_1, 'k-', 'linewidth', 2)
% hold on
% plot(x_2, f_2, 'r--', 'linewidth', 2)
% hold on
% plot(x_3, f_3, 'b-.', 'linewidth', 2);
% hold on
% plot(x_4, f_4, 'g-', 'linewidth', 2);
% % hold on
% % plot(x_5, f_5, 'g-', 'linewidth', 1);
% h = legend('BOX:CPU-VM:CPU ', 'BOX:MEM-VM:MEM', 'BOX:CPU-VM:MEM', 'BOX:MEM-VM:CPU');%, 'VM:CPU-VM:CPU', ...
%           % 'VM:MEM-VM:MEM', 'VM:CPU-VM:MEM');
% set(h, 'box','on','location','northwest','fontsize',18);
% set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1]);
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% % set(gca, 'xscale', 'log');
% ylabel('CDF', 'fontsize', 18); 
% xlabel('XCF', 'fontsize', 18);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'VM_BOX_CPU_MEM_MAX_XCF'));
% 

% CDF of mean
% path = '../New_Data_VM_BOX_Figure/';
% fig = figure;
% set(fig,'Position',[200, 200, 600, 400]);
% set(gca, 'fontsize', 18);
% [f_5, x_5] = ecdf(vm_cpu_vm_cpu_all(:,1));
% [f_6, x_6] = ecdf(vm_mem_vm_mem_all(:,1));
% [f_7, x_7] = ecdf(vm_cpu_vm_mem_all(:,1));
% [f_8, x_8] = ecdf(vm_cpu_vm_mem_pair_all(:,1));
% plot(x_5, f_5, 'k-', 'linewidth', 2)
% hold on
% plot(x_6, f_6, 'r--', 'linewidth', 2)
% hold on
% plot(x_7, f_7, 'b-.', 'linewidth', 2);
% hold on
% plot(x_8, f_8, 'm:', 'linewidth', 2);
% h = legend('VM:CPU-VM:CPU', 'VM:RAM-VM:RAM', 'VM:CPU-VM:RAM', 'VM:CPU-VM:RAM (Pair)');
% set(h, 'box','on','location','southeast','fontsize',15);
% set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1]);
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% % set(gca, 'xscale', 'log');
% ylabel('CDF', 'fontsize', 18); 
% xlabel('Mean XCF', 'fontsize', 18);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'VM_VM_CPU_MEM_mean_XCF'));
% 
% % CDF of median
% fig = figure;
% set(fig,'Position',[200, 200, 600, 400]);
% set(gca, 'fontsize', 18);
% [f_5, x_5] = ecdf(vm_cpu_vm_cpu_all(:,2));
% [f_6, x_6] = ecdf(vm_mem_vm_mem_all(:,2));
% [f_7, x_7] = ecdf(vm_cpu_vm_mem_all(:,2));
% [f_8, x_8] = ecdf(vm_cpu_vm_mem_pair_all(:,2));
% plot(x_5, f_5, 'k-', 'linewidth', 2)
% hold on
% plot(x_6, f_6, 'r--', 'linewidth', 2)
% hold on
% plot(x_7, f_7, 'b-.', 'linewidth', 2);
% hold on
% plot(x_8, f_8, 'm:', 'linewidth', 2);
% h = legend('VM:CPU-VM:CPU', 'VM:RAM-VM:RAM', 'VM:CPU-VM:RAM', 'VM:CPU-VM:RAM (Pair)');
% set(h, 'box','on','location','southeast','fontsize',15);
% set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1]);
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% % set(gca, 'xscale', 'log');
% ylabel('CDF', 'fontsize', 18); 
% xlabel('Median XCF', 'fontsize', 18);
% xlabel('Mean XCF', 'fontsize', 18);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'VM_VM_CPU_MEM_mean_XCF'));
% 
% % CDF of median
% fig = figure;
% set(fig,'Position',[200, 200, 600, 400]);
% set(gca, 'fontsize', 18);
% [f_5, x_5] = ecdf(vm_cpu_vm_cpu_all(:,2));
% [f_6, x_6] = ecdf(vm_mem_vm_mem_all(:,2));
% [f_7, x_7] = ecdf(vm_cpu_vm_mem_all(:,2));
% [f_8, x_8] = ecdf(vm_cpu_vm_mem_pair_all(:,2));
% plot(x_5, f_5, 'k-', 'linewidth', 2)
% hold on
% plot(x_6, f_6, 'r--', 'linewidth', 2)
% hold on
% plot(x_7, f_7, 'b-.', 'linewidth', 2);
% hold on
% plot(x_8, f_8, 'm:', 'linewidth', 2);
% h = legend('VM:CPU-VM:CPU', 'VM:RAM-VM:RAM', 'VM:CPU-VM:RAM', 'VM:CPU-VM:RAM (Pair)');
% set(h, 'box','on','location','southeast','fontsize',15);
% set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1]);
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% % set(gca, 'xscale', 'log');
% ylabel('CDF', 'fontsize', 18); 
% xlabel('Median XCF', 'fontsize', 18);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'VM_VM_CPU_MEM_median_XCF'));
% 
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
