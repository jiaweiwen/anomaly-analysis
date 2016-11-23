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
var_thres = 95;

percentile_to_check = 90;

mkdir('../New_Data_time_series_linear_subset_pca_figure');
path = '../New_Data_time_series_linear_subset_pca_figure/';

BOX_ID = []; ORIGINAL_VM_NUM = []; REDUCED_VM_NO = [];
ALL_ABS_ERROR = []; ALL_APE = []; ALL_R_SQUARE = [];

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
%     
    if numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1)) <= 70
        continue;
    end
    
    % First, extract all of the time series
    time = box_vm_time_series_summary{1, box_id}{1, 1}(:,1) - box_vm_time_series_summary{1, box_id}{1, 1}(1,1);
    pm_id = box_vm_time_series_summary{1, box_id}{1, 1}(1,2);
    all_time_series = []; %box_vm_time_series_summary{1, box_id}{1, 1}(:,4);
    box_time_series = box_vm_time_series_summary{1, box_id}{1, 1}(:,4);

    pca_used_time_series_idx = [];
    for vm_id = 2 : size_box
        if var(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5)) ~= 0
            pca_used_time_series_idx(end+1) = vm_id -1;
        end
        all_time_series = [all_time_series, ...
                    box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5)];
    end
    
    % Second, applying PCA to get the principal component scores, which
    % will help to reduce the dimension of all VMs.
    [wcoeff, score, latent, ~, explained] = pca(all_time_series(:,pca_used_time_series_idx), ...
                                            'VariableWeights', 'variance');
    sum_explained = cumsum(explained);
    [~, num_components] = max(sum_explained >= var_thres);
   
    time_len = numel(score(:,1));
    
    % Thirdly, fitting all the time series using the representatives
    X = [ones(time_len, 1), score(:, 1:num_components)];
    abs_error = []; ape = []; r_square = [];
    for vm_no = 1 : size_box -1
        Y = all_time_series(:, vm_no);
        [b,bint,r,rint,stats] = regress(Y,X);
        
        abs_error(end+1) = mean(abs(r));
        ape(end+1) = nanmean(abs(r) ./ Y);
        r_square(end+1) = stats(1); 
    end
    
    % Step 5: fit box time series using the representative
    [b, bint, r, rint, stats] = regress(box_time_series, X);
    box_abs_error = mean(abs(r));
    box_ape = nanmean(abs(r) ./ box_time_series);
    box_r_square = stats(1);
    
    % write the summary of this box information after linear fitting
    BOX_ID(end+1) = pm_id; 
    ORIGINAL_VM_NUM(end+1) = size_box - 1; 
    REDUCED_VM_NO(end+1) = num_components;
    
    ALL_ABS_ERROR(:, end+1) = [box_abs_error; mean(abs_error); median(abs_error); prctile(abs_error,percentile_to_check)]; 
    ALL_APE(:,end+1) = [box_ape; mean(ape); median(ape); prctile(ape,percentile_to_check)]; 
    ALL_R_SQUARE(:,end+1) = [box_r_square; mean(r_square); median(r_square); prctile(r_square,percentile_to_check)];
end

SUMMARY_METRICS = [BOX_ID', ORIGINAL_VM_NUM', REDUCED_VM_NO', ALL_ABS_ERROR(1:2,:)', ALL_APE(1:2,:)', ALL_R_SQUARE(1:2,:)', ALL_R_SQUARE(3,:)'];
SUMMARY_METRICS = sortrows(SUMMARY_METRICS, [-2, -3]);
% metric_name = {'BOX ID', 'Original # of VMs', 'Used # of VMs', 'Mean (mean) Abs Error', ...
%                'Median (mean) Abs Error', '90%ile (mea) Abs Error', 'Mean (mean) APE', ...
%                'Median (mean) APE', '90%ile (mea) APE', 'Mean R-square', ...
%                'Median R-square', '90%ile R-square'};
% dlmwrite('All_metrics.txt', metric_name);
% dlmwrite('All_metrics.txt', SUMMARY_METRICS, '-append');

% vm: ape cdf
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box, x_box] = ecdf(ALL_APE(1,:)*100);
[f_ave, x_ave] = ecdf(ALL_APE(2,:) * 100);
[f_median, x_median] = ecdf(ALL_APE(3,:) * 100);
[f_90, x_90] = ecdf(ALL_APE(4,:) * 100);
plot(x_ave, f_ave, 'k-', 'linewidth', 2)
% hold on
% plot(x_median, f_median, 'r-.', 'linewidth', 2)
% hold on
% plot(x_90, f_90, 'b:', 'linewidth', 2)
% h = legend('Mean', 'Median');
% set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
title('VM Fitting');
ylabel('CDF', 'fontsize', 18); 
xlabel('Aboslute Percentage Error of CPU USED PCT (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_vm_ape'));

% box: ape cdf
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
plot(x_box, f_box, 'k-', 'linewidth', 2)
% hold on
% plot(x_90, f_90, 'b:', 'linewidth', 2)
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
title('BOX Fitting');
ylabel('CDF', 'fontsize', 18); 
xlabel('Aboslute Percentage Error of CPU USED PCT (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_box_ape'));

% vm: r-square cdf
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box, x_box] = ecdf(ALL_R_SQUARE(1,:));
[f_ave, x_ave] = ecdf(ALL_R_SQUARE(2,:));
[f_median, x_median] = ecdf(ALL_R_SQUARE(3,:));
[f_90, x_90] = ecdf(ALL_R_SQUARE(4,:));
plot(x_ave, f_ave, 'k-', 'linewidth', 2)
hold on
plot(x_median, f_median, 'r-.', 'linewidth', 2)
hold on
plot(x_90, f_90, 'b:', 'linewidth', 2)
h = legend('Mean', 'Median','90%ile');
set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.2 : 1]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
title('VM Fitting');
ylabel('CDF', 'fontsize', 18); 
xlabel('Coefficient of Determination for Linear Fitting (R^2)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_vm_R2'));

% box: r-square cdf
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
plot(x_box, f_box, 'b:', 'linewidth', 2)
set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.2 : 1]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
title('BOX Fitting');
ylabel('CDF', 'fontsize', 18); 
xlabel('Coefficient of Determination for Linear Fitting (R^2)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_box_R2'));

% Check the percentage of reduction
reduced_pct = (ORIGINAL_VM_NUM-REDUCED_VM_NO) ./ ORIGINAL_VM_NUM * 100;
box_size = floor(ORIGINAL_VM_NUM / 5) * 5;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(reduced_pct', box_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
%set(gca, 'xlim', [0 120]); set(gca, 'xtick', [0 : 20 : 120], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('Reduced Percentage of Used VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_vm_percent'));

% % Need to notice that: Let's only look at the huge vm numbers BOX
% big_box_idx = ORIGINAL_VM_NUM >= 20;
% fig = figure;
% set(fig,'Position',[200, 200, 600, 400]);
% set(gca, 'fontsize', 18);
% reduced_pct = (ORIGINAL_VM_NUM(big_box_idx)-REDUCED_VM_NO(big_box_idx)) ...
%                ./ ORIGINAL_VM_NUM(big_box_idx) * 100;
% [f_perc, x_perc] = ecdf(reduced_pct);
% plot(x_perc, f_perc, 'k-', 'linewidth', 2);
% % h = legend('Original # of VMs', 'Reduced # of VMs');
% % set(h, 'box','on','location','northwest','fontsize',18);
% set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% % set(gca, 'xscale', 'log');
% ylabel('CDF', 'fontsize', 18); 
% xlabel('Percentage of Reduced VMs (%)', 'fontsize', 18);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'cdf_reduce_vm_percent_big_box'));

