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
max_lag = 1;
xcf_thres = 0.7;

percentile_to_check = 90;

mkdir('../New_Data_time_series_linear_subset_figure');
path = '../New_Data_time_series_linear_subset_figure/';

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

    for vm_id = 2 : size_box
        all_time_series = [all_time_series, ...
                    box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5)];
    end
    
    % Second, calculate the dissimilarities among these time series using
    % cross-correlation. Need to notice that, here we change the xcf into
    % the range of [0,1], which makes it more sense to measure the
    % dissimilarity

    total_metric_num = size(all_time_series, 2);
    d = ones(total_metric_num); d_half = zeros(total_metric_num);
    cross_all = ones(total_metric_num);
    for row = 2 : total_metric_num
        for col = 1 : row -1    
            [xcf, lags]= crosscorr(all_time_series(:, row), all_time_series(:, col), max_lag);
            d(row, col) = abs(xcf(2));
            d(col, row) = d(row, col);
            d_half(row, col) = d(row, col);                     
        end
    end    
           
    labels = {};
    for vm_id = 1 : size_box
        vm_label_no = {strcat(mat2str(vm_id))};
        labels = {labels{:}, vm_label_no{:}};
    end
    
    % Step 1: Do clustering and pick up the representative for each cluster 
    great_than_thres = d >= xcf_thres;
    overall_xcf = [];
    for vm_no = 1 : total_metric_num
        chosen_xcf = d(vm_no, great_than_thres(vm_no, :));
        overall_xcf(end+1, 1:2) = [sum(great_than_thres(vm_no,:)), mean(chosen_xcf)];
    end
    
    [overall_xcf, idx] = sortrows(overall_xcf, [-1,-2]);
    
    % Test
    original_overall_xcf = overall_xcf; original_idx = idx;
    
    xcf_cluster = {};
    representative_cluster = [];
    representative_cluster_time_series = [];
    vm_no = 1;
    while vm_no <= total_metric_num
        if overall_xcf(vm_no, 1) == 0 
            vm_no = vm_no + 1; continue;
        end
        
        xcf_cluster{end+1} = [idx(vm_no)];
        representative_cluster(end+1) = idx(vm_no);
        representative_cluster_time_series(:,end+1) = all_time_series(:, idx(vm_no));
        % exclude itself from the candidates in the cluster
        d(idx(vm_no),idx(vm_no)) = 0;
        cand_in_cluster = find(d(idx(vm_no),:) >= xcf_thres);
        xcf_cluster{end} = [xcf_cluster{end}, cand_in_cluster];
        
        % change all the related column to '0' in xcf matrix
        d(idx(vm_no), :) = 0; d(:, idx(vm_no)) = 0;
        
        % change all the chosen candidate as '0'
        for cand_id = 1 : numel(cand_in_cluster)
            idx_in_overall_xcf = find(idx == cand_in_cluster(cand_id));
            overall_xcf(idx_in_overall_xcf, 1:2) = 0;
            d(cand_in_cluster(cand_id), :) = 0; 
            d(:, cand_in_cluster(cand_id)) = 0;
        end
        
        vm_no = vm_no + 1;
    end
    
    % Step 2: Do stepwise choose, the representative in each cluster are
    % grouped in a big representative cluster. Then get the best subset of
    % the representative cluster
    final_representative_cluster = [];
    cand_id = numel(representative_cluster);
    time_len = numel(all_time_series(:,1));
    while cand_id >= 1 
        Y = representative_cluster_time_series(:, cand_id);
        X = [representative_cluster_time_series(:, 1: cand_id -1), ...
             representative_cluster_time_series(:, cand_id+1:end)];
        [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X, Y, 'display', 'off');
        
        if sum(inmodel) > 0
            representative_cluster_time_series(:, cand_id) = 0;
        else
            final_representative_cluster(end+1) = representative_cluster(cand_id);
        end    
        cand_id = cand_id - 1;
    end
    % final_representative_cluster = fliplr(final_representative_cluster);
    
    
    % Step 3: fitting all the time series using the representatives
    X = [ones(time_len, 1), all_time_series(:, final_representative_cluster)];
    abs_error = []; ape = []; r_square = [];
    for vm_no = 1 : total_metric_num
        % if the time series is in the representative
        if numel(find(final_representative_cluster == vm_no)) == 1
            continue;
        end
        
        Y = all_time_series(:, vm_no);
        [b,bint,r,rint,stats] = regress(Y,X);
        
        ssresid = sum(r.^2);
        sstotal = (length(Y) - 1) * 
        
        abs_error(end+1) = mean(abs(r));
        ape(end+1) = nanmean(abs(r) ./ Y);
        r_square(end+1) = stats(1); 
    end
    
    % write the summary of this box information after linear fitting
    BOX_ID(end+1) = pm_id; 
    ORIGINAL_VM_NUM(end+1) = total_metric_num; 
    REDUCED_VM_NO(end+1) = numel(final_representative_cluster);
    
    ALL_ABS_ERROR(:, end+1) = [mean(abs_error); median(abs_error); prctile(abs_error,percentile_to_check)]; 
    ALL_APE(:,end+1) = [mean(ape); median(ape); prctile(ape,percentile_to_check)]; 
    ALL_R_SQUARE(:,end+1) = [mean(r_square); median(r_square); prctile(r_square,percentile_to_check)];
end

SUMMARY_METRICS = [BOX_ID', ORIGINAL_VM_NUM', REDUCED_VM_NO', ALL_ABS_ERROR(1,:)', ALL_APE(1,:)', ALL_R_SQUARE(1,:)', ALL_R_SQUARE(3,:)'];
SUMMARY_METRICS = sortrows(SUMMARY_METRICS, [-2, -3]);
metric_name = {'BOX ID', 'Original # of VMs', 'Used # of VMs', 'Mean (mean) Abs Error', ...
               'Median (mean) Abs Error', '90%ile (mea) Abs Error', 'Mean (mean) APE', ...
               'Median (mean) APE', '90%ile (mea) APE', 'Mean R-square', ...
               'Median R-square', '90%ile R-square'};
dlmwrite('All_metrics.txt', metric_name);
dlmwrite('All_metrics.txt', SUMMARY_METRICS, '-append');


fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_ave, x_ave] = ecdf(ALL_APE(1,:) * 100);
[f_median, x_median] = ecdf(ALL_APE(2,:) * 100);
[f_90, x_90] = ecdf(ALL_APE(3,:) * 100);
plot(x_ave, f_ave, 'k-', 'linewidth', 2)
hold on
plot(x_median, f_median, 'r-.', 'linewidth', 2)
% hold on
% plot(x_90, f_90, 'b:', 'linewidth', 2)
h = legend('Mean', 'Median');
set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('Aboslute Percentage Error of CPU USED PCT (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_ape'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_ave, x_ave] = ecdf(ALL_R_SQUARE(1,:));
[f_median, x_median] = ecdf(ALL_R_SQUARE(2,:));
[f_90, x_90] = ecdf(ALL_R_SQUARE(3,:));
plot(x_ave, f_ave, 'k-', 'linewidth', 2)
hold on
plot(x_median, f_median, 'r-.', 'linewidth', 2)
hold on
plot(x_90, f_90, 'b:', 'linewidth', 2)
h = legend('Mean', 'Median','90%ile');
set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.2 : 1]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('Coefficient of Determination for Linear Fitting (R^2)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_R2'));


fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_origin, x_origin] = ecdf(ORIGINAL_VM_NUM);
[f_reduce, x_reduce] = ecdf(REDUCED_VM_NO);
plot(x_origin, f_origin, 'k-', 'linewidth', 2)
hold on
plot(x_reduce, f_reduce, 'r-.', 'linewidth', 2)
h = legend('Original # of VMs', 'Actual Used # of VMs');
set(h, 'box','on','location','southeast','fontsize',18);
set(gca,'xlim', [0 80]); set(gca, 'xtick', [0 : 20 : 80]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('Number of VMs w/ PRACTISE', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_reduce_vm'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_perc, x_perc] = ecdf((ORIGINAL_VM_NUM-REDUCED_VM_NO) ./ ORIGINAL_VM_NUM * 100);
plot(x_perc, f_perc, 'k-', 'linewidth', 2);
% h = legend('Original # of VMs', 'Reduced # of VMs');
% set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('Percentage of Reduced VMs (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_reduce_vm_percent'));
