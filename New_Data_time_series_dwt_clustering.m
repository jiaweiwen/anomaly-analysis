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
max_lag = 4;

mkdir('../New_Data_time_series_cluster_dwt_figure');
path = '../New_Data_time_series_cluster_dwt_figure/';

all_coeff_determine = [];
all_error = [];

all_cluster = []; all_box_size = []; all_used_pct = [];
for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary{1, box_id});
    
    % If we don't have time series
    if size_box <= 2
        continue;
    end
%     
%     if box_id > 50
%         break;
%     end
    
%     if numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1)) <= 70
%         continue;
%     end
    
    % First, extract all of the time series
    time = box_vm_time_series_summary{1, box_id}{1, 1}(:,1) - box_vm_time_series_summary{1, box_id}{1, 1}(1,1);
    pm_id = box_vm_time_series_summary{1, box_id}{1, 1}(1,2);
    all_time_series = []; %box_vm_time_series_summary{1, box_id}{1, 1}(:,4);

    for vm_id = 2 : size_box
        all_time_series = [all_time_series, ...
                    box_vm_time_series_summary{1, box_id}{1, vm_id}(:, end), ...
                    box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, end)];
    end
    
    % Second, calculate the dissimilarities among these time series using
    % cross-correlation. Need to notice that, here we change the xcf into
    % the range of [0,1], which makes it more sense to measure the
    % dissimilarity

    total_metric_num = size(all_time_series, 2);
    d = zeros(total_metric_num); d_half = zeros(total_metric_num);
    cross_all = ones(total_metric_num);
    for row = 2 : total_metric_num
        for col = 1 : row -1
            d(row, col) = dtw(all_time_series(:, row), ...
                              all_time_series(:, col), max_lag);
            d(col, row) = d(row, col);
            d_half(row, col) = d(row, col);
            % in case of totally same case
            if d_half(row, col) == 0
                d_half(row, col) = -1;
            end
        end
    end    
    
       
    % in case of too many clustering
%     if total_metric_num < 20
%         continue;
%     end
%                   
    labels = {};
    for vm_id = 1 : size_box
        vm_label_no = strcat('VM', strcat(mat2str(vm_id)));
        labels = {labels{:}, vm_label_no};
    end
%     
%     % Plot the heat map of distance matrix
%     fig = figure;
%     set(fig,'Position',[200,200,1200,1000]);
%     heatmap(cross_all, labels, labels, '%0.2f', 'TextColor', 'w', ...
%             'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
%             'TickAngle', 45, 'TickFontSize',6, 'FontSize', 6);
%     caxis([-1 1]);
%     title_name = strcat('BOX ID=', mat2str(pm_id));
%     title(title_name);
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path1, 'PM_', mat2str(pm_id), '_XCF'));
% 
    % Third, build the hierarchical tree for clustering
    d_half_vec = reshape(d_half, 1, total_metric_num^2);
    d_vec = nonzeros(d_half_vec)';
    d_vec(find(d_vec == -1)) = 0;
    z = linkage(d_vec, 'single'); % hierarchical clustering
    
    % Verify dissimilarity
%     c = cophenet(z, d_vec);
%     disp(strcat('The accuracy of clustering ([0,1]) is ', mat2str(c)));
    
     % Create clusters and pick up the optimal number of clusters
    s_diff_cluster = [];
    for num_cluster = 2 : max(2, ceil(total_metric_num/2))
        cidx = cluster(z, 'maxclust', num_cluster);
        s_diff_cluster(num_cluster-1) = mean(silhouette([], cidx, d_vec));
    end
    
%     disp('The silhuette results are as follows');
%     disp(s_diff_cluster);
    [max_s, idx] = max(s_diff_cluster);

    all_cluster(end+1) = idx + 1;
    all_box_size(end+1) = total_metric_num;
    all_used_pct(end+1) = (idx+1) / total_metric_num;
    
    %%%%%% attention: when we determine to do the optimal clustering %%%%%%
    % Do the optimal clustering 
    % cidx = cluster(z, 'maxclust', idx + 1);   
       
    % Display dendrogram
%     fig = figure; 
%     set(fig,'Position',[200, 200, 600, 400]);
%     set(gca,'fontsize', 10);
%     [H, T] = dendrogram(z);
%     hAxis = get(H(1),'parent')
%     % Get the permutation of the nodes
%     perm = str2num(get(hAxis,'XtickLabel'));
%     title(strcat('PM = ', mat2str(pm_id), ', cophenetic correlation coefficient = ', mat2str(floor(c*1000)/1000)));
%     set(gca,'xticklabel',labels(perm),'fontsize', 10);
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, 'dendrogram_pm_', mat2str(pm_id)));
    
%     disp('The clustering results are as follows');
%     disp(z);
   
    % Plot these time series (each cluster in one figure)
    
    if pm_id == 1021577
        cc = hsv(total_metric_num/2);
        line_style = {'-', '--', '-.', ':'};
        fig = figure;
        set(fig,'Position',[200, 200, 900, 400]);
        cand_legend = {};
        for time_series_idx = 1 : total_metric_num/2
            cand_legend = {cand_legend{:}, labels{time_series_idx}};
            plot(time/3600, all_time_series(:, time_series_idx), ...
                        'color', cc(time_series_idx, :), 'linewidth',2, ...
                        'linestyle', line_style{time_series_idx});
            hold on
        end
        h = legend(cand_legend);
        set(h,'location','northwest','box','on','fontsize',18)
        set(gca, 'xlim', [time(1)/3600 24]);
        set(gca, 'xtick', [time(1)/3600 : 2 : 24],'fontsize',18);
        xlabel('Time (hour)');
        ylabel('CPU USED PCT (%)');
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'PM_', mat2str(pm_id)));  
        break;
    end
    
%     cc = hsv(total_metric_num);
%     fig = figure;
%     set(fig,'Position',[200, 200, 1500, 1500]);
%     for fig_idx = 1 : idx + 1
%         subplot(idx+1, 2, (fig_idx-1)*2 + 1);
%         cand_legend = {};
%         time_series_set = {};
%         for time_series_idx = 1 : total_metric_num
%             if cidx(time_series_idx) == fig_idx
%                 cand_legend = {cand_legend{:}, labels{time_series_idx}};
%                 plot(time/3600, all_time_series(:, time_series_idx), ...
%                     'color', cc(time_series_idx, :));
%                 time_series_set{end + 1} = all_time_series(:, time_series_idx);
%                 hold on                
%             end     
%         end
%         h = legend(cand_legend);
%         set(h,'location','eastoutside','box','off','fontsize',10)
%         if fig_idx == 1
%             title(strcat('PM = ', mat2str(pm_id), ', Cluster Index = ', mat2str(fig_idx)));
%         else
%             title(strcat('Cluster Index = ', mat2str(fig_idx)));
%         end
%         set(gca, 'xlim', [time(1)/3600 ceil(time(end)/3600)+4]);
%         set(gca, 'xtick', [time(1)/3600 : 2 : ceil(time(end)/3600)+4],'fontsize',10);
%         xlabel('Time (hour)');
%         ylabel('CPU USED PCT (%)');
%         
%         if numel(time_series_set) == 1
%             continue;
%         end
%         
%         % Pick up the representatives of the time series set
%         [error_summary, time_series_rank, final_error_summary, ...
%          final_coeff_determine, final_new_time_series_set] = PickRepresent(time_series_set, max_lag);
%         all_error(end+1) = final_error_summary(1);
%         all_coeff_determine(end+1) = final_coeff_determine(1);
%         subplot(idx+1, 2, fig_idx * 2);
%         x_tick = 1 : numel(time_series_set);
%         plot(x_tick, error_summary(:,1)*100, 'r*-');
%         hold on
%         plot(x_tick, error_summary(:,2)*100, 'bo-');
%         hold on
%         plot(x_tick(final_new_time_series_set), error_summary(final_new_time_series_set,1)*100, 'k^', 'markersize', 10)
%         hold on
%         plot(x_tick(final_new_time_series_set), error_summary(final_new_time_series_set,2)*100, 'k^', 'markersize', 10)
%         h = legend('Mean APE','Worst APE','Picked Points');
%         set(h,'location','eastoutside','box','on','fontsize',10)
%         set(gca,'xticklabel',{time_series_rank},'fontsize', 10);
%         xlabel('VM ID');
%         ylabel('Abosulte Percentage Error (%)');
%         max_ape = max(error_summary(:,2)) * 100;
%         set(gca, 'xtick', x_tick);
%         set(gca, 'ylim', [0 max_ape+1]);
%     end   
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path, 'PM_', mat2str(pm_id),'_Cluster'));
%     close all
end

save('../New_Data/all_cluster_dwt', 'all_cluster');
save('../New_Data/all_box_size_dwt', 'all_box_size');
save('../New_Data/all_used_pct', 'all_used_pct');


fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
bins = 10 : 10 : 80;
[N, X] = hist(all_cluster, bins);
bar(X, N/sum(N) * 100);
set(gca, 'xlim', [0 85]); set(gca, 'xtick', [0 : 10 : 80], 'fontsize', 18);
set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0 : 10: 100]);
ylabel('Percentage (%) of Boxes'); xlabel('Number of Clusters');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'dwt_number_cluster_pdf'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
bins = 24 : 24 : 240;
[N, X] = hist(all_box_size, bins);
bar(X, N/sum(N) * 100);
set(gca, 'xlim', [0 252]); set(gca, 'xtick', [0 : 24 : 240], 'fontsize', 18);
set(gca, 'ylim', [0 70]); set(gca, 'ytick', [0 : 10: 70]);
ylabel('Percentage (%) of Boxes'); xlabel('Box Size');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'dwt_box_size_pdf'));

all_box_size = ceil(all_box_size / 24) * 24;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
boxplot(all_cluster', all_box_size');
set(gca, 'ylim', [0 80]); set(gca, 'ytick', [0 : 10 : 80], 'fontsize', 18);
% set(gca, 'yscale', 'log');
ylabel('Number of Clusters'); xlabel('Box Size');
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'dwt_number_cluster_box_size'));


% fig = figure;
% set(fig,'Position',[200, 200, 600, 400]);
% set(gca, 'fontsize', 18);
% [f_ave, x_ave] = ecdf(all_error * 100);
% plot(x_ave, f_ave, 'k', 'linewidth', 2)
% % h = legend('Mean');
% % set(h, 'box','on','location','northwest','fontsize',18);
% set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% % set(gca, 'xscale', 'log');
% ylabel('CDF', 'fontsize', 18); 
% xlabel('Aboslute Percentage Error of CPU USED PCT (%)', 'fontsize', 18);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'cdf_ape'));
% 
% 
% fig = figure;
% set(fig,'Position',[200, 200, 600, 400]);
% set(gca, 'fontsize', 18);
% [f_ave, x_ave] = ecdf(all_coeff_determine);
% plot(x_ave, f_ave, 'k', 'linewidth', 2)
% % h = legend('Mean');
% % set(h, 'box','on','location','northwest','fontsize',18);
% set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1]);
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% % set(gca, 'xscale', 'log');
% ylabel('CDF', 'fontsize', 18); 
% xlabel('Coefficient of Determination for Linear Fitting (R^2)', 'fontsize', 18);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'cdf_R2'));
