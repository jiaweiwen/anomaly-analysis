% This script is used to do time series clustering
% 1. Apply DTW (dynamic time warping) to measure the dissimalrity between
%    different time serier.
% 2. Use Hierachical Tree to do clustering based on DTW
% 3. Apply silhuette to determine the optimal number of clusters.

close all; clear; clc

load ../box_vm_time_series_summary

size_box_vm = size(box_vm_time_series_summary);

% Determine the maximum time length
max_time = 96 * 3;
grat_small = 900;
max_lag = 2;

mkdir('../time_series_cluster_xcf_figure');
path = '../time_series_cluster_xcf_figure/';

mkdir('../XCF_for_clustering_figure');
path1 = '../XCF_for_clustering_figure/';

for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary{1, box_id});
    
    % If we don't have time series
    if size_box == 0
        continue;
    end
    
    % First, extract all of the time series
    time = box_vm_time_series_summary{1, box_id}{1, 1}(1:max_time,2);
    pm_id = box_vm_time_series_summary{1, box_id}{1, 1}(1,1);
    all_time_series = box_vm_time_series_summary{1, box_id}{1, 1}(1:max_time,3:4);
    for vm_id = 2 : size_box
        all_time_series = [all_time_series, ...
                  box_vm_time_series_summary{1, box_id}{1, vm_id}(1:max_time,4:5)];
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
            % assume there is no shift among the time series
            temp = crosscorr(all_time_series(:,row), all_time_series(:, col), max_lag);
            d(row, col) = 1 - max(abs(temp));
            cross_all(row, col) = max(temp);
            cross_all(col, row) = max(temp);
            d(col, row) = d(row, col);
            d_half(row, col) = d(row, col);
        end
    end 
           
    labels = {'BOX-CPU','BOX-MEM'};
    for vm_id = 1 : size_box
        vm_label_no = {strcat('VM',mat2str(vm_id),'-CPU'), ...
                       strcat('VM',mat2str(vm_id),'-MEM')};
        labels = {labels{:}, vm_label_no{:}};
    end
    
    % Plot the heat map of distance matrix
    fig = figure;
    set(fig,'Position',[200,200,500,300]);
    heatmap(cross_all, labels, labels, '%0.2f', 'TextColor', 'w', ...
            'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
            'TickAngle', 45, 'TickFontSize',10, 'FontSize', 10);
    caxis([-1 1]);
    title_name = strcat('BOX ID=', mat2str(pm_id));
    title(title_name);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path1, 'PM_', mat2str(pm_id), '_XCF'));

    % Third, build the hierarchical tree for clustering
    d_half_vec = reshape(d_half, 1, total_metric_num^2);
    d_vec = nonzeros(d_half_vec)';
    z = linkage(d_vec, 'single'); % hierarchical clustering
    
    % Verify dissimilarity
    c = cophenet(z, d_vec);
    disp(strcat('The accuracy of clustering ([0,1]) is ', mat2str(c)));
    
    % Display dendrogram
    fig = figure; 
    set(fig,'Position',[200, 200, 500, 300]);
    set(gca,'fontsize', 15);
    dendrogram(z);
    % set(gca,'yscale','log');
    title(strcat('PM = ', mat2str(pm_id), ', cophenetic correlation coefficient = ', mat2str(floor(c*1000)/1000)));
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'dendrogram_pm_', mat2str(pm_id)));
    
    % Create clusters and pick up the optimal number of clusters
    s_diff_cluster = [];
    for num_cluster = 2 : min(4, total_metric_num/2)
        cidx = cluster(z, 'maxclust', num_cluster);
        s_diff_cluster(num_cluster-1) = mean(silhouette([], cidx, d_vec));
    end
    
    disp('The silhuette results are as follows');
    disp(s_diff_cluster);
    [max_s, idx] = max(s_diff_cluster);
    
    % Do the optimal clustering 
    cidx = cluster(z, 'maxclust', idx + 1);
    
%     disp('The clustering results are as follows');
%     disp(z);
   
    % Plot these time series (each cluster in one figure)
    
    cc = hsv(total_metric_num);
    fig = figure;
    set(fig,'Position',[200, 200, 2000, 2000]);
    set(gca,'fontsize', 15);
    for fig_idx = 1 : idx + 1
        subplot(idx+1, 1, fig_idx);
        cand_legend = {};
        for time_series_idx = 1 : total_metric_num
            if cidx(time_series_idx) == fig_idx
                cand_legend = {cand_legend{:}, labels{time_series_idx}};
                plot(time/3600, all_time_series(:, time_series_idx), ...
                    'color', cc(time_series_idx, :));
                hold on                
            end     
        end
        h = legend(cand_legend);
        set(h,'location','northeast','box','off')
        if fig_idx == 1
            title(strcat('PM = ', mat2str(pm_id), ', Cluster Index = ', mat2str(fig_idx)));
        else
            title(strcat('Cluster Index = ', mat2str(fig_idx)));
        end
        set(gca, 'xlim', [time(1)/3600 time(end)/3600]);
        set(gca, 'xtick', [time(1)/3600 : max_time*grat_small/3600/24 : time(end)/3600]);
    end   
    xlabel('Time (hour)');
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'PM_', mat2str(pm_id),'_Cluster'));
    close all
end