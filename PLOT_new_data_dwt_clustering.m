% this script is to plot the figure of DWT clustering
close all; clear; clc

load ../New_Data/all_cluster_dwt
load ../New_Data/all_box_size_dwt

% Just change into the number of clusters
all_box_size = all_box_size / 2;
path = '../New_Data_time_series_cluster_dwt_figure/';

% we have cut the cluster number, due to > 32 only counts for less than 1
% percent
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
bins = [[2,3,4,5,6], [8, 16, 32]];
xlabel_name = {'2'; '2-3';'3-4'; '4-5';'5-6';'6-8';'8-16';'16-32'};
[N, X] = hist(all_cluster, bins);
X = 1 : numel(X);
bar(X, N/sum(N) * 100);
set(gca, 'xlim', [0 numel(X)+1]); set(gca, 'xticklabel',xlabel_name , 'fontsize', 18);
set(gca, 'ylim', [0 60]); set(gca, 'ytick', [0 : 10: 60]);
ylabel('Percentage (%) of Boxes'); xlabel('Number of Clusters');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'dwt_number_cluster_pdf'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
ecdf(all_box_size);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
set(gca, 'xscale', 'log');
ylabel('CDF'); xlabel('Number of VMs Per Box');
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'dwt_box_size_cdf'));

% remove the tails of box size
upper_bound = prctile(all_box_size, 97);
used_box_size_idx = find(all_box_size < upper_bound);

cand_all_box_size = ceil(log2(all_box_size(used_box_size_idx)));
cand_all_cluster = all_cluster(used_box_size_idx);

max_box_size = max(cand_all_box_size);
max_power_idx = ceil(log2(max_box_size));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
boxplot(cand_all_cluster', cand_all_box_size');
set(gca, 'ylim', [0 40]); set(gca, 'ytick', [0,2,4,8,16,32,40], 'fontsize', 18);
set(gca, 'xticklabel',{'2','2-4','4-8','8-16','16-32','32-64'});
% set(gca, 'yscale', 'log');
ylabel('Number of Clusters'); xlabel('Number of VMs Per Box');
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'dwt_number_cluster_box_size'));

