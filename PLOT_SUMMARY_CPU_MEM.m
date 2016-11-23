% this script is used to generate what we are interesting metrics
clear
close all
clc

% Algorithm1: Use all VM Predictors
path1 = '../New_Data_time_series_linear_subset_vif_cpu_mem_figure/';
load (strcat(path1, 'BOX_APE_BASELINE'))
load (strcat(path1, 'BOX_R_SQUARE_BASELINE'))
load (strcat(path1, 'PRIO_REDUCED_TICKET_PCT'))
box_ape_1 = BOX_APE_ALL_VM; % store mean APE for each box
box_rsquare_1 = BOX_R_SQUARE_ALL_VM; % store R2 for each box
ticket_reduced_pct_1 = PRIO_REDUCED_PCT; % store ticket reduction for each box

% Algorithm2: COX + VIF + STEPWISE
load (strcat(path1, 'BOX_VM_APE_COX_VIF_STEPWISE'))
load (strcat(path1, 'BOX_VM_R_SQUARE_COX_VIF_STEPWISE'))
load (strcat(path1, 'VM_REDUCED_PCT_COX_VIF_STEPWISE'))
load (strcat(path1, 'REDUCED_TICKET_PCT'))
load (strcat(path1, 'ORIGINAL_VM_NUM'))
box_ape_2 = ALL_APE(1, :); % store mean APE for each box
vm_ape_2 = ALL_APE(2, :); % store mean (mean) APE of VMs for each box
box_rsquare_2 = ALL_R_SQUARE(1, :); % store R2 for each box
vm_rsquare_2 = ALL_R_SQUARE(2, :); % store mean (mean) R2 of VMs for each box
vm_reduced_pct_2 = VM_REDUCED_PCT; % store each box vm reduction pct
ticket_reduced_pct_2 = REDUCED_PCT; % store each box ticket reduction pct

% Algorithm3: DWT + VIF + STEPWISE
% path2 = '../New_Data_time_series_linear_subset_dwt_figure/';
% load (strcat(path2, 'BOX_VM_APE_DWT_STEPWISE'))
% load (strcat(path2, 'BOX_VM_R_SQUARE_DWT_STEPWISE'))
% load (strcat(path2, 'VM_REDUCED_PCT_DWT_STEPWISE'))
% load (strcat(path2, 'REDUCED_TICKET_PCT'))
% box_ape_3 = ALL_APE(1, :); % store mean APE for each box
% vm_ape_3 = ALL_APE(2, :); % store mean (mean) APE of VMs for each box
% box_rsquare_3 = ALL_R_SQUARE(1, :); % store R2 for each box
% vm_rsquare_3 = ALL_R_SQUARE(2, :); % store mean (mean) R2 of VMs for each box
% vm_reduced_pct_3 = VM_REDUCED_PCT; % store each box vm reduction pct
% ticket_reduced_pct_3 = REDUCED_PCT; % store each box ticket reduction pct

%%%%%%%%%%%%%%%%%%%%%% Plot the comparison results%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir('../New_Data_VM_Ticket_Reduction_Performance_CPU_MEM_Figure');
path = '../New_Data_VM_Ticket_Reduction_Performance_CPU_MEM_Figure/';

% Figure 1: BOX APE & R2 comparison
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box_1, x_box_1] = ecdf(box_ape_1*100);
[f_box_2, x_box_2] = ecdf(box_ape_2*100);
% [f_box_3, x_box_3] = ecdf(box_ape_3*100);
plot(x_box_1, f_box_1, 'k-', 'linewidth', 2)
hold on
plot(x_box_2, f_box_2, 'r--', 'linewidth', 2)
% hold on
% plot(x_box_3, f_box_3, 'b-.', 'linewidth', 2)
h = legend('Alg-1', 'Alg-2');
set(h, 'box','on','location','southeast','fontsize',18);
set(gca,'xlim', [0 50]); set(gca, 'xtick', [0 : 5 : 50]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
% title('Fitting Error');
ylabel('CDF', 'fontsize', 18); 
xlabel('BOX: Aboslute Percentage Error of CPU USED PCT (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_box_ape'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box_1, x_box_1] = ecdf(box_rsquare_1);
[f_box_2, x_box_2] = ecdf(box_rsquare_2);
% [f_box_3, x_box_3] = ecdf(box_rsquare_3);
plot(x_box_1, 1-f_box_1, 'k-', 'linewidth', 2)
hold on
plot(x_box_2, 1-f_box_2, 'r--', 'linewidth', 2)
% hold on
% plot(x_box_3, 1-f_box_3, 'b-.', 'linewidth', 2)
h = legend('Alg-1', 'Alg-2');
set(h, 'box','on','location','southwest','fontsize',18);
set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
% title('Fitting Error');
ylabel('CCDF', 'fontsize', 18); 
xlabel('BOX: Coeffecient of Determination (R^2)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_box_r2'));

% Figure 2: VM ape and R2
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box_2, x_box_2] = ecdf(vm_ape_2*100);
% [f_box_3, x_box_3] = ecdf(vm_ape_3*100);
plot(x_box_2, f_box_2, 'r--', 'linewidth', 2)
% hold on
% plot(x_box_3, f_box_3, 'b-.', 'linewidth', 2)
h = legend('Alg-2');
set(h, 'box','on','location','southeast','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 10 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
% title('Fitting Error');
ylabel('CDF', 'fontsize', 18); 
xlabel('VM: Mean Aboslute Percentage Error of CPU USED PCT (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_vm_ape'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box_2, x_box_2] = ecdf(vm_rsquare_2);
[f_box_3, x_box_3] = ecdf(vm_rsquare_3);
plot(x_box_2, 1-f_box_2, 'r--', 'linewidth', 2)
hold on
plot(x_box_3, 1-f_box_3, 'b-.', 'linewidth', 2)
h = legend('Alg-2', 'Alg-3');
set(h, 'box','on','location','southwest','fontsize',18);
set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
% title('Fitting Error');
ylabel('CCDF', 'fontsize', 18); 
xlabel('VM: Mean Coeffecient of Determination (R^2)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_vm_r2'));

% Figure 3: vm reduction
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box_2, x_box_2] = ecdf(vm_reduced_pct_2 * 100);
[f_box_3, x_box_3] = ecdf(vm_reduced_pct_3 * 100);
plot(x_box_2, f_box_2, 'r--', 'linewidth', 2)
hold on
plot(x_box_3, f_box_3, 'b-.', 'linewidth', 2)
h = legend('Alg-2', 'Alg-3');
set(h, 'box','on','location','southeast','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 10 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
% title('Fitting Error');
ylabel('CDF', 'fontsize', 18); 
xlabel('VM Reduction Percentage Per BOX (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_vm_reduction'));

% let's only focus on the big box
big_idx = ORIGINAL_VM_NUM >= 40;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box_2, x_box_2] = ecdf(vm_reduced_pct_2(big_idx) * 100);
[f_box_3, x_box_3] = ecdf(vm_reduced_pct_3(big_idx) * 100);
plot(x_box_2, f_box_2, 'r--', 'linewidth', 2)
hold on
plot(x_box_3, f_box_3, 'b-.', 'linewidth', 2)
h = legend('Alg-2', 'Alg-3');
set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 10 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
% title('Fitting Error');
ylabel('CDF', 'fontsize', 18); 
xlabel('VM Reduction Percentage Per Big BOX (>= 40 VMs) (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_big_box_vm_reduction'));

% Figure 4: Ticket Reduction
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box_1, x_box_1] = ecdf(ticket_reduced_pct_1 * 100);
[f_box_2, x_box_2] = ecdf(ticket_reduced_pct_2 * 100);
[f_box_3, x_box_3] = ecdf(ticket_reduced_pct_3 * 100);
plot(x_box_1, 1-f_box_1, 'k-', 'linewidth', 2)
hold on
plot(x_box_2, 1-f_box_2, 'r--', 'linewidth', 2)
hold on
plot(x_box_3, 1-f_box_3, 'b-.', 'linewidth', 2)
h = legend('Alg-1','Alg-2', 'Alg-3');
set(h, 'box','on','location','southwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 10 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
% title('Fitting Error');
ylabel('CCDF', 'fontsize', 18); 
xlabel('Ticket Reduction Percentage Per BOX (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_ticket_reduction'));


%%%%%%%%%%%%%%%%%%%%%%%%%%Display the table%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
per_to_check = 90;
% Algorithm1:
box_ape_summary_1 = [mean(box_ape_1), median(box_ape_1), prctile(box_ape_1, per_to_check)];
vm_ape_summary_1 = [0, 0, 0];
vm_reduction_summary_1 = [0, 0, 0];
ticket_reduction_summary_1 = [mean(ticket_reduced_pct_1), median(ticket_reduced_pct_1), ...
                              prctile(ticket_reduced_pct_1, per_to_check)];
                          
% Algorithm2:
box_ape_summary_2 = [mean(box_ape_2), median(box_ape_2), prctile(box_ape_2, per_to_check)];
vm_ape_summary_2 = [nanmean(vm_ape_2), nanmedian(vm_ape_2), prctile(vm_ape_2, per_to_check)];
vm_reduction_summary_2 = [nanmean(vm_reduced_pct_2(big_idx)), nanmedian(vm_reduced_pct_2(big_idx)), ...
                              prctile(vm_reduced_pct_2(big_idx), per_to_check)];
ticket_reduction_summary_2 = [nanmean(ticket_reduced_pct_2), nanmedian(ticket_reduced_pct_2), ...
                              prctile(ticket_reduced_pct_2, per_to_check)];
                          
% Algorithm3:
box_ape_summary_3 = [mean(box_ape_3), median(box_ape_3), prctile(box_ape_3, per_to_check)];
vm_ape_summary_3 = [nanmean(vm_ape_3), nanmedian(vm_ape_3), prctile(vm_ape_3, per_to_check)];
vm_reduction_summary_3 = [nanmean(vm_reduced_pct_3(big_idx)), nanmedian(vm_reduced_pct_3(big_idx)), ...
                              prctile(vm_reduced_pct_3(big_idx), per_to_check)];
ticket_reduction_summary_3 = [nanmean(ticket_reduced_pct_3), nanmedian(ticket_reduced_pct_3), ...
                              prctile(ticket_reduced_pct_3, per_to_check)];
