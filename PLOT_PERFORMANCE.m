box_size = floor(ORIGINAL_VM_NUM / grat + 1) * grat;

% box+vm: ape cdf
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box, x_box] = ecdf(ALL_APE(1,:)*100);
[f_ave, x_ave] = ecdf(ALL_APE(2,:) * 100);
plot(x_box, f_box, 'k-', 'linewidth', 2)
hold on
plot(x_ave, f_ave, 'r-.', 'linewidth', 2)
h = legend('BOX', 'VM');
set(h, 'box','on','location','southeast','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
title('Fitting Error');
ylabel('CDF', 'fontsize', 18); 
xlabel('Aboslute Percentage Error of CPU USED PCT (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'cdf_vm_box_ape'));


% box+vm: r-square cdf
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_box, x_box] = ecdf(ALL_R_SQUARE(1,:));
[f_ave, x_ave] = ecdf(ALL_R_SQUARE(2,:));
plot(x_box, 1-f_box, 'k-', 'linewidth', 2)
hold on
plot(x_ave, 1-f_ave, 'r-.', 'linewidth', 2)
h = legend('BOX', 'VM');
set(h, 'box','on','location','southwest','fontsize',18);
set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.2 : 1]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
title('Fitting Goodness');
ylabel('CCDF', 'fontsize', 18); 
xlabel('Coefficient of Determination for Linear Fitting (R^2)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'ccdf_vm_box_R2'));

% original VM and used VM CDF 
vm_reduced_pct = (ORIGINAL_VM_NUM-REDUCED_VM_NO) ./ ORIGINAL_VM_NUM * 100;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(vm_reduced_pct', box_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
%set(gca, 'xlim', [0 120]); set(gca, 'xtick', [0 : 20 : 120], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('Reduced Percentage of Used VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_vm_percent'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(vm_reduced_pct', box_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
%set(gca, 'xlim', [0 120]); set(gca, 'xtick', [0 : 20 : 120], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('Reduced Percentage of Used VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_vm_percent'));

%%%%%%%%%%%%%%%%%%%%%%%look into ticket reduction%%%%%%%%%%%%%%%%%%%%%%%%%%
ticket_big_idx = ORIGINAL_TICKET > 0;
ORIGINAL_TICKET = ORIGINAL_TICKET(ticket_big_idx);
PRIO_RESIZE_TICKET = PRIO_RESIZE_TICKET(ticket_big_idx);
RESIZE_TICKET = RESIZE_TICKET(ticket_big_idx);
ORIGINAL_VM_NUM = ORIGINAL_VM_NUM(ticket_big_idx);
box_size = floor(ORIGINAL_VM_NUM / grat + 1) * grat;

% Ticket reduction
prio_reduced_pct = (ORIGINAL_TICKET - PRIO_RESIZE_TICKET) ./ ORIGINAL_TICKET * 100;
reduced_pct = (ORIGINAL_TICKET - RESIZE_TICKET) ./ ORIGINAL_TICKET * 100;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_only, x_only] = ecdf(prio_reduced_pct);
[f_both, x_both] = ecdf(reduced_pct);
plot(x_only, 1-f_only, 'k-', 'linewidth', 2)
% h = legend('Ticket Reduction Only', 'VM Reduction + Ticket Reduction');
% set(h, 'box','on','location','southwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
title('Ticket Reduction Only');
ylabel('CCDF', 'fontsize', 18); 
xlabel('Ticket Reduction Compared w/ Original Number of Tickets (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'ccdf_ticket_reduction_only'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
plot(x_both, 1-f_both, 'k-', 'linewidth', 2)
% h = legend('Ticket Reduction Only', 'VM Reduction + Ticket Reduction');
% set(h, 'box','on','location','southwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
title('VM Reduction + Ticket Reduction');
ylabel('CCDF', 'fontsize', 18); 
xlabel('Ticket Reduction Compared w/ Original Number of Tickets (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'ccdf_ticket_reduction_both'));

% Ticket 3D PLOT (VM+TICKET) V.S. TICKET ONLY
box_size = floor(ORIGINAL_VM_NUM / grat + 1) * grat;
ticket_size = floor(ORIGINAL_TICKET / grat + 1) * grat;
unique_box_size = unique(box_size');
max_box_size = max(box_size);
unique_ticket_num = unique(ticket_size');
max_ticket_size = max(ticket_size);

X = [];
for i = 1 : max_box_size
    X(i, 1 : max_ticket_size) = i;
end

Y = [];
for i = 1 : max_ticket_size
    Y(1 : max_box_size, i) = i;
end

mean_case_both = accumarray([box_size', ticket_size'], reduced_pct', [], @mean);
mean_case_only = accumarray([box_size', ticket_size'], prio_reduced_pct', [], @mean);

lower_case_both = accumarray([box_size', ticket_size'], reduced_pct', [], @(t) prctile(t, 10));
lower_case_only = accumarray([box_size', ticket_size'], prio_reduced_pct', [], @(t) prctile(t, 10));

upper_case_both = accumarray([box_size', ticket_size'], reduced_pct', [], @(t) prctile(t, 90));
upper_case_only = accumarray([box_size', ticket_size'], prio_reduced_pct', [], @(t) prctile(t, 90));

best_case_both = accumarray([box_size', ticket_size'], reduced_pct', [], @max);
best_case_only = accumarray([box_size', ticket_size'], prio_reduced_pct', [], @max);

% assume that best case could have reduction
best_to_check_idx = best_case_only == 0;

upper_case_only(best_to_check_idx) = -1;
lower_case_only(best_to_check_idx) = -1;
mean_case_only(best_to_check_idx) = -1;
upper_case_both(best_to_check_idx) = -1;
lower_case_both(best_to_check_idx) = -1;
mean_case_both(best_to_check_idx) = -1;

fig = figure;
set(fig,'Position',[200, 200, 600, 400]); 
set(gca, 'fontsize', 18);
surf(X, Y, mean_case_both);
set(gca,'xlim', [0 120]); set(gca, 'xtick', [0 : 20 : 120], 'fontsize', 15);
set(gca,'ylim', [0 400]); set(gca, 'ytick', [0 : 40 : 400], 'fontsize', 15);
set(gca,'zlim', [0 100]); set(gca, 'ztick', [0 : 20 : 100], 'fontsize', 15);
colormap([1  1  0; 0  1  1])
title('VM Reduction + Ticket Reduction')
zlabel('Percent of Reduced Tickets(%)', 'fontsize', 15); 
xlabel('Original Number of VMs per BOX', 'fontsize', 15);
ylabel('Original Number of Tickets', 'fontsize' ,15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'3dplot_reduce_ticket_pct_both'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]); 
set(gca, 'fontsize', 18);
surf(X, Y, mean_case_only);
set(gca,'xlim', [0 120]); set(gca, 'xtick', [0 : 20 : 120], 'fontsize', 15);
set(gca,'ylim', [0 400]); set(gca, 'ytick', [0 : 40 : 400], 'fontsize', 15);
set(gca,'zlim', [0 100]); set(gca, 'ztick', [0 : 20 : 100], 'fontsize', 15);
colormap([1  1  0; 0  1  1])
title('Ticket Reduction Only')
zlabel('Percent of Reduced Tickets(%)', 'fontsize', 15); 
xlabel('Original Number of VMs per BOX', 'fontsize', 15);
ylabel('Original Number of Tickets', 'fontsize' ,15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'3dplot_reduce_ticket_pct_only'));

% Ticket BOX PLOT (VM+TICKET)
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(reduced_pct', box_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('VM Reduction + Ticket Reduction')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_both'));

% Ticket BOX PLOT (TICKEY ONLY)
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(prio_reduced_pct', box_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Ticket Reduction Only')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_only'));
