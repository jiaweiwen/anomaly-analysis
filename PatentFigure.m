close all; clear; clc

load ../box_vm_time_series_summary_cpu_only.mat

test_box_index = 7;
util_thres = 8;

begin_idx = 1;
test_length = 96*1.5;
grat_small = 900;

mkdir('../Patent_figures');
path = '../Patent_figures/';

for machine_id = 1 : numel(box_vm_time_series_summary{test_box_index})
    fig = figure;
    set(fig,'Position',[200, 200, 500, 200]);
    time = box_vm_time_series_summary{test_box_index}{1,machine_id}(begin_idx:begin_idx+test_length,2)/3600;
    metric = box_vm_time_series_summary{test_box_index}{1,machine_id}(begin_idx:begin_idx+test_length,end-1);
    max_util = ceil(max(metric)/10)*10;
    plot(time, metric, 'ko-','markersize',3);
    hold on
    line([time(1) time(end)], [util_thres util_thres], 'color', 'r', 'linestyle', '-.')
    xlabel('Time (hour)', 'fontsize', 15); ylabel('CPU USED PCT (%)', 'fontsize', 15)
    set(gca, 'xtick', [time(1) : 3 : time(end)],'fontsize', 15);
    set(gca, 'xlim', [time(1) time(end)]);
    set(gca, 'ytick',[0 : max_util/5 : max_util],'fontsize', 15);
    set(gca, 'ylim', [0 max_util]);
    set(gcf, 'paperpositionmode', 'auto');
    print('-dpng','-r300', strcat(path,strcat('Machine_Id_', mat2str(machine_id), '_overtime_all')));
    
    fig = figure;
    set(fig,'Position',[200, 200, 500, 200]);
    
    stairs(time, metric > util_thres, 'bo-','markersize',3);
    hold on
    xlabel('Time (hour)', 'fontsize', 15); ylabel('Ticket (1) or Not (0)', 'fontsize', 15)
    set(gca, 'xtick', [time(1) : 3 : time(end)],'fontsize', 15);
    set(gca, 'xlim', [time(1) time(end)]);
    set(gca, 'ytick',[0 1],'fontsize', 15);
    set(gca, 'ylim', [0 1]);
    set(gcf, 'paperpositionmode', 'auto');
    print('-dpng','-r300', strcat(path,strcat('Machine_Id_', mat2str(machine_id), '_ticket_overtime_all')));
    
end
