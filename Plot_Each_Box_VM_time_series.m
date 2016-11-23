% Pick up the interesting BOX with highly correlated VMs inside

close all; clear; clc

load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat

size_box_vm = size(box_vm_time_series_summary);

% Determine the maximum time length
grat_small = 900;

% mkdir('../DSN_section1_figures');
path = '../DSN_section1_figures/';

for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary{1, box_id});
%     
    if box_id >= 50
        break;
    end
    
    if size_box ~= 4
        continue;
    end
    
    pm_id = box_vm_time_series_summary{1, box_id}{1, 1}(1,2);
    
%     if pm_id ~= 1021577
%         continue;
%     end
    
    fig = figure;
    set(fig, 'Position', [200 200 600 400]);
    legend_name = {};
    cc = {'r', 'b', 'y','m', 'k','g'};
    time = 1 : numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
    time = time * grat_small / 3600;
    max_cpu = 0; line_style = {'-', '--', '-.', ':','-^','-o'};
    for vm_id = 2 : size_box
        
        legend_name{end+1} = strcat('VM',mat2str(vm_id-1),'-CPU');
        legend_name{end+1} = strcat('VM',mat2str(vm_id-1),'-RAM');
        if max(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5)) > max_cpu 
            max_cpu = max(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5));
        end
        plot(time, box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5), ...
             'color', cc{vm_id-1}, 'linewidth', 1.5, ...
                        'linestyle', line_style{vm_id -1});
        hold on
        plot(time, box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5), ...
             'color', cc{size_box - 1 + vm_id - 1}, 'linewidth', 1.5, ...
                        'linestyle', line_style{vm_id -1});
        hold on
    end
    
    set(gca,'ylim', [0 ceil(max_cpu/10)*10]); set(gca, 'ytick', [0 : ceil(max_cpu/10)*2 : ceil(max_cpu/10)*10]);
    set(gca, 'xlim', [0 24]); set(gca, 'xtick', [0 : 3 : 24], 'fontsize', 24);
    % plot(get(gca, 'xlim'), [60 60], 'k--', 'linewidth',2);
    h = legend(legend_name);
    set(h, 'box', 'on', 'location', 'northwest', 'fontsize', 24)
    ylabel('CPU USED PCT (%)', 'fontsize', 24); 
    xlabel('Time (hour)', 'fontsize', 26);
    % title(strcat('PM ID = ', {' '}, mat2str(pm_id)))
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path,'box_', mat2str(pm_id)));
    
%     if pm_id == 1021599 || pm_id == 611527
%         max_cpu = 100;
%         for vm_id = 1 : size_box
%             fig = figure;
%             set(fig, 'Position', [200, 200, 600, 300]);
%             time = 1 : numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
%             time = time * grat_small / 3600;
%             if vm_id ~= 1
%                 plot(time, box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5), 'r-o');
%             else
%                 plot(time, box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4), 'r-o');
%             end
%             hold on
%             set(gca,'ylim', [0 ceil(max_cpu/10)*10]); set(gca, 'ytick', [0 : ceil(max_cpu/10)*2 : ceil(max_cpu/10)*10], 'fontsize', 18);
%             set(gca, 'xlim', [0 floor(max(time))]); set(gca, 'xtick', [0 : 3 : floor(max(time))], 'fontsize', 18);
%             plot(get(gca, 'xlim'), [60 60], 'k--', 'linewidth',1);
%             ylabel('CPU USED PCT (%)', 'fontsize', 18); 
%             xlabel('Time (hour)', 'fontsize', 18);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-dpng2','-r300', strcat('../DSN_PPT/','box_', mat2str(pm_id), '_', mat2str(vm_id)));
%             
%             fig = figure;
%             set(fig, 'Position', [200, 200, 600, 300]);
%             time = 1 : numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
%             time = time * grat_small / 3600;
%             if vm_id ~= 1
%                 stairs(time, box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5) > 60, 'b');
%             else
%                 stairs(time, box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4) > 60, 'b');
%             end
%             hold on
%             set(gca,'ylim', [0 1]); set(gca, 'ytick', [0 1]);
%             set(gca, 'xlim', [0 floor(max(time))]); set(gca, 'xtick', [0 : 3 : floor(max(time))], 'fontsize', 18);
%             plot(get(gca, 'xlim'), [60 60], 'k--', 'linewidth',1);
%             ylabel('CPU USED PCT (%)', 'fontsize', 18); 
%             xlabel('Time (hour)', 'fontsize', 18);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-dpng2','-r300', strcat('../DSN_PPT/','box_ticket_', mat2str(pm_id), '_', mat2str(vm_id)));
%         end
%         
%     end

    if pm_id == 1021577
        break;
    end
end