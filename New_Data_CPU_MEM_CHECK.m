% this script is to check *SUM(capacity of VMs) and Sum(capacity of BOX)

close all; clear; clc


load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat
box_vm_time_series_summary_cpu = box_vm_time_series_summary;

mkdir('../New_Data_CPU_MEM_VM_PM_RELATION');
path = '../New_Data_CPU_MEM_VM_PM_RELATION/';

size_box_vm = size(box_vm_time_series_summary_cpu);
measure_cpu = [];
usage_capacity = [];
suspected_box_cpu = [];
vm_cap_greater_than_box = {};
vm_cap_greater_than_box_all = [];
for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary_cpu{1, box_id});

    % If we don't have time series
    if size_box <= 2
        continue;
    end
    
    len = numel(box_vm_time_series_summary_cpu{1, box_id}{1,1}(:,3));
    box_idx = box_vm_time_series_summary_cpu{1, box_id}{1,1}(1,2);
    
    % Check cpu
    measure_cpu(end+1, 1) = mean(box_vm_time_series_summary_cpu{1, box_id}{1,1}(:,3) ./ ...
                                 box_vm_time_series_summary_cpu{1, box_id}{1,1}(:,4) * 100); 
    temp_used_cpu_box = box_vm_time_series_summary_cpu{1, box_id}{1,1}(:,3);
    
    sum_cpu = 0; sum_cpu_used = 0;
    temp_used_cpu = zeros(len, 1);
    signal = 0;
    for vm_id = 2 : size_box
        vm_cap = mean(box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(:,4) ./ ...
                      box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(:,5) * 100);
        sum_cpu = sum_cpu + vm_cap;
        time_idx = 1;
        if vm_cap > measure_cpu(end, 1)
            if time_idx == 1
                vm_cap_greater_than_box{end+1}(time_idx, 1:2) = [box_id, vm_cap/ measure_cpu(end,1)];
            else
                vm_cap_greater_than_box{end}(time_idx, 1:2) = [box_id, vm_cap/ measure_cpu(end,1)];
            end
            vm_cap_greater_than_box_all(end+1) = vm_cap/ measure_cpu(end,1);
            time_idx = time_idx + 1;
        end
        
        temp_used_cpu = temp_used_cpu + box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(:,4);
    end
    measure_cpu(end, 2) = sum_cpu; 

    % Check the total usage 
    usage_capacity(end+1, 1:3) = [sum_cpu/measure_cpu(end, 1), nanmean(temp_used_cpu ./ temp_used_cpu_box), prctile(temp_used_cpu ./ temp_used_cpu_box, 90)]; 

    if nanmean(temp_used_cpu ./ temp_used_cpu_box) > 1
        suspected_box_cpu(end+1, 1:2) = [box_idx,nanmean(temp_used_cpu ./ temp_used_cpu_box)] ;
    end
end
suspected_box_cpu = sortrows(suspected_box_cpu, -2);

fig = figure;
set(fig, 'Position', [200, 200, 600, 400]);
hist(vm_cap_greater_than_box_all);
ylabel('Number of Outliers', 'fontsize', 15); 
xlabel('Per VM CPU Capacity / PM CPU Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'hist_per_vm_greater_than_box_cpu'));

fig = figure;
set(fig,'Position',[200, 200, 1200, 400]);
subplot(1,2,1);
scatter(usage_capacity(:,1), usage_capacity(:,2));
%set(gca, 'yscale', 'log');
set(gca,'ylim', [0 1]); 
%set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
% 15)
ylabel('VM CPU PCT / PM CPU PCT (MEAN)', 'fontsize', 15); 
xlabel('VM CPU Capacity / PM CPU Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);

subplot(1,2,2);
scatter(usage_capacity(:,1), usage_capacity(:,3));
%set(gca, 'yscale', 'log');
set(gca,'ylim', [0 100]); 
%set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
% 15)
ylabel('VM CPU PCT / PM CPU PCT (90%ile)', 'fontsize', 15); 
xlabel('VM CPU Capacity / PM CPU Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'scatter_box_vm_cpu'));

fig = figure;
usage_capacity(:,1) = ceil(usage_capacity(:,1));
unique_usage = unique(usage_capacity(:,1));
set(fig,'Position',[200, 200, 1200, 400]);
subplot(1,2,1);
boxplot(usage_capacity(:,2), usage_capacity(:,1));
%set(gca, 'yscale', 'log');
set(gca,'ylim', [0 1]); 
%set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
% 15)
ylabel('VM CPU PCT / PM CPU PCT (MEAN)', 'fontsize', 15); 
xlabel('VM CPU Capacity / PM CPU Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);

subplot(1,2,2);
boxplot(usage_capacity(:,3), usage_capacity(:,1));
%set(gca, 'yscale', 'log');
set(gca,'ylim', [0 1]); 
%set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
% 15)
ylabel('VM CPU PCT / PM CPU PCT (90%ile)', 'fontsize', 15); 
xlabel('VM CPU Capacity / PM CPU Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'boxplot_box_vm_cpu'));

size_box_vm = size(box_vm_time_series_summary_mem);
measure_mem = [];
usage_capacity_mem = [];
suspected_box_mem = [];
for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary_mem{1, box_id});
    
    % If we don't have time series
    if size_box <= 2
        continue;
    end

    len = numel(box_vm_time_series_summary_mem{1, box_id}{1,1}(:,3));
    box_idx = box_vm_time_series_summary_mem{1, box_id}{1,1}(1,2);
    
    % Check mem
    measure_mem(end+1, 1) = mean(box_vm_time_series_summary_mem{1, box_id}{1,1}(:,3));
    temp_used_mem_box = box_vm_time_series_summary_mem{1, box_id}{1,1}(:,3) .* ...
                        box_vm_time_series_summary_mem{1, box_id}{1,1}(:,4) / 100;
    
    sum_mem = 0;
    temp_used_mem = zeros(len, 1);
    for vm_id = 2 : size_box
        sum_mem = sum_mem + mean(box_vm_time_series_summary_mem{1, box_id}{1,vm_id}(:,4));
        temp_used_mem = temp_used_mem + box_vm_time_series_summary_mem{1, box_id}{1,vm_id}(:,4) .* ...
                                        box_vm_time_series_summary_mem{1, box_id}{1,vm_id}(:,5) / 100;
    end
    measure_mem(end, 2) = sum_mem;
    
    % Check the total usage 
    usage_capacity_mem(end+1, 1:3) = [sum_mem/measure_mem(end, 1), nanmean(temp_used_mem ./ temp_used_mem_box), prctile(temp_used_mem ./ temp_used_mem_box, 90)]; 
    
    if nanmean(temp_used_mem ./ temp_used_mem_box) > 1
        suspected_box_mem(end+1, 1:2) = [box_idx, nanmean(temp_used_mem ./ temp_used_mem_box)];
    end
end
suspected_box_mem = sortrows(suspected_box_mem, -2);

fig = figure;
set(fig,'Position',[200, 200, 1200, 400]);

subplot(1,2,1)
scatter(usage_capacity_mem(:,1), usage_capacity_mem(:,2));
set(gca, 'yscale', 'log'); set(gca, 'xscale', 'log');
% set(gca,'ylim', [0 1]); 
%set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
ylabel('VM RAM PCT / PM RAM PCT (Mean)', 'fontsize', 15); 
xlabel('VM RAM Capacity / PM RAM Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);

subplot(1,2,2)
scatter(usage_capacity_mem(:,1), usage_capacity_mem(:,3));
set(gca, 'yscale', 'log'); set(gca, 'xscale', 'log');
% set(gca,'ylim', [0 1]); %set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
ylabel('VM RAM PCT / PM RAM PCT (90%ile)', 'fontsize', 15); 
xlabel('VM RAM Capacity / PM RAM Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'scatter_box_vm_mem'));

fig = figure;
usage_capacity_mem(:,1) = ceil(usage_capacity_mem(:,1));
set(fig,'Position',[200, 200, 1200, 400]);
subplot(1,2,1); 
boxplot(usage_capacity_mem(:,2), usage_capacity_mem(:,1));
%set(gca, 'yscale', 'log');
set(gca,'ylim', [0 1]); %set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
% 15)
ylabel('VM RAM PCT / PM RAM PCT (MEAN)', 'fontsize', 15); 
xlabel('VM RAM Capacity / PM RAM Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);

subplot(1,2,2);
boxplot(usage_capacity_mem(:,3), usage_capacity_mem(:,1));
%set(gca, 'yscale', 'log');
set(gca,'ylim', [0 1]); %set(gca, 'ytick', [0 : 20 : 100], 'fontsize',
% 15)
ylabel('VM RAM PCT / PM RAM PCT (90%ile)', 'fontsize', 15); 
xlabel('VM RAM Capacity / PM RAM Capacity', 'fontsize', 15);
set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'boxplot_box_vm_mem'));


