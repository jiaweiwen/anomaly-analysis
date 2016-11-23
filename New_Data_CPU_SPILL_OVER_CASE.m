% This script is used to checked how often box is over-utilized
close all; clear; clc

load ../New_Data/box_vm_time_series_summary_cpu_only.mat
box_vm_time_series_summary_cpu = box_vm_time_series_summary;

% Need to check the following three things:
% 1. Sum of VM demands are alomost equal to BOX capacity;
% 2. util of BOX reaches 100%
% 3. how 1. same or different from 2.

size_box_vm = size(box_vm_time_series_summary_cpu);
vm_demands_equal_box = [];
record_vm_demands_no_less_than_box = [];
box_equal_to_100 = [];
vm_equal_to_100 = []; % how much time and the total time length

total_vm_number = 0;
spill_over_vm_number = 0;

total_box_number = 0;
spill_over_box_vm_number = 0;

equal_thres = 0.001;
for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary_cpu{1, box_id});

    % If we don't have time series
    if size_box <= 2
        continue;
    end
    
    total_box_number = total_box_number + 1;
    total_vm_number = total_vm_number + size_box -1;
    
    box_cap = box_vm_time_series_summary_cpu{1, box_id}{1,1}(:,3) ./ ...
                   box_vm_time_series_summary_cpu{1, box_id}{1,1}(:,4) * 100; 
               
    no_len = numel(box_vm_time_series_summary_cpu{1, box_id}{1,1}(:, 1));
    pm_id = box_vm_time_series_summary_cpu{1, box_id}{1,1}(1, 2);
    
    sum_vm_demands = zeros(no_len, 1);
    signal_box_spill_over = 0;
    for vm_id = 2 : size_box
        sum_vm_demands = sum_vm_demands + box_vm_time_series_summary_cpu{1, box_id}...
                                          {1,vm_id}(:, 4);
        vm_reach_100_idx = box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(:, 5) >= 100;                              
        vm_reach_100 = sum(vm_reach_100_idx);
        
        v_cap = mean(box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(:, 4) ./ ...
                     box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(:, 5) * 100);
        
        if vm_reach_100 > 0
            spill_over_vm_number = spill_over_vm_number + 1;
            vm_equal_to_100(end+1:end+vm_reach_100, 1:4) = ...
                [ones(vm_reach_100,1)*vm_reach_100, ...
                 box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(vm_reach_100_idx, 4), ...
                 ones(vm_reach_100,1)*v_cap, ones(vm_reach_100,1)* box_vm_time_series_summary_cpu{1, box_id}{1,vm_id}(1,2)];
            signal_box_spill_over = 1;
        end       
    end
    
    if signal_box_spill_over == 1
        spill_over_box_vm_number = spill_over_box_vm_number + 1;
    end
    
    box_vm_diff = (box_cap - sum_vm_demands) ./ box_cap;
    great_than_thres_idx = box_vm_diff <= equal_thres;
    sum_vm_demands_same_as_box_cap = sum(great_than_thres_idx);
    
    if sum_vm_demands_same_as_box_cap > 0
        vm_demands_equal_box(end+1, 1:3) = [sum_vm_demands_same_as_box_cap, no_len, pm_id];
        record_vm_demands_no_less_than_box(end+1 : end+sum_vm_demands_same_as_box_cap, 1:3) = ... 
        [sum_vm_demands(great_than_thres_idx), box_cap(great_than_thres_idx), ones(sum_vm_demands_same_as_box_cap,1) * pm_id];
    end
    
    box_reach_100 = sum(box_vm_time_series_summary_cpu{1, box_id}{1,1}(:,4) >= 100);
    
    if box_reach_100 > 0
        box_equal_to_100(end+1, 1:3) = [box_reach_100, no_len, pm_id];
    end
    
end

% compare *vm_demands_equal_box* and *box_equal_to_100*
if numel(box_equal_to_100) ~= 0
    [same_pm, vm_idx, box_idx] = intersect(vm_demands_equal_box(:,3) , box_equal_to_100(:,3));

    % first show the possibility of 'intersection'
    ratio_intersect_vm = numel(same_pm) / numel(vm_demands_equal_box(:,1));
    disp(strcat('The ratio of intersection for VM is ', mat2str(ratio_intersect_vm)));

    ratio_intersect_box = numel(same_pm) / numel(box_equal_to_100(:,1));
    disp(strcat('The ratio of intersection for BOX is ', mat2str(ratio_intersect_box)));

    % second check if the ratio between *the number of VM equal box cap* and *box reaches 100%*
    compare_vm_box = [vm_demands_equal_box(vm_idx, 1), box_equal_to_100(box_idx,:) ];
    compare_ratio = [compare_vm_box(:,1) ./ compare_vm_box(:,2), compare_vm_box(:,2)];

    fig = figure;
    set(fig, 'Position', [200, 200, 600, 400]);
    scatter(compare_ratio(:,2), compare_ratio(:,1));
end

fig = figure;
set(fig, 'Position', [200, 200,1200, 400]);
subplot(1,2,1);
vm_cap_ratio = vm_equal_to_100(:,2) ./ vm_equal_to_100(:,3);
vm_bin = prctile(vm_cap_ratio, 5) : (prctile(vm_cap_ratio, 95) - prctile(vm_cap_ratio, 5))/100 : prctile(vm_cap_ratio, 95);
bincounts = histc(vm_cap_ratio, vm_bin); 
plot(vm_bin, bincounts/sum(bincounts),'b-', 'linewidth', 1.5);
%set(gca, 'xtick', vm_bin);
xlabel('VM Demdans / VCap'); ylabel('PDF');
set(gca,'fontsize',15);

subplot(1,2,2);
[f_eq, x_eq] = ecdf(vm_cap_ratio); 
plot(x_eq, f_eq, 'b-', 'linewidth', 1.5);
set(gca, 'xlim', [prctile(vm_cap_ratio, 5) prctile(vm_cap_ratio, 95)]);
xlabel('VM Demdans / VCap'); ylabel('CDF');
set(gca,'fontsize',15);
