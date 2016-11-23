% Show an example of how our greedy algorithm works
close all; clc; clear
load ../BOX_PREDICTION_DATA/all_prediction.mat
load ../BOX_PREDICTION_DATA/all_target.mat
box_prediction = all_prediction;
box_target = all_target;

load ../VM_PREDICTION_DATA/all_prediction_bagging.mat
load ../VM_PREDICTION_DATA/all_target_bagging.mat
vm_prediction = all_prediction;
vm_target = all_target;

load ../measure_capacity_box_vm_cpu.mat

box_id_prediction = 54;
box_id_measure = 53;
font_size = 15;
path = '../Patent_figures/';

no_machine = numel(measure_capacity_box_vm{box_id_measure});

%% Get the common part
time_thres = 2*96;
common_index = [box_target{box_id_measure, 1}(2,1:time_thres)]; 
for machine_id = 1 : no_machine-1
    [value, ia, ib] = intersect(common_index, vm_target{box_id_prediction,machine_id}(2, :));
    common_index = common_index(ia);
end

%% Get demands for each machine
util_thres = 0.2;
box_vm_time_series_cpu_pct_demands = {};
box_vm_time_series_cpu_capacity_demands = {};
box_vm_time_series_cpu_pct_demands_prediction = {};
box_vm_time_series_cpu_capacity_demands_prediction = {};
[value, ia, ib] = intersect(box_target{box_id_measure, 1}(2,:), common_index);
% Actual demands
box_vm_time_series_cpu_pct_demands{box_id_measure}{1} = box_target{box_id_measure,1}(1, ia')';
box_vm_time_series_cpu_capacity_demands{box_id_measure}{1} = box_target{box_id_measure,1}(1, ia')' ...
                  / 100  * measure_capacity_box_vm{box_id_measure}{1,1}(1) / util_thres;
% Predicted demands
box_vm_time_series_cpu_pct_demands_prediction{box_id_measure}{1} = box_prediction{box_id_measure,1}(1, ia')';
box_vm_time_series_cpu_capacity_demands_prediction{box_id_measure}{1} = box_prediction{box_id_measure,1}(1, ia')' ...
                  / 100  * measure_capacity_box_vm{box_id_measure}{1,1}(1) / util_thres;
              
sum_vm = zeros(numel(common_index), 1);
sum_vm_prediction = zeros(numel(common_index), 1);
for machine_id = 1 : no_machine-1
    [value, ia, ib] = intersect(vm_target{box_id_prediction, 1}(2,:), common_index);
    % Actual demands
    box_vm_time_series_cpu_pct_demands{box_id_measure}{machine_id+1} = vm_target{box_id_prediction,machine_id}(1,ia')';
    box_vm_time_series_cpu_capacity_demands{box_id_measure}{machine_id+1} = vm_target{box_id_prediction,machine_id}(1,ia')' ...
                   / 100  * measure_capacity_box_vm{box_id_measure}{1,machine_id+1}(1) / util_thres;
               
    % Predicted demands
    box_vm_time_series_cpu_pct_demands_prediction{box_id_measure}{machine_id+1} = vm_prediction{box_id_prediction,machine_id}(1,ia');
    box_vm_time_series_cpu_capacity_demands_prediction{box_id_measure}{machine_id+1} = vm_prediction{box_id_prediction,machine_id}(1,ia')' ...
                   / 100  * measure_capacity_box_vm{box_id_measure}{1,machine_id+1}(1) / util_thres;
               
    sum_vm = sum_vm + box_vm_time_series_cpu_capacity_demands{box_id_measure}{machine_id+1};
    sum_vm_prediction = sum_vm_prediction + box_vm_time_series_cpu_capacity_demands_prediction{box_id_measure}{machine_id+1};
end
box_vm_time_series_cpu_capacity_demands{box_id_measure}{1} = box_vm_time_series_cpu_capacity_demands{box_id_measure}{1}...
    - sum_vm;
box_vm_time_series_cpu_capacity_demands_prediction{box_id_measure}{1} =  box_vm_time_series_cpu_capacity_demands_prediction{box_id_measure}{1} ...
     - sum_vm_prediction;
 
fig = figure;
set(fig, 'Position', [200 200 1200 900]);
time = [1 : numel(common_index)]*900/(3600);
for machine_id = 1 : no_machine
    subplot(no_machine, 1, machine_id)
    start = time(1); end_time = time(end);
    plot(time, box_vm_time_series_cpu_pct_demands{box_id_measure}{machine_id}, 'ro-', ...
         'markersize', 2);
    hold on
    plot(time, box_vm_time_series_cpu_pct_demands_prediction{box_id_measure}{machine_id}, 'k*-', ...
         'markersize', 2);
    hold on
    plot([start end_time], [util_thres * 100 util_thres*100], 'm--')
    xlabel('Time (Hour)'); ylabel('CPU USED PCT (%)');
    h = legend('Actual Workload', 'Predicted Workload', 'Ticket Threshold');
    set(h, 'box', 'off', 'location', 'northeast');
    %set(gca, 'xlabel', 'fontsize', font_size); set(gca, 'ylabel', 'fontsize', font_size);
    set(gca, 'xlim', [0 end_time]);
    set(gca, 'xtick', [0 : 3 : end_time]);
    set(gca, 'ylim', [0 100]);
    if machine_id == 1
        title('BOX');
    else
        title(strcat('VM', mat2str(machine_id-1)));
    end
end
set(gcf, 'paperpositionmode', 'auto');
print('-dpng','-r300', strcat(path,'prediction_overtime_all'));

% Step 2.3: Call the optimal (brute) solution to solve this problem
sample = {}; counts = {};
for cand_id = 1 : no_machine
    sample{cand_id} = unique(box_vm_time_series_cpu_capacity_demands_prediction{box_id_measure}{cand_id});
    if sample{cand_id}(1) ~= 0
        sample{cand_id} = [0; sample{cand_id}];
    end
    counts{cand_id} = hist(box_vm_time_series_cpu_capacity_demands_prediction{box_id_measure}{cand_id}, sample{cand_id});

    sample{cand_id} = flipud(sample{cand_id});
    counts{cand_id} = cumsum(fliplr(counts{cand_id})) - fliplr(counts{cand_id});
end

actual_case = {}; actual_counts = {};
for cand_id = 1 : no_machine
    actual_case{cand_id} = unique(box_vm_time_series_cpu_capacity_demands{box_id_measure}{cand_id});
    if actual_case{cand_id}(1) ~= 0
        actual_case{cand_id} = [0; actual_case{cand_id}];
    end
    actual_counts{cand_id} = hist(box_vm_time_series_cpu_capacity_demands{box_id_measure}{cand_id}, actual_case{cand_id});

    actual_case{cand_id} = flipud(actual_case{cand_id});
    actual_counts{cand_id} = cumsum(fliplr(actual_counts{cand_id})) - fliplr(actual_counts{cand_id});
end

% Call the Optimal solution to solve this problem
% time_begin = clock;
% all_combination = find_optimal(1, measure_capacity_box_vm{box_id_measure}{1,1}(1) , no_machine, sample, counts, {}, []);
% time_end = clock;
% time_cost = time_end - time_begin;
% 
% optimal_tickets{box_id_measure} = all_combination;
% 
% sum_ticket = 0;
% for cand_id = 1 : no_machine
%     sum_ticket = sum_ticket + counts{cand_id}(all_combination{end}(cand_id));     
% end
% 
% optimal_number_ticket(box_id_measure) = sum_ticket;
% 
% disp('********************************************************');
% disp(strcat('Optimal: PM', mat2str(box_id_measure), ': Cost Time is ', ...
%      mat2str(time_cost(end-1)), ' min', mat2str(time_cost(end)), ' sec'));

% Call the greedy algorithm to reduce the tickets
time_begin = clock;
[candidate_greedy, solution_or_not] = greedy_find_approximate(sample, ...
                                             counts, measure_capacity_box_vm{box_id_measure}{1,1}(1));
time_end = clock;
time_cost = time_end - time_begin;

if solution_or_not
    greedy_tickets{box_id_measure} = candidate_greedy;
end

sum_ticket = 0;
sum_allocation = 0;
current_allocation = [];
for cand_id = 1 : no_machine
    sum_ticket = sum_ticket + counts{cand_id}(candidate_greedy(cand_id));
    sum_allocation = sum_allocation + sample{cand_id}(candidate_greedy(cand_id));
    current_allocation = [current_allocation, sample{cand_id}(candidate_greedy(cand_id))];
end
left_cpu = measure_capacity_box_vm{box_id_measure}{1,1}(1) - sum_allocation;
re_allo = left_cpu * current_allocation / sum(current_allocation);

% one more step, re-allocate the left part of CPU from greedy algorithm
sum_actual_ticket = 0;
for cand_id = 1 : no_machine
    index_ticket = find(actual_case{cand_id} <= re_allo(cand_id) + sample{cand_id}(candidate_greedy(cand_id)),1);
    sum_actual_ticket = sum_actual_ticket + actual_counts{cand_id}(index_ticket)   
end
disp(strcat('The actual ticket number after reduction is ', mat2str(sum_actual_ticket)));

% Previously, the tickets number 
sum_original_ticket = 0;
for cand_id = 1 : no_machine
    sum_original_ticket = sum_original_ticket + ...
        sum(box_vm_time_series_cpu_pct_demands{box_id_measure}{cand_id} > (util_thres *  100));
end
disp(strcat('The original ticket number is ', mat2str(sum_original_ticket)));

greedy_number_ticket(box_id_measure) = sum_ticket;
disp(strcat('Greedy: PM', mat2str(box_id_measure), ': Cost Time is ', ...
     mat2str(time_cost(end-1)), ' min', mat2str(time_cost(end)), ' sec'));
 
fig = figure;
set(fig, 'Position', [200 200 1200 900]);
time = [1 : numel(common_index)]*900/(3600);
for machine_id = 1 : no_machine
    subplot(no_machine, 1, machine_id)
    start = time(1); end_time = time(end);
    if machine_id == 1       
        plot(time, box_vm_time_series_cpu_pct_demands{box_id_measure}{machine_id} , 'bo-', ...
            'markersize', 2);
    else
        plot(time, min(100,(box_vm_time_series_cpu_pct_demands{box_id_measure}{machine_id} * measure_capacity_box_vm{box_id_measure}{1,machine_id}(1)) ...
            / (re_allo(machine_id) + sample{machine_id}(candidate_greedy(machine_id)))) , 'bo-', ...
            'markersize', 2);
    end
    hold on
    plot([start end_time], [util_thres * 100 util_thres*100], 'm--')
    xlabel('Time (Hour)'); ylabel('CPU USED PCT (%)');
    h = legend('Actual Workload', 'Ticket Threshold');
    set(h, 'box', 'off', 'location', 'northeast');
    %set(gca, 'xlabel', 'fontsize', font_size); set(gca, 'ylabel', 'fontsize', font_size);
    set(gca, 'xlim', [0 end_time]);
    set(gca, 'xtick', [0 : 3 : end_time]);
    set(gca, 'ylim', [0 100]);
    if machine_id == 1
        title('BOX');
    else
        title(strcat('VM', mat2str(machine_id-1)));
    end
end
set(gcf, 'paperpositionmode', 'auto');
print('-dpng','-r300', strcat(path,'reallocation_overtime_all'));
 




