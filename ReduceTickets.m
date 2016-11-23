clear; close all; clc

%% Step 1: Get the threshold for the 
load ../box_differentID
% Grat = 15 min; CPU_USED_PCT 	CPU_USED_MHZ	CPU_COUNT	MEM_USED_PCT	MEM_USED_MB	
box = machine_cell;

load ../vm_differentID
% Grat = 15 min; CPU_USED_PCT	CPU_USED_MHZ	CPU_COUNT	MEM_USED_PCT	MEM_USED_MB	
% VM contains the BOX information, so we also need to divide based on BOX
vm = machine_cell;

% Use the 80 percentile of PCT to serve as the threshold
all_box_cpu_pct = [];
for box_id = 1 : numel(box)
    all_box_cpu_pct = [all_box_cpu_pct; box{box_id}(:,3)];
end

all_vm_cpu_pct = [];
for vm_id = 1 : numel(vm)
    all_vm_cpu_pct = [all_vm_cpu_pct; vm{vm_id}(:,4)];
end

all_cpu_pct = [all_box_cpu_pct; all_vm_cpu_pct];

pct = 30;
cpu_ticket_thres = prctile(all_cpu_pct, pct)

%% Step 2: Get the full capacity of CPU for each box, then do allication
load ../box_vm_time_series_summary_cpu_only.mat

mean_cpu = []; 
box_vm_time_series_ticket_summary = {};
box_vm_time_series_cpu_capacity_demands = {};
optimal_tickets = {}; optimal_number_ticket = [];
greedy_tickets = {}; greedy_number_ticket = [];
original_number_ticket = [];
sample_day = 1; grat_small = 900; 
size_sample = sample_day * 24 * 3600 / grat_small;
for box_id = 1 : numel(box_vm_time_series_summary)
    size_box = numel(box_vm_time_series_summary{1, box_id});
    if size_box == 0
        continue
    end
    
    % Get the mean of CPU Capacity for this box
    % It's safe to do so, we have the validation for the CPU Capacity for
    % each box
    non_zeros_index = find(box_vm_time_series_summary{box_id}{1,1}(:,4));
    cpu_box = box_vm_time_series_summary{box_id}{1,1}(non_zeros_index, 4) ./...
              box_vm_time_series_summary{box_id}{1,1}(non_zeros_index,3) * 100;
    mean_cpu(box_id) = mean(cpu_box);
    
    % Step 2.1: Get the original tickets
    original_number_ticket(box_id) = 0;
    for cand_id = 1 : size_box
        all_metric_no = numel(box_vm_time_series_summary{box_id}{cand_id}(1,:));
        box_vm_time_series_ticket_summary{box_id}{cand_id}(:, 1:all_metric_no-2) = ...
            box_vm_time_series_summary{box_id}{cand_id}(:,1:all_metric_no-2);
        box_vm_time_series_ticket_summary{box_id}{cand_id}(:, all_metric_no-1) = ...
            box_vm_time_series_summary{box_id}{cand_id}(:,all_metric_no-1) >= cpu_ticket_thres;  
        box_vm_time_series_ticket_summary{box_id}{cand_id}(:, all_metric_no) = ...
            box_vm_time_series_summary{box_id}{cand_id}(:,all_metric_no); 
        original_number_ticket(box_id) = original_number_ticket(box_id) + ...
            sum(box_vm_time_series_ticket_summary{box_id}{cand_id}(1:size_sample, all_metric_no-1)); 
    end
    
    % Step 2.2: Get the capacity demands for each machine, in the optimal
    % case that don't exceed the threshold of CPU Util, which could be
    % expressed as the CPU Capacity Demands
    vm_capacity_demands = zeros(numel(box_vm_time_series_summary{box_id}{1,1}(:,1)),1);
    for cand_id = 2 : size_box
        box_vm_time_series_cpu_capacity_demands{box_id}{1,cand_id} ...
            = box_vm_time_series_summary{box_id}{cand_id}(:,end)/cpu_ticket_thres * 100;
        vm_capacity_demands = vm_capacity_demands + box_vm_time_series_summary{box_id}{cand_id}(:,end);
    end    
    % For the box 
    box_vm_time_series_cpu_capacity_demands{box_id}{1,1} ...
        = (box_vm_time_series_summary{box_id}{1,1}(:, end) - vm_capacity_demands) ...
          / cpu_ticket_thres * 100;
     
    % Step 2.3: Call the optimal (brute) solution to solve this problem
    sample = {}; counts = {};
    for cand_id = 1 : size_box
        sample{cand_id} = unique(box_vm_time_series_cpu_capacity_demands{box_id}{cand_id}(1:size_sample));
        if sample{cand_id}(1) ~= 0
            sample{cand_id} = [0; sample{cand_id}];
        end
        counts{cand_id} = hist(box_vm_time_series_cpu_capacity_demands{box_id}{cand_id}(1:size_sample), sample{cand_id});
        
        sample{cand_id} = flipud(sample{cand_id});
        counts{cand_id} = cumsum(fliplr(counts{cand_id})) - fliplr(counts{cand_id});
    end
   
    % Call the Optimal solution to solve this problem
    time_begin = clock;
    all_combination = find_optimal(1, mean_cpu(box_id) , size_box, sample, counts, {}, []);
    time_end = clock;
    time_cost = time_end - time_begin;
    
    optimal_tickets{box_id} = all_combination;
    
    sum_ticket = 0;
    for cand_id = 1 : size_box
        sum_ticket = sum_ticket + counts{cand_id}(all_combination{end}(cand_id));     
    end
    
    optimal_number_ticket(box_id) = sum_ticket;
    
    disp('********************************************************');
    disp(strcat('Optimal: PM', mat2str(box_id), ': Cost Time is ', ...
         mat2str(time_cost(end-1)), ' min', mat2str(time_cost(end)), ' sec'));
     
    % Call the greedy algorithm to reduce the tickets
    time_begin = clock;
    [candidate_greedy, solution_or_not] = greedy_find_approximate(sample, ...
                                                 counts, mean_cpu(box_id));
    time_end = clock;
    time_cost = time_end - time_begin;
    
    if solution_or_not
        greedy_tickets{box_id} = candidate_greedy;
    end
    
    sum_ticket = 0;
    for cand_id = 1 : size_box
        sum_ticket = sum_ticket + counts{cand_id}(candidate_greedy(cand_id));     
    end
    greedy_number_ticket(box_id) = sum_ticket;
    disp(strcat('Greedy: PM', mat2str(box_id), ': Cost Time is ', ...
         mat2str(time_cost(end-1)), ' min', mat2str(time_cost(end)), ' sec'));
    
end
% save('optimal_tickets','optimal_tickets');
% save('optimal_number_ticket','optimal_number_ticket')

fig = figure;
set(fig, 'Position', [200 200 600 400]); num_box = 1;
markersize = 8;
set(gca, 'fontsize', 18); ticket_name = {};
for box_id = 1 : numel(box_vm_time_series_summary)
    size_box = numel(box_vm_time_series_summary{1, box_id});
    if size_box == 0
        continue
    end
    tick_name{num_box} = mat2str(box_vm_time_series_summary{box_id}{1,1}(1,1));
    plot(num_box, original_number_ticket(box_id), 'ko','markersize',markersize);
    hold on
    plot(num_box, optimal_number_ticket(box_id), 'r*','markersize',markersize);
    hold on
    plot(num_box, greedy_number_ticket(box_id), 'b^','markersize',markersize);
    hold on
    num_box = num_box + 1;
end

h = legend('Origninal Number of Tickets', 'Optimal Number of Ticktes', 'Greedy Results of Tickts');
set(h, 'location', 'northeast', 'box','on');
set(gca,'xtick', 1 : num_box -1, 'xticklabel', tick_name);
set(gca, 'yscale', 'log');
xlabel('PM ID'); ylabel('Number of Tickets');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat('Reduce_Tickets'));


