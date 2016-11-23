close all; clear; clc

load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat

time_grat = 900;

dom_vm_thres = 0.6;

ticket_thres = [60, 70, 80];

size_box_vm = size(box_vm_time_series_summary);

BOX_WITH_TICKETS = {[], [], []};

for box_id = 1 : size_box_vm(2)

    size_box = numel(box_vm_time_series_summary{1, box_id});

    % If we don't have time series
    if size_box < 2
        continue;
    end

    % hope they have the same the data at the same time
    % which is better for our evaluation
    total_points = numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
    if total_points ~= 96 
        continue;
    end
    
    % Then we could assume all these boxes have the same data
    box_util = box_vm_time_series_summary{1, box_id}{1,1}(:,end);
    box_demands = box_vm_time_series_summary{1, box_id}{1,1}(:,end-1);
    box_cap = nanmean(box_demands ./ box_util * 100);  
    for ticket_id = 1 : numel(ticket_thres)
        box_ticket = box_util > ticket_thres(ticket_id);
        
        if sum(box_ticket) == 0
            break;
        end
        
        % find the dom_vm_thres
        vm_ticket_id = find(box_ticket == 1);
        % box_chosen_demands = box_demands(vm_ticket_id);
        vm_cand = [];
        for vm_no = 2 : size_box
            vm_cand(end+1) = sum(box_vm_time_series_summary{box_id}{1, vm_no}(vm_ticket_id,end-1));
        end
        vm_cand = sort(vm_cand, 'descend');
        vm_cand_cumsum = cumsum(vm_cand) / sum(vm_cand);
        [~, idx] = max(vm_cand_cumsum >= dom_vm_thres);
        
        BOX_WITH_TICKETS{ticket_id} =  [BOX_WITH_TICKETS{ticket_id}; ...
                                        box_id, sum(box_ticket), mean(box_util), ...
                                        box_cap, idx, idx/(size_box-1)] ;
    end 
end

% sorting the box by their ticket number and mean box utilization, and
% number of dominant VMs and ratios
for ticket_id = 1 : numel(ticket_thres)
    BOX_WITH_TICKETS{ticket_id} = sortrows(BOX_WITH_TICKETS{ticket_id}, [5, -2, -3, 6]);
    
    % check how often the dominant VMs are less than 1 or 2
    fig = figure;
    set(fig, 'Position', [200 200 300 200])
    unique_dom_vm = unique(BOX_WITH_TICKETS{ticket_id}(:, 5));
    histogram(BOX_WITH_TICKETS{ticket_id}(:, 5), unique_dom_vm);
    title(mat2str(ticket_thres(ticket_id)))
    set(gca, 'fontsize', 10);
    
    fig = figure;
    set(fig, 'Position', [200 200 300 200])
    unique_dom_vm = unique(BOX_WITH_TICKETS{ticket_id}(:, 6));
    histogram(BOX_WITH_TICKETS{ticket_id}(:, 6), unique_dom_vm);
    title(mat2str(ticket_thres(ticket_id)))
    set(gca, 'fontsize', 10);
    
end


% Pick up the similar BOXes
ticket_number_upper = 50; ticket_number_lower = 10;
vm_number_upper = 3; vm_ratio_upper = 0.35;
vm_number_lower = 2;
for ticket_id = 2 : 2
    box_cand = BOX_WITH_TICKETS{ticket_id};
    
    cand_idx1 = find(box_cand(:, 2) < ticket_number_upper);
    cand_idx2 = find(box_cand(:, 2) > ticket_number_lower);
    cand_idx3 = find(box_cand(:, 5) <= vm_number_upper);
    cand_idx4 = find(box_cand(:, 5) >= vm_number_lower);
    
    cand_idx5 = find(box_cand(:, 6) > vm_ratio_upper);
    
    idx = intersect(cand_idx5, intersect(intersect(cand_idx1, cand_idx4), intersect(cand_idx2, cand_idx3)));
    
    test_box = box_cand(idx, :);    
    
    if numel(test_box(:,1)) >= 3
        test_box = test_box(1:3,:);
    end
    
    % set the ticket_threshold == 75% to see if migration could works
    [original_ticket, reduced_ticket, original_vm_position, final_vm_position, ...
     optimal_ticket, optimal_weight, optimal_copy] = ...
                                migration(box_vm_time_series_summary_mem, ...
                                          box_vm_time_series_summary, test_box, 75);
    
end

