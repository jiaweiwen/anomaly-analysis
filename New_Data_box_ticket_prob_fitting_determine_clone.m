% this script aims at charactering BOX tickets per box

close all; clear; clc

load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros

path = '../New_Data_box_tickets_fit/';
mkdir(path);

mkdir(strcat(path, 'Fig/'));

ticket_thres = [60, 70, 80];

time_grat = 900;
day_point = 96;

% Prob(box)             % Prob(vm|box)
PROB_BOX_TICKET = []; PROB_VM_BOX_TICKET = [];
% related features
BOX_TICKET_FEATURE = []; VM_BOX_TICKET_FEATURE = [];

BOX_ID = [];

box_num = 1;
size_box_vm = size(box_vm_time_series_summary);

test_lag = 0;

fuzzy_bound = 0;

sig_vm_thres = 0.8;

USED_BOX = 0;

for box_id = 1 : size_box_vm(2)
    
    % feature 0: box size 
    size_box = numel(box_vm_time_series_summary{1, box_id});

    % If we don't have time series
    if size_box < 2
        continue;
    end
    
%     if box_num >= 1000
%         break;
%     end

    % at least we should have 3 days data
    total_points = numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
    if total_points < day_point * 3
        continue;
    end

    % First extract the number of tickets for CPU and RAM for different
    % thresholds
    USED_BOX = USED_BOX + 1;
    
    BOX_ID(end+1) = box_id;
    
    box = box_vm_time_series_summary{1, box_id}{1,1}(:,4);
    box_demands = box_vm_time_series_summary{1, box_id}{1,1}(:,3);
    
    % feature 1: capaciy of box
    box_capacity = nanmean(box_demands ./ box * 100);
    
    % feature set 2: mean usage and std of usage
    mean_box_usage = nanmean(box); std_box_usage = nanstd(box);
    median_box_usage = nanmedian(box); 
    box_usage_25 = prctile(box, 25); box_usage_75 = prctile(box, 75);
    
    % feature set 3: the vm total v_cap
    vm_cap = []; vm_on = []; vm_cap_on = []; vm_mean_usage = [];
    for vm_id = 2 : size_box
        non_zero_idx = box_vm_time_series_summary{1,box_id}{vm_id}(:, 4) ~= 0;
        vm_on(:, end+1) = non_zero_idx;
        vm_demands = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx, 4);
        vm_usage = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx,5);
        vm_mean_usage(end+1) = nanmean(vm_usage);
        vm_cap(end+1) = nanmean(vm_demands ./ vm_usage * 100);        
        vm_cap_on(:,end+1) = vm_on(:, end) * vm_cap(end); 
    end
    sum_vcap = sum(vm_cap);
    mean_vm_usage = nanmean(vm_mean_usage); std_vm_usage 
    
%     vm_on_sum = sum(vm_on');
%     vm_on_cap_sum = sum(vm_cap_on');
%     
%     mean_on_vm = nanmean(vm_on_sum); std_on_vm = nanstd(vm_on_sum);
%     median_on_vm = nanmedian(vm_on_sum); 
%     on_vm_25 = prctile(vm_on_sum, 25); on_vm_75 = prctile(vm_on_sum, 75);
%     
%     mean_on_cap_vm = nanmean(vm_on_cap_sum); std_on_cap_vm = nanstd(vm_on_cap_sum);
%     median_on_cap_vm = nanmedian(vm_on_cap_sum); 
%     on_cap_vm_25 = prctile(vm_on_cap_sum, 25); on_cap_vm_75 = prctile(vm_on_cap_sum, 75);
      
    for ticket_id = 1 : numel(ticket_thres)
        box_ticket = box > ticket_thres(ticket_id);
        sum_ticket = sum(box_ticket);
        % RECORD THE PROBABILITY OF BOX
        PROB_BOX_TICKET(end+1) = sum_ticket / total_points;
        BOX_TICKET_FEATURE(:, end+1) = [ticket_thres(ticket_id); ...
                                       size_box - 1; box_capacity; sum_vcap;...
                                       mean_box_usage; std_box_usage; ....
                                       median_box_usage; box_usage_25; box_usage_75];
        for vm_id = 2 : size_box
            vm = box_vm_time_series_summary{1,box_id}{vm_id}(:, end);
            if sum_ticket == 0
                PROB_VM_BOX_TICKET(end+1,1:2) = [0, box_id];
            else
                vm_ticket = vm(box_ticket) > ticket_thres(ticket_id);
                PROB_VM_BOX_TICKET(end+1,1:2) = [sum(vm_ticket) / sum_ticket, box_id];
            end
            mean_vm_usage = nanmean(vm); std_vm_usage = nanstd(vm);
            median_vm_usage = nanmedian(vm); 
            vm_usage_25 = prctile(vm, 25); vm_usage_75 = prctile(vm, 75);
            VM_BOX_TICKET_FEATURE(:, end+1) = [ticket_thres(ticket_id); ...
                                               size_box - 1; box_capacity; sum_vcap;...
                                               mean_box_usage; std_box_usage; ....
                                               median_box_usage; box_usage_25; box_usage_75; ...
                                               vm_cap(vm_id-1); mean_vm_usage; std_vm_usage; ....
                                               median_vm_usage; vm_usage_25; vm_usage_75];
        end
    end
end

save(strcat(path, 'PROB_BOX_TICKET'), 'PROB_BOX_TICKET');
save(strcat(path, 'BOX_TICKET_FEATURE'), 'BOX_TICKET_FEATURE');
save(strcat(path, 'PROB_VM_BOX_TICKET'), 'PROB_VM_BOX_TICKET');
save(strcat(path, 'VM_BOX_TICKET_FEATURE'), 'VM_BOX_TICKET_FEATURE');

disp('The training begins');

% choose the interested freatures
interet_idx = [1, 5, 6, 7, 8, 9];
BOX_TICKET_FEATURE = BOX_TICKET_FEATURE(interet_idx, :);

% First, generate the scatter plot among features;
N = numel(BOX_TICKET_FEATURE(1,:)) / 3;
mean_usage = BOX_TICKET_FEATURE(2, 1:N);
% y_name = {'std','median','25_pct','75_pct'};
for feature_id = 3 : numel(interet_idx)
    % first do the polynomial fitting
    if feature_id > 3
        degree = 1;
    else
        degree = 2;
    end
    p = polyfit(mean_usage, BOX_TICKET_FEATURE(feature_id, 1:N), degree);
    fit_val = polyval(p, mean_usage);
end

% Decision Tree for Prob(BOX)
[rtree_box, error_box, ape_box, imp_box] = DecisionTree(BOX_TICKET_FEATURE', PROB_BOX_TICKET);
[rtree_vm, error_vm, ape_vm, imp_vm] = DecisionTree(VM_BOX_TICKET_FEATURE', PROB_VM_BOX_TICKET(:,1)');

% check the high state probability error (0.5 as threshold)
prob_thres = 0.4;
vm_idx = PROB_VM_BOX_TICKET(:,1) > prob_thres;
figure;
ecdf(abs(error_vm(vm_idx)));
figure; 
ecdf(abs(ape_vm(vm_idx)));

prob_thres = 0.4;
box_idx = PROB_BOX_TICKET > prob_thres;
figure;
ecdf(abs(error_box(box_idx)));
figure; 
ecdf(abs(ape_box(box_idx)));