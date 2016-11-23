% this script aims at charactering BOX tickets per box

close all; clear; clc

tenant_file = importdata('../New_Data_7days/server_to_tenant.csv');

box_tenant_map = tenant_file;

tenant_id_all = unique(box_tenant_map(:,2));

tenant_box_vm_time_series = cell(1, numel(tenant_id_all));
tenant_test_box_num = zeros(1, numel(tenant_id_all));

load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros

path = '../New_Data_box_vm_tickets_characterization_per_box_per_tenant_with_zeros/';
mkdir(path);


ticket_thres = [60, 70, 80];

time_grat = 900;
day_point = 96;

size_box_vm = size(box_vm_time_series_summary);

test_lag = 0;

fuzzy_bound = 0;

sig_vm_thres = 0.8;

BOX_TOTAL_CAP = 0;

NUMBER_TICKETS = 0;

tenant_box_num_thres = 30;

for box_id = 1 : size_box_vm(2)

    size_box = numel(box_vm_time_series_summary{1, box_id});

    % If we don't have time series
    if size_box < 3
        continue;
    end
    
%     if box_num >= 1000
%         break;
%     end

    % at least we should have 3 days data
    total_points = numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
    if total_points < 96
        continue;
    end
    
    pm_id = box_vm_time_series_summary{1, box_id}{1,1}(1,2);
    tenant_id = find(pm_id == box_tenant_map(:,1));
    if isempty(tenant_id)
        continue;
    else
        result_tenant_idx = find(box_tenant_map(tenant_id,2) == tenant_id_all);
    end
    % record the box_id corresponding to the tenant id
    tenant_box_vm_time_series{result_tenant_idx}(end+1) = box_id;
    tenant_test_box_num(result_tenant_idx) = tenant_test_box_num(result_tenant_idx) + 1;
end

% sort the tenant based on the number of servers
[~, IDX] = sort(tenant_test_box_num, 'descend');
final_box_vm_time_series = {};
idx = 1;
while tenant_test_box_num(IDX(idx)) ~= 0
    final_box_vm_time_series{idx} = tenant_box_vm_time_series{IDX(idx)};
    idx = idx + 1;
end

tenant_box_vm_time_series = final_box_vm_time_series;

%%%%%%%%%%%%%%%%% For each tenant, calculate the probability %%%%%%%%%%%%%%
number_of_tenant = numel(tenant_box_vm_time_series);
BOX_TICKET_ALL = {{},{},{}};
BOX_TICKET_ALL_GROUPED = {{},{},{}};

BOX_TICKET_INTER_ARRIVAL = {{}, {}, {}};
BOX_TICKET_INTER_ARRIVAL_ACF = {{}, {}, {}};

BOX_TICKET_INTER_ARRIVAL_GROUPED = {{}, {}, {}};
BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED = {{}, {}, {}};

BOX_TICKET_DURATION = {{}, {}, {}};

PROB_BOX_VM_TICKET = {{}, {}, {}};
BOX_TICKET_VM_NO_TICKET = {{}, {}, {}};
BOX_TICKET_VM_TICKET = {{}, {}, {}};

PROB_BOX_VM_TICKET_DIFF_NUM = {{}, {}, {}};

TENANT_BOX_SIZE = [];

TENANT_BOX_VM_TICKET = {};

for tenant_id = 1 : number_of_tenant
    
    num_box = numel(tenant_box_vm_time_series{tenant_id});
    
    if num_box < tenant_box_num_thres
        break;
    end
    
    for ticket_id = 1 : numel(ticket_thres)
        BOX_TICKET_ALL{ticket_id}{end+1} = [];
        BOX_TICKET_ALL_GROUPED{ticket_id}{end+1} = [];

        BOX_TICKET_INTER_ARRIVAL{ticket_id}{end+1} = [];
        BOX_TICKET_INTER_ARRIVAL_ACF{ticket_id}{end+1} = [];

        BOX_TICKET_INTER_ARRIVAL_GROUPED{ticket_id}{end+1} = [];
        BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED{ticket_id}{end+1} = [];

        BOX_TICKET_DURATION{ticket_id}{end+1} = [];

        PROB_BOX_VM_TICKET{ticket_id}{end+1} = [];
        BOX_TICKET_VM_NO_TICKET{ticket_id}{end+1} = [];
        BOX_TICKET_VM_TICKET{ticket_id}{end+1} = [];

        PROB_BOX_VM_TICKET_DIFF_NUM{ticket_id}{end+1} = [];
    end
    
    box_size_summary = []; box_vm_ticket_summary = {[], [], []};
    for box_idx = 1 : num_box
        % First extract the number of tickets for CPU and RAM for different
        % thresholds

        box_id = tenant_box_vm_time_series{tenant_id}(box_idx);
         
        total_points = numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
         
        box = box_vm_time_series_summary{1, box_id}{1,1}(:,4);
        box_demands = box_vm_time_series_summary{1, box_id}{1,1}(:,3);
        box_capacity = nanmean(box_demands ./ box * 100);

        time_stamp = box_vm_time_series_summary{1, box_id}{1,1}(:, 1);
        
        size_box = numel(box_vm_time_series_summary{1, box_id});
        box_size_summary(end+1) = size_box - 1;

        box_ticket_overtime = [];

        for ticket_id = 1 : numel(ticket_thres)
            box_ticket_overtime(:, ticket_id) = box > ticket_thres(ticket_id);
            sum_ticket = sum(box_ticket_overtime(:, ticket_id));
            % RECORD HOW MANY TICKETS ON EACH BOX
            BOX_TICKET_ALL{ticket_id}{tenant_id}(end+1, :) = ...
                [sum_ticket, total_points, sum_ticket/total_points * day_point, size_box - 1];

            % RECORD THE INTER-ARRIVAL TIME DISTRIBUTION ON EACH BOX
            non_zero_idx = find(box_ticket_overtime(:, ticket_id) == 1); 
            total_col = numel(non_zero_idx);
            if total_col == 0
                continue;
            end

            % group applied
            box_bursty_ticket = []; box_bursty_inter_arrival = [];
            col_num = 1; ticket_num = 1; 

            while col_num <= total_col - 1 
                time_delta = time_stamp(non_zero_idx(col_num + 1)) - time_stamp(non_zero_idx(col_num));
                if time_delta/time_grat == 1
                    ticket_num = ticket_num + 1;
                    col_num = col_num + 1;
                    continue;
                else
                    box_bursty_ticket(end+1,1:2) = [ticket_num, time_delta];
                end
                col_num = col_num + 1;
                ticket_num = 1;                
            end
            % Consider the last ticket
            box_bursty_ticket(end+1, 1:2) = [ticket_num, 0];

            BOX_TICKET_ALL_GROUPED{ticket_id}{tenant_id}(end+1, :) = ...
                       [numel(box_bursty_ticket(:,1)), total_points, ...
                        numel(box_bursty_ticket(:,1))/total_points * day_point, size_box - 1];

            % this value need to be checked the counts distribution on each box

            chosen_time = time_stamp(non_zero_idx);
            inter_arrival = (chosen_time(2:end) - chosen_time(1:end-1) - time_grat)/time_grat;
            % No group applied
            % form: mean, std, 5%ile, 95%ile, 25%ile, 75%ile
            BOX_TICKET_INTER_ARRIVAL{ticket_id}{tenant_id} = [BOX_TICKET_INTER_ARRIVAL{ticket_id}{tenant_id}; ...
                                                   nanmean(inter_arrival), ...
                                                   nanstd(inter_arrival),...
                                                   prctile(inter_arrival, 5), prctile(inter_arrival, 95), ...
                                                   prctile(inter_arrival, 25), prctile(inter_arrival, 75)];

            % group applied
            inter_arrival_group = box_bursty_ticket(:,2);
            duration = box_bursty_ticket(:,1);
            BOX_TICKET_INTER_ARRIVAL_GROUPED{ticket_id}{tenant_id} = [BOX_TICKET_INTER_ARRIVAL_GROUPED{ticket_id}{tenant_id}; ... 
                                                           nanmean(inter_arrival_group), ...
                                                           nanstd(inter_arrival_group),...
                                                           prctile(inter_arrival_group, 5), prctile(inter_arrival_group, 95), ...
                                                           prctile(inter_arrival_group, 25), prctile(inter_arrival_group, 75)];
            BOX_TICKET_DURATION{ticket_id}{tenant_id} = [BOX_TICKET_DURATION{ticket_id}{tenant_id}; ... 
                                               nanmean(duration), ...
                                               nanstd(duration),...
                                               prctile(duration, 5), prctile(duration, 95), ...
                                               prctile(duration, 25), prctile(duration, 75)];
            if numel(non_zero_idx) >= 5     
                acf_inter_arrival = autocorr(inter_arrival, min(20, floor(1/3*numel(inter_arrival))));
                [val, idx] = max(abs(acf_inter_arrival(2:end)));
                % form: max acf, lag
                BOX_TICKET_INTER_ARRIVAL_ACF{ticket_id}{tenant_id} = [BOX_TICKET_INTER_ARRIVAL_ACF{ticket_id}{tenant_id}; ...
                                                           val, idx];

                % consider grouping
                if numel(inter_arrival_group) > 3
                    acf_inter_arrival = autocorr(inter_arrival_group, min(20, floor(1/3*numel(inter_arrival_group))));
                    [val, idx] = max(abs(acf_inter_arrival(2:end)));
                    BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED{ticket_id}{tenant_id} = [BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED{ticket_id}{tenant_id}; ...
                                                                       val, idx];                  
                end            
            end

        end  

        vm_ticket_overtime = {}; 
        for vm_id = 2 : size_box
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%consider VM%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vm = box_vm_time_series_summary{1, box_id}{1,vm_id}(:,5);
            for ticket_id = 1 : numel(ticket_thres)
                temp = vm > ticket_thres(ticket_id) - fuzzy_bound;       
                vm_ticket_overtime{ticket_id}(:, vm_id - 1) = temp;
            end
        end

        vm_ticket_overtime_all = {};
        for ticket_id = 1 : numel(ticket_thres)
            vm_ticket_overtime_all{ticket_id} = sum(vm_ticket_overtime{ticket_id}');
        end

        % conditional probability
        for ticket_id = 1 : numel(ticket_thres)
            % only consider box with tickets           
            non_zero_idx = find(box_ticket_overtime(:, ticket_id) == 1);  
            
            % Prob(vm | box) and Prob(vm | box no ticket)
            if numel(non_zero_idx) ~= 0
                box_vm_ticket_summary{ticket_id}(end+1,1:2) = [numel(non_zero_idx); 0];
                
                % Prob(vm | box)
                sum_vm_ticket = sum(vm_ticket_overtime_all{ticket_id}(non_zero_idx) > 0);
                prob_vm_box = sum_vm_ticket / numel(non_zero_idx);
                
                % Prob(vm | box no ticket)
                zero_idx = find(box_ticket_overtime(:, ticket_id) == 0);  
                sum_vm_ticket_box_no = sum(vm_ticket_overtime_all{ticket_id}(zero_idx) > 0);
                prob_vm_box_no = sum_vm_ticket_box_no / numel(zero_idx);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%% Special Issue %%%%%%%%%%%%%%%%%%%%%%%%%
                % box has tickets and no VM has tickets
                vm_no_ticket_id = find(vm_ticket_overtime_all{ticket_id}(non_zero_idx) == 0);
                if numel(vm_no_ticket_id) ~= 0
                    util = []; demands = [];
                    mean_vcap = [];
                    for vm_id = 2 : size_box
                        util(:, vm_id - 1) = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx(vm_no_ticket_id), 5);
                        demands(:, vm_id - 1) = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx(vm_no_ticket_id), 4);

                        mean_vcap(vm_id - 1) = nanmean(box_vm_time_series_summary{1, box_id}{vm_id}(:,4) ./ ...
                                                       box_vm_time_series_summary{1, box_id}{vm_id}(:,5) * 100);
                    end

                    % vm demands ranking
                    dom_vm_num = [];
                    for col_id = 1 : numel(vm_no_ticket_id)
                        test_demand = sort(demands(col_id, :), 'descend');
                        cumsum_demand = cumsum(test_demand);
                        cumsum_ratio = cumsum_demand / cumsum_demand(end); 
                        [~, idx] = max(cumsum_ratio >= sig_vm_thres);
                        dom_vm_num(end+1) = idx;
                    end

                    % several value
                    mean_util = nanmean(util'); std_util = nanstd(util');
                    max_util = nanmax(util'); min_util = nanmin(util');
                    upper_util = prctile(util', 90); lower_util = prctile(util', 10);
                    total_demands = sum(demands');
                    vm_box_demands_ratio = total_demands ./ box_vm_time_series_summary{1, box_id}{1,1}(non_zero_idx(vm_no_ticket_id), 3)';

                    % one value
                    total_vcap = sum(mean_vcap); 
                    vm_box_cap_ratio = total_vcap / box_capacity;

                    % FORMAT: mean util, std util, V_d / B_d,
                    %         V_c / B_c, dominant VM number, Number of VM,
                    %         max util, max/min, upper util, upper/lower
                    BOX_TICKET_VM_NO_TICKET{ticket_id}{tenant_id} = [BOX_TICKET_VM_NO_TICKET{ticket_id}{tenant_id}; ...
                                                          mean_util', std_util', vm_box_demands_ratio', ...
                                                          ones(numel(vm_no_ticket_id),1) * vm_box_cap_ratio, ...
                                                          dom_vm_num', ones(numel(vm_no_ticket_id),1) * (size_box-1), ...
                                                          max_util', max_util' ./ min_util', upper_util', upper_util' ./ lower_util']; 
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%% Common Case %%%%%%%%%%%%%%%%%%%%%%%%%
                % box has tickets and VM has tickets
                vm_ticket_id = find(vm_ticket_overtime_all{ticket_id}(non_zero_idx) ~= 0);
                if numel(vm_ticket_id) ~= 0
                    util = []; demands = [];
                    mean_vcap = [];
                    for vm_id = 2 : size_box
                        util(:, vm_id - 1) = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx(vm_ticket_id), 5);
                        demands(:, vm_id - 1) = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx(vm_ticket_id), 4);

                        mean_vcap(vm_id - 1) = nanmean(box_vm_time_series_summary{1, box_id}{vm_id}(:,4) ./ ...
                                                       box_vm_time_series_summary{1, box_id}{vm_id}(:,5) * 100);
                    end

                    % vm demands ranking
                    dom_vm_num = [];
                    for col_id = 1 : numel(vm_ticket_id)
                        test_demand = sort(demands(col_id, :), 'descend');
                        cumsum_demand = cumsum(test_demand);
                        cumsum_ratio = cumsum_demand / cumsum_demand(end); 
                        [~, idx] = max(cumsum_ratio >= sig_vm_thres);
                        dom_vm_num(end+1) = idx;
                    end

                    % several value
                    mean_util = nanmean(util'); std_util = nanstd(util');
                    max_util = nanmax(util'); min_util = nanmin(util');
                    upper_util = prctile(util', 90); lower_util = prctile(util', 10);
                    total_demands = sum(demands');
                    vm_box_demands_ratio = total_demands ./ box_vm_time_series_summary{1, box_id}{1,1}(non_zero_idx(vm_ticket_id), 3)';

                    % one value
                    total_vcap = sum(mean_vcap); 
                    vm_box_cap_ratio = total_vcap / box_capacity;

                    % FORMAT: mean util, std util, V_d / B_d,
                    %         V_c / B_c, dominant VM number, Number of VM
                    BOX_TICKET_VM_TICKET{ticket_id}{tenant_id} = [BOX_TICKET_VM_TICKET{ticket_id}{tenant_id}; ...
                                                      mean_util', std_util', vm_box_demands_ratio', ...
                                                      ones(numel(vm_ticket_id),1) * vm_box_cap_ratio, ...
                                                      dom_vm_num', ones(numel(vm_ticket_id),1) * (size_box-1), ...
                                                      max_util', max_util' ./ min_util', upper_util', upper_util' ./ lower_util'];
                end


                % calculate Prob(box | vm) and Prob(box | vm no ticket)
                % Prob(box | vm)
                non_zero_idx = find(vm_ticket_overtime_all{ticket_id} > 0); 
                box_vm_ticket_summary{ticket_id}(end, 2) = numel(non_zero_idx);
                sum_box_ticket = sum(box_ticket_overtime(non_zero_idx, ticket_id) == 1);
                prob_box_vm = sum_box_ticket / numel(non_zero_idx);
                
                % Prob(box | vm no ticket)
                zero_idx = find(vm_ticket_overtime_all{ticket_id} == 0);  
                sum_box_ticket_no = sum(box_ticket_overtime(zero_idx, ticket_id) == 1);
                prob_box_vm_no = sum_box_ticket_no / numel(zero_idx);

                % record the sum of VM with tickets at the stamp for box ticket
                sum_vm_ticket = ceil(vm_ticket_overtime_all{ticket_id}(non_zero_idx) / (size_box-1) * 10) * 10;
                box_ticket_or_not = box_ticket_overtime(non_zero_idx, ticket_id);
                vm_ratio = 10:10:100; mean_prob_box_vm = -ones(1,numel(vm_ratio));
                for vm_ratio_id = 1 : numel(vm_ratio)
                    idx = find(sum_vm_ticket == vm_ratio(vm_ratio_id));
                    if numel(idx) ~= 0
                        mean_prob_box_vm(vm_ratio_id) = sum(box_ticket_or_not(idx)) / numel(idx);   
                        mean_prob_box_vm(vm_ratio_id) = ceil(mean_prob_box_vm(vm_ratio_id) * 50) /50;
                    end
                end

                

                % format: Prob(vm|box), Prob(box|vm), size_vm
                PROB_BOX_VM_TICKET{ticket_id}{tenant_id}(end+1, :) = [prob_vm_box, prob_vm_box_no, ...
                                                                      prob_box_vm, prob_box_vm_no, size_box - 1];

                PROB_BOX_VM_TICKET_DIFF_NUM{ticket_id}{tenant_id}(end+1, :) = mean_prob_box_vm;
            end
        end  
    end
    
    TENANT_BOX_SIZE(end+1, 1:2) = [mean(box_size_summary), std(box_size_summary)];
    for ticket_id = 1 : numel(ticket_thres)
        TENANT_BOX_VM_TICKET{ticket_id}(end+1, 1:4) = [mean(box_vm_ticket_summary{ticket_id}(:,1));...
                                                       std(box_vm_ticket_summary{ticket_id}(:,1));...
                                                       mean(box_vm_ticket_summary{ticket_id}(:,2));...
                                                       std(box_vm_ticket_summary{ticket_id}(:,2))];
    end
end

save(strcat(path, 'TENANT_BOX_SIZE'), 'TENANT_BOX_SIZE');

save(strcat(path, 'BOX_TICKET_ALL'), 'BOX_TICKET_ALL');
save(strcat(path, 'BOX_TICKET_ALL_GROUPED'), 'BOX_TICKET_ALL_GROUPED');
save(strcat(path, 'BOX_TICKET_INTER_ARRIVAL'), 'BOX_TICKET_INTER_ARRIVAL');
save(strcat(path, 'BOX_TICKET_INTER_ARRIVAL_ACF'), 'BOX_TICKET_INTER_ARRIVAL_ACF');

save(strcat(path, 'BOX_TICKET_INTER_ARRIVAL_GROUPED'), 'BOX_TICKET_INTER_ARRIVAL_GROUPED');
save(strcat(path, 'BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED'), 'BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED');
save(strcat(path, 'BOX_TICKET_DURATION'), 'BOX_TICKET_DURATION');

save(strcat(path, 'PROB_BOX_VM_TICKET'), 'PROB_BOX_VM_TICKET');
save(strcat(path, 'BOX_TICKET_VM_NO_TICKET'), 'BOX_TICKET_VM_NO_TICKET');
save(strcat(path, 'BOX_TICKET_VM_TICKET'), 'BOX_TICKET_VM_TICKET');

save(strcat(path, 'PROB_BOX_VM_TICKET_DIFF_NUM'), 'PROB_BOX_VM_TICKET_DIFF_NUM');