% VM cloning + migration
close all; clear; clc

load ../New_Data_7days/box_vm_time_series_summary_mem_only_with_zeros
box_vm_time_series_mem = box_vm_time_series_summary;
load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros

server_to_tenant = importdata('../New_Data_7days/server_to_tenant.csv');

moment_path = '../New_Data_box_moments_fit/moment_coeff/';
load (strcat(moment_path, 'MOMENT_CONF_COEFF_UPPER'));
load (strcat(moment_path, 'MOMENT_CONF_COEFF_LOWER'));

% % structure: server-tenant
% server_to_tenant = sortrows(server_to_tenant, 2);
% 
% % structure - index: tenant_id; each entry: box + vm series 
% unique_tenant = unique(server_to_tenant(:,2));
% tenant_box_series = cell(1, numel(unique_tenant));
% for box_id = 1 : numel(box_vm_time_series_summary)
%     % if no time series
%     if numel(box_vm_time_series_summary{1, box_id}) == 0
%         continue;
%     end
%     
%     pm_id = box_vm_time_series_summary{1, box_id}{1,1}(1,2);
%     
%     tenant_id = find(server_to_tenant(:,1) == pm_id);
%     if numel(tenant_id) == 0
%         continue;
%     end
%     
%     % merge memory into the cpu series
%     for vm_id = 1 : numel(box_vm_time_series_summary{1, box_id})
%         box_vm_time_series_summary{1, box_id}{1, vm_id}(:,end+1:end+2) = ...
%                   box_vm_time_series_mem{1, box_id}{1, vm_id}(:,end-1:end);
%     end
%     
%     tenant_id = find(unique_tenant == server_to_tenant(tenant_id,2));
%     % each tenant has the 
%     tenant_box_series{tenant_id}{end+1} = ...
%                                      box_vm_time_series_summary{1, box_id};   
% end
% 
% % characterize on the tenant:
% % 1. sort the gathered tenant-server by the number of servers
% [~, I] = sort(cellfun(@length, tenant_box_series), 'descend');
% tenant_box_series = tenant_box_series(I);
% unique_tenant = unique_tenant(I);

clone_path = '../New_Data_box_moments_fit/tenant_server/';
% mkdir(clone_path)
% save(strcat(moment_path, 'unique_tenant'), 'unique_tenant');
% save(strcat(moment_path, 'tenant_box_series'), 'tenant_box_series');

load (strcat(moment_path, 'unique_tenant'));
load (strcat(moment_path, 'tenant_box_series'));

load (strcat(moment_path, 'BOX_MEAN'));
load (strcat(moment_path, 'BOX_MOMENTS'));


%%%%%%%%%%%%%%%%%%% check the box memory distribution %%%%%%%%%%%%%%%%%%%%%
mean_mem_check = [];
for box_id = 1 : numel(box_vm_time_series_mem)
    if numel(box_vm_time_series_mem{box_id}) == 0
        continue;
    end
    
    mean_mem_check(end+1) = nanmean(box_vm_time_series_mem{box_id}{1}(:,end));

end
font_size = 10;
fig = figure;
set(fig, 'Position', [200 200 400 300]);
[f, x] = ecdf(mean_mem_check);
plot(x, f, 'r-', 'linewidth', 1.5);
xlabel('Mean RAM Usage'); 
ylabel('CDF');
set(gca, 'xlim', [0 100]); set(gca, 'xtick', [0:10:100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0:0.1:1]); 
set(gca, 'fontsize', font_size);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat('../New_Data_box_moments_fit/Fig/cdf_mean_ran_usage'));

%%%%%%%%%%%%%%%%%%% VM Cloning: test by the tenant %%%%%%%%%%%%%%%%%%%%%%%%
moments = 5 : 5 : 95;
moments = [moments, 98];

server_num_thres = 10;
ticket_thres = 60;
usage_tail_idx = numel(moments) - 1;

fig_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Target_Fig/');
mkdir(fig_path);

result_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Target_Results/'); 
mkdir(result_path);

% determine the threshold for mean usage
d_usage = 5;
mean_usage = 0 : d_usage : ticket_thres;
conf_upper = [];
for mean_id = 1 : numel(mean_usage)
    used_mean = mean_usage(mean_id);
    
    p2 = MOMENT_CONF_COEFF_UPPER(usage_tail_idx, 1);
    p1 = MOMENT_CONF_COEFF_UPPER(usage_tail_idx, 2);
    p0 = MOMENT_CONF_COEFF_UPPER(usage_tail_idx, 3);       
    CONF_UPPER(mean_id) = max(min(100, used_mean^2 * p2 + used_mean * p1 + p0), 0);
end

% derermine the threshold for the mean usage
idx = max(find(CONF_UPPER <= ticket_thres));
mean_usage_thres = mean_usage(idx);

% Plot the upper bound of tail prediction
font_size = 13;
fig = figure;
set(fig, 'Position', [200 200 400 300]);
plot(mean_usage, CONF_UPPER, 'ro-', 'linewidth', 1.5);
xlabel('Mean Usage'); 
ylabel(strcat('Upper Bound of Usage Tail (=', mat2str(moments(usage_tail_idx)), '%ile)'));
set(gca, 'xlim', [0 ticket_thres]); set(gca, 'xtick', mean_usage);
set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0:10:100]); 
hold on; plot(get(gca,'xlim'), [ticket_thres ticket_thres], 'b-', 'linewidth', 1.5);
h = legend('Tail Upper Bound', 'Usage Threshold');
set(h, 'Location', 'southeast');
set(gca, 'fontsize', font_size);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'tail_upper_prediction_',mat2str(moments(usage_tail_idx))));

% things to check
% 1. distribution changes or not (x)
% 2. tail changes (x)
% 3. Number of migrations + Number of clones (x)
% 4. Violation of Memory (x)
% 5. dependency still there or not (this could be deleted if not)

% todo later:
% 1. memory constraint needs to be strict (x)
% 2. clones needs to be back to the original Box
% 3. introduce the randomness in workload splitting (x)

day_thres = 5;
test_day_thres = 1;

ticket_compare = []; 
tail_compare = []; tail_violation = [];
num_of_migrate = []; num_of_clone = [];
memory_target = 95; memory_violation_compare = [];

target_usage_tail = moments(usage_tail_idx);

% assume that the splitting has randomness
% assume it follows normal distribution
split_mean = 0;
split_sigma = 0.15; 

wind_length = 96;

tenant_size = {};

BOX_TRAIN_MEAN = BOX_MEAN; BOX_TRAIN_MOMENT = BOX_MOMENTS;
% BOX_TRAIN_MEAN = []; BOX_TRAIN_MOMENT = [];
BOX_TEST_MEAN = []; BOX_TEST_MOMENT = [];
BOX_PREDICT_MEAN = []; BOX_PREDICT_MOMENT = [];

BOX_VIOLATION_RATIO_TRAIN = [];

USED_TENANT_CLUSTER = [];

for tenant_id = 1 : numel(unique_tenant)
    num_box = numel(tenant_box_series{tenant_id});
    if  num_box < server_num_thres - 6
        break;
    end
    
    problem_box = [];
    safe_box_for_clone = []; safe_box_for_nothing = [];
    
    common_time_series = tenant_box_series{tenant_id}{1}{1}(:,1);
    for box_id = 2 : num_box
        temp_time = tenant_box_series{tenant_id}{box_id}{1}(:,1);
        common_time_series = intersect(common_time_series, temp_time);
    end
    
    % we only pick the long enough data
    if numel(common_time_series) < wind_length * day_thres
        % disp(strcat('tenant id is ', mat2str(tenant_id), 'num of box is ', mat2str(num_box)));
        continue;
    end
    
    box_vm_train_series = cell(1, num_box);
    box_vm_test_series = cell(1, num_box);
    box_vm_clone_or_migrate = cell(1, num_box);  
    box_vm_test_cpu_cap = zeros(1, num_box);
    box_vm_test_mem_cap = zeros(1, num_box);
    for box_id = 1 : num_box
        temp_time = tenant_box_series{tenant_id}{box_id}{1}(:,1);
        [~, idx, ~] = intersect(temp_time, common_time_series);
        num_vm = numel(tenant_box_series{tenant_id}{box_id});
        box_vm_train_series{box_id} = {};
        box_vm_test_series{box_id} = {};
        for vm_id = 1 : num_vm
            box_vm_train_series{box_id}{vm_id} = tenant_box_series{tenant_id}{box_id}{vm_id}...
                                                                  (idx(1 : (day_thres-test_day_thres)*wind_length), :);
            box_vm_test_series{box_id}{vm_id} = tenant_box_series{tenant_id}{box_id}{vm_id}...
                                                                  (idx((day_thres-test_day_thres)*wind_length+1 : ...
                                                                       day_thres*wind_length), :);                                                 
        end
        
        %%%%%%%%%%%%%%%%%%%%% Training Tells Problems %%%%%%%%%%%%%%%%%%%%%
        mean_box_usage = round(nanmean(box_vm_train_series{box_id}{1}(:,4)) / d_usage) * d_usage;
        mean_box_cap = nanmean(box_vm_train_series{box_id}{1}(:,3) ...
                                   ./ box_vm_train_series{box_id}{1}(:,4) * 100);
        mean_box_mem_cap = nanmean(box_vm_train_series{box_id}{1}(:,5) ...
                                   ./ box_vm_train_series{box_id}{1}(:,6) * 100);                      
                               
        box_vm_test_cpu_cap(box_id) = mean_box_cap;  
        box_vm_test_mem_cap(box_id) = mean_box_mem_cap;      
            
        % we need to add the threhold for the memory usage, e.g., 90%
        rest_mem_cap = max(0, nanmean(box_vm_train_series{box_id}{1}(:,end-1) ...
                                      ./ box_vm_train_series{box_id}{1}(:,end) * 100) * memory_target / 100 ...
                                      - nanmax(box_vm_train_series{box_id}{1}(:,end-1)));
                                                  
        great_than_thres = box_vm_train_series{box_id}{1}(:, 4) > ticket_thres;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ji: LOOK HERE!! 
        % To be fair, we set the threshold to determine if the
        % Migration/Clone is needed 
        % The threshold is 85%ile
        
        % Updated Version:
        % Threshold is set based on the demands
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if sum(great_than_thres)/numel(great_than_thres) >= 1 - moments(usage_tail_idx)/100
            % determine if box need to be remove VMs                  
            if mean_box_usage > mean_usage_thres
                problem_box(end+1, 1:2) = [box_id, (mean_box_usage - mean_usage_thres) / 100 * mean_box_cap];
            else
                % attention!!!!!! some boxes may cannot decrease by d_usage
                problem_box(end+1, 1:2) = [box_id, min(d_usage, mean_box_usage) / 100 * mean_box_cap];
            end
            
            % check the ticketing points and the corresponding VM usage
            % Next Step: ticket should have several levels
            vm_ticket_point = []; 
            vm_ticket_rank_temp = zeros(sum(great_than_thres), num_vm-1);
            vm_ticket_rank = zeros(1, num_vm-1); 
            vm_mem_cost = []; vm_cpu_mean_usage_cap = [];
            for vm_id = 2 : num_vm
                % get the actual usage
                vm_ticket_point(:, end+1) = box_vm_train_series{box_id}{vm_id}(great_than_thres, end-3);
                
                % in case that: VM migration happen in real time
                non_zero_idx_mem = find(box_vm_train_series{box_id}{vm_id}(:, end-1) ~= 0);
                if numel(non_zero_idx_mem) == 0
                    vm_mem_cost(end+1) = 0;
                    vm_cpu_mean_usage_cap(end+1) = 0;
                else
                    vm_mem_cost(end+1) = nanmean(box_vm_train_series{box_id}{vm_id}(non_zero_idx_mem, end-1) ./ ...
                                                 box_vm_train_series{box_id}{vm_id}(non_zero_idx_mem, end) * 100); 
                    vm_cpu_mean_usage_cap(end+1) = nanmean(box_vm_train_series{box_id}{vm_id}(non_zero_idx_mem, end-3));                         
                end

            end
            
            % in case of some VM memory cost is zeros
            nan_idx_mem = isnan(vm_mem_cost);
            vm_mem_cost(nan_idx_mem) = 0;
%             [~, mem_idx] = sort(vm_mem_cost);
%             vm_mem_cost(mem_idx) = 1 : num_vm-1;
           
%             % rank the VM by the CPU usage
%             for ticket_id = 1 : sum(great_than_thres)
%                 [~, idx] = sort(vm_ticket_point(ticket_id,:));
%                 vm_ticket_rank_temp(ticket_id, idx) = 1 : num_vm - 1;
%             end
%             [~, idx] = sort(mean(vm_ticket_point));
%             vm_ticket_rank(idx) = 1 : num_vm-1;
            vm_ticket_rank = nanmean(vm_ticket_point);
            
            % usage low, rank low
            % lower cloning ratio, better to be cloned/migrated
            vm_clone_ratio = vm_mem_cost ./ vm_ticket_rank;
            [~, vm_clone_idx] = sort(vm_clone_ratio);
            
            % check which vm need to be cloned or migrated
            needed_migrate_cpu_cap = problem_box(end, 2);
            total_mem_cost = 0;
            for cand_vm_id = 1 : num_vm - 1
                vm_id = vm_clone_idx(cand_vm_id);
                if vm_cpu_mean_usage_cap(vm_id) <= needed_migrate_cpu_cap
                    box_vm_clone_or_migrate{box_id}(end+1, :) = [1, 1, vm_cpu_mean_usage_cap(vm_id), ...
                                                                 vm_id, vm_mem_cost(vm_id)]; % migrate all
                    needed_migrate_cpu_cap = needed_migrate_cpu_cap - vm_cpu_mean_usage_cap(vm_id);
                    total_mem_cost = total_mem_cost + vm_mem_cost(vm_id);
                else
                    box_vm_clone_or_migrate{box_id}(end+1, :) = ...
                        [2, needed_migrate_cpu_cap / vm_cpu_mean_usage_cap(vm_id), ...
                         needed_migrate_cpu_cap, vm_id, vm_mem_cost(vm_id)]; % clone
                    needed_migrate_cpu_cap = 0;
                    total_mem_cost = total_mem_cost + vm_mem_cost(vm_id);
                end
                
                if needed_migrate_cpu_cap == 0
                    break;
                end
            end
            
            % the total memory cost from this box
            problem_box(end, 3) = total_mem_cost;
                      
        else
            if mean_box_usage >= mean_usage_thres
                safe_box_for_nothing(end+1, :) = [box_id, 0];
            else
                safe_box_for_clone(end+1, :) = [box_id, (mean_usage_thres - mean_box_usage) / 100 * mean_box_cap, rest_mem_cap];   
            end
        end
        
    end    
        
    % compare the tickets
    if numel(problem_box) == 0 || numel(safe_box_for_clone) == 0
        continue;
    end
    
    % record the ratio of box with tail violation in training process
    disp(numel(problem_box(:,1)) / num_box);
    BOX_VIOLATION_RATIO_TRAIN(end+1) = numel(problem_box(:,1)) / num_box;
    
    % sort the safe box by availabel demands
    safe_box_for_clone = sortrows(safe_box_for_clone,[-2, -3]);
    
    % do the migration/cloning on the *testing* data sets
    % first, do the testing copies
    box_vm_predict_series = box_vm_test_series;
    
    problem_box = problem_box(:,1); num_prob_box = numel(problem_box);
    tag_prob_box = ones(1, num_prob_box);
    
    % initial the number of migration and cloning
    num_of_clone(end+1) = 0; num_of_migrate(end+1) = 0;
    
    % record the tested tenant cluster id
    USED_TENANT_CLUSTER(end+1) = tenant_id;
    
    while sum(tag_prob_box) ~= 0
        % everytime work on one problematic box
        for prob_box_id = 1 : num_prob_box
            box_id = problem_box(prob_box_id);       

            % if all the problematic VMs have been solved
            num_clone_vm = numel(box_vm_clone_or_migrate{box_id});
            if num_clone_vm == 0
                tag_prob_box(prob_box_id) = 0;
                continue;
            end
            
            safe_box_for_clone = sortrows(safe_box_for_clone,[-2, -3]);

            % everytime, we only clone/migrate one VM for the fairness concern
            vm_cand = 1;
            migrate_cpu_cap = box_vm_clone_or_migrate{box_id}(vm_cand, 3);
            migrate_mem_cost = box_vm_clone_or_migrate{box_id}(vm_cand, 5);

            clone_or_migrate = box_vm_clone_or_migrate{box_id}(vm_cand, 1);
            migrate_ratio = box_vm_clone_or_migrate{box_id}(vm_cand, 2);

            migrate_vm_id = box_vm_clone_or_migrate{box_id}(vm_cand, 4);

            for safe_box_id = 1 : numel(safe_box_for_clone(:,1))
                if migrate_cpu_cap <= safe_box_for_clone(safe_box_id, 2) && ...
                   migrate_mem_cost <= safe_box_for_clone(safe_box_id, 3) 

                    % do the cloning or migration
                    if clone_or_migrate == 1 % migrate
                        num_of_migrate(end) = num_of_migrate(end) + 1;
                    elseif clone_or_migrate == 2
                        num_of_clone(end) = num_of_clone(end) + 1;
                        split_factor = max(1, migrate_ratio * (1 + normrnd(split_mean, split_sigma, [test_day_thres * wind_length, 1])));
                    end                           

                    % step1: update the safe box
                    safe_box_for_clone(safe_box_id, 2) = safe_box_for_clone(safe_box_id, 2) - migrate_cpu_cap;
                    safe_box_for_clone(safe_box_id, 3) = safe_box_for_clone(safe_box_id, 3) - migrate_mem_cost;

                    % step2: update the safe box predicted trace
                    used_box_id = safe_box_for_clone(safe_box_id, 1);
                    % step2.1:
                        % add the time series first
                        % assume only the demands are cut in ratio:
                        % both CPU and RAM
                    box_vm_predict_series{used_box_id}{end+1} = box_vm_test_series{box_id}{migrate_vm_id+1};
                    if clone_or_migrate == 1
                        box_vm_predict_series{used_box_id}{end}(:, end-3) = box_vm_predict_series{used_box_id}{end}(:, end-3) * migrate_ratio; 
                        box_vm_predict_series{used_box_id}{end}(:, end-1) = box_vm_predict_series{used_box_id}{end}(:, end-1) * migrate_ratio; 
                    else
                        box_vm_predict_series{used_box_id}{end}(:, end-3) = box_vm_predict_series{used_box_id}{end}(:, end-3) .* split_factor; 
                        box_vm_predict_series{used_box_id}{end}(:, end-1) = box_vm_predict_series{used_box_id}{end}(:, end-1) .* split_factor; 
                    end

                    % step2.2:
                        % update the used box cpu and memory usage information
                    % CPU
                    box_vm_predict_series{used_box_id}{1}(:, end-3) =  box_vm_predict_series{used_box_id}{1}(:, end-3) + ...
                                                                    box_vm_predict_series{used_box_id}{end}(:, end-3);                                                                
                    box_vm_predict_series{used_box_id}{1}(:, end-2) = min(box_vm_predict_series{used_box_id}{1}(:, end-3) ...
                                                                       / box_vm_test_cpu_cap(used_box_id) * 100, 100);
                    % MEMORY             
                    %%%%% Question: Is memory splitted in ratio? %%%%%%
                    box_vm_predict_series{used_box_id}{1}(:, end-1) =  box_vm_predict_series{used_box_id}{1}(:, end-1) + ...
                                                                    box_vm_predict_series{used_box_id}{end}(:, end-1);                                                                
                    box_vm_predict_series{used_box_id}{1}(:, end) = box_vm_predict_series{used_box_id}{1}(:, end-1) ...
                                                                       / box_vm_test_mem_cap(used_box_id) * 100; 


                    % step3: update the problem box predicted trace
                    % step3.1:
                        % update the box cpu and memory usage information
                    % CPU
                    box_vm_predict_series{box_id}{1}(:, end-3) = box_vm_predict_series{box_id}{1}(:, end-3) - ...
                                                                    box_vm_predict_series{used_box_id}{end}(:, end-3); 
                    box_vm_predict_series{box_id}{1}(:, end-2) = max(0, box_vm_predict_series{box_id}{1}(:, end-3) ...
                                                                     / box_vm_test_cpu_cap(box_id) * 100);
                    % MEMORY
                    box_vm_predict_series{box_id}{1}(:, end-1) = box_vm_predict_series{box_id}{1}(:, end-1) - ...
                                                                    box_vm_predict_series{used_box_id}{end}(:, end-1);
                    box_vm_predict_series{box_id}{1}(:, end) = max(0, box_vm_predict_series{box_id}{1}(:, end-1) ...
                                                                     / box_vm_test_mem_cap(box_id) * 100);


                    % step3.2:
                        % add the time series first
                        % assume only the demands are cut in ratio: CPU
                        % and RAM
                    if clone_or_migrate == 1
                        box_vm_predict_series{box_id}{migrate_vm_id+1}(:, end-3 : end) = 0;
                    else
                        box_vm_predict_series{box_id}{migrate_vm_id+1}(:, end-3) = ...
                                              box_vm_predict_series{box_id}{migrate_vm_id+1}(:, end-3) .* (1-split_factor);
                        box_vm_predict_series{box_id}{migrate_vm_id+1}(:, end-1) = ...
                                              box_vm_predict_series{box_id}{migrate_vm_id+1}(:, end-1) .* (1-split_factor);  
                    end

                    % get rid of the loop
                    break
                end 

            end

            % step4: remove the cloned/migrated VM from the problem box     
            box_vm_clone_or_migrate{box_id} = box_vm_clone_or_migrate{box_id}(2:end, :);

        end
    end
    
    ticket_compare(end+1, 1:2) = 0;
    memory_violation_compare(end+1, 1:2) = 0;
    for box_id = 1 : num_box
        ticket_compare(end, 1) = ticket_compare(end, 1) + sum(box_vm_test_series{box_id}{1}(:, end-2) > ticket_thres);
        ticket_compare(end, 2) = ticket_compare(end, 2) + sum(box_vm_predict_series{box_id}{1}(:, end-2) > ticket_thres);
        
        % record the mean and percentile
%         BOX_TRAIN_MEAN(end+1) = nanmean(box_vm_train_series{box_id}{1}(:,4));
        BOX_TEST_MEAN(end+1) = nanmean(box_vm_test_series{box_id}{1}(:,4));
        BOX_PREDICT_MEAN(end+1) = nanmean(box_vm_predict_series{box_id}{1}(:,4));
        train_moment = []; 
        test_moment = []; predict_moment = [];
        for moment_id = 1 : numel(moments)
%             train_moment(moment_id) = prctile(box_vm_train_series{box_id}{1}(:,4), moments(moment_id));
            test_moment(moment_id) = prctile(box_vm_test_series{box_id}{1}(:,4), moments(moment_id));
            predict_moment(moment_id) = prctile(box_vm_predict_series{box_id}{1}(:,4), moments(moment_id));
        end
%         BOX_TRAIN_MOMENT(end+1, :) = train_moment;
        BOX_TEST_MOMENT(end+1, :) = test_moment;
        BOX_PREDICT_MOMENT(end+1, :) = predict_moment;    
  
        % record the memory changes
        memory_violation_compare(end, 1) = memory_violation_compare(end, 1) + sum(box_vm_test_series{box_id}{1}(:, end) > memory_target);
        memory_violation_compare(end, 2) = memory_violation_compare(end, 2) + sum(box_vm_predict_series{box_id}{1}(:, end) > memory_target);
        
    end
    
    % record the target tail
    temp_tail = [BOX_TEST_MOMENT(end - num_box + 1 : end, usage_tail_idx), ...
                 BOX_PREDICT_MOMENT(end - num_box + 1 : end, usage_tail_idx)];
    tail_compare(end+1, 1:2) = prctile(temp_tail, 90);
    tail_violation(end+1, 1:2) = sum(temp_tail > ticket_thres) / num_box;
    
    ticket_compare(end, 1:2) = ticket_compare(end, 1:2)/num_box;
    tenant_size{end+1} = mat2str(num_box);
    
    % disp('Check me pls');
%     disp('hahahahaha')
%     disp(numel(BOX_VIOLATION_RATIO_TRAIN));
%     disp(numel(ticket_compare(:,1)))
%     disp('hahahahaha')
end

save(strcat(result_path, 'USED_TENANT_CLUSTER_', moments(usage_tail_idx)), 'USED_TENANT_CLUSTER');
save(strcat(result_path, 'BOX_VIOLATION_RATIO_TRAIN_', moments(usage_tail_idx)), 'BOX_VIOLATION_RATIO_TRAIN');
save(strcat(result_path, 'tail_violation_', moments(usage_tail_idx)), 'tail_violation');
save(strcat(result_path, 'ticket_compare_', moments(usage_tail_idx)), 'ticket_compare');

%%%%%%%%%%%%%%%%%%%%%%% Check 1: Distribution Changes %%%%%%%%%%%%%%%%%%%%%
usage_size = 5;
BOX_TRAIN_MEAN_TEST = round(BOX_TRAIN_MEAN / usage_size) * usage_size;
BOX_TEST_MEAN_TEST = round(BOX_TEST_MEAN / usage_size) * usage_size;
BOX_PREDICT_MEAN_TEST = round(BOX_PREDICT_MEAN /usage_size) * usage_size;

font_size = 10;

% check the histogram of the Box usage
bins = 0 : 10 : 100;
fig = figure;
set(fig, 'Position', [200 200 450 300]);
[n_train, edges_train] = histcounts(BOX_TRAIN_MEAN_TEST, bins);
[n_test, edges_test] = histcounts(BOX_TEST_MEAN_TEST, bins);
[n_predict, edges_predict] = histcounts(BOX_PREDICT_MEAN_TEST, bins);
usage_summary = [n_test'/sum(n_test), ...
                 n_predict'/sum(n_predict)];
h = bar(usage_summary, 'grouped');
set(h(1), 'facecolor', 'r');
set(h(2), 'facecolor', 'k');
xlabel('Mean Box CPU Usage (%)'); 
ylabel('Histogram');
set(gca, 'xticklabel', bins);
l = legend('Testing Set Before', 'Testing Set After');
set(l, 'location', 'northeast');
set(gca, 'ylim', [0 0.3]); set(gca, 'ytick', [0: 0.05 : 0.3]); 
set(gca, 'fontsize', font_size);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'pdf_mean_train_test_predict_', ...
                                mat2str(moments(moment_id))));
 
% check the moments distribution
test_usage = [usage_size, 20, 25, 30, 35, 40, ticket_thres - usage_size];% 0 : usage_size : ticket_thres - usage_size;
sig_level = 0.05;
for test_id = 1 : numel(test_usage)
    train_idx = find(BOX_TRAIN_MEAN_TEST == test_usage(test_id));
    test_idx = find(BOX_TEST_MEAN_TEST == test_usage(test_id));
    predict_idx = find(BOX_PREDICT_MEAN_TEST == test_usage(test_id));
    
    % check the distribution of the moments for train/test/predict
    for moment_id = usage_tail_idx : usage_tail_idx
        [f_train, x_train] = ecdf(BOX_TRAIN_MOMENT(train_idx, moment_id));
        [f_test, x_test] = ecdf(BOX_TEST_MOMENT(test_idx, moment_id));
        [f_predict, x_predict] = ecdf(BOX_PREDICT_MOMENT(predict_idx, moment_id));
        
        [h1, p1] = kstest2(BOX_TRAIN_MOMENT(train_idx, moment_id), ...
                   BOX_TEST_MOMENT(test_idx, moment_id), 'Alpha', sig_level);
        [h2, p2] = kstest2(BOX_TRAIN_MOMENT(train_idx, moment_id), ...
                   BOX_PREDICT_MOMENT(predict_idx, moment_id), 'Alpha', sig_level);  
        [h3, p3] = kstest2(BOX_TEST_MOMENT(test_idx, moment_id), ...
                   BOX_PREDICT_MOMENT(predict_idx, moment_id), 'Alpha', sig_level); 
            
                         
        fig = figure;
        set(fig, 'Position', [200 200 450 300]);
        plot(x_train, f_train, 'b-', 'linewidth', 1.5);
        hold on
        plot(x_test, f_test, 'r--', 'linewidth', 1.5);
        hold on
        plot(x_predict, f_predict, 'k:', 'linewidth', 1.5);
        xlabel(strcat(mat2str(moments(moment_id)), '%ile')); 
        ylabel('CDF');
%         if test_id == 2
%             title1 = strcat('Mean Box CPU Usage: ', {' '}, mat2str(mean_usage_thres));
%         else
        title1 = strcat('Mean Box CPU Usage: ', {' '}, mat2str(test_usage(test_id)));

%         title2 = strcat('p_{train:test1} =', {' '}, mat2str(round(p1, 4)), ', ', ...
%                       'p_{train:test2} =', {' '}, mat2str(round(p2, 4)), ', ', ...
%                       'p_{test1:test2} =', {' '}, mat2str(round(p3, 4)));
        title2 = strcat('h_{train:test1} =', {' '}, mat2str(h1), ', ', ...
                      'h_{train:test2} =', {' '}, mat2str(h2), ', ', ...
                      'h_{test1:test2} =', {' '}, mat2str(h3));   
        title(title1{1});
        %title({title1{1}; title2{1}});
        %set(gca, 'xlim', [0 100]); set(gca, 'xtick', [0:10:100]); 
        set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0:0.1:1]);
        if test_id >= 2
            
            hold on
            plot([ticket_thres ticket_thres], get(gca,'ylim'), 'm-')
            h = legend(strcat('# of Training Set = ', mat2str(numel(train_idx))), ...
                       strcat('# of Testing Set Before = ', mat2str(numel(test_idx))), ...
                       strcat('# of Testing Set After = ', mat2str(numel(predict_idx))), ...
                       'Ticket Threshold');
        else
            h = legend(strcat('# of Training Set = ', mat2str(numel(train_idx))), ...
                       strcat('# of Testing Set Before = ', mat2str(numel(test_idx))), ...
                       strcat('# of Testing Set After = ', mat2str(numel(predict_idx))));
        end
        set(h, 'Location', 'southeast');
        set(gca, 'fontsize', font_size + 3);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(fig_path, 'dist_compare_mean_', ...
                                        mat2str(test_usage(test_id)), '_moment_', mat2str(moments(moment_id))));
    end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check 1: Finished %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Check 2: tail comparison %%%%%%%%%%%%%%%%%%%%%%%
% 90%ile of all the distri
upper = ceil(max(max(tail_compare)) / usage_size) * usage_size;
lower = floor(min(min(tail_compare)) / usage_size) * usage_size;
bins = lower : (upper - lower) / 10 : upper;
fig = figure;
set(fig, 'Position', [200 200 450 300]);
[n_test, edge_test] = histcounts(tail_compare(:,1), bins);
[n_predict, edge_predict] = histcounts(tail_compare(:,2), bins);
tail_test = [n_test'/sum(n_test), n_predict'/sum(n_predict)];
h = bar(tail_test, 'grouped');
set(h(1), 'facecolor', 'r');
set(h(2), 'facecolor', 'k');
xlabel(strcat('90%ile of ', {' '}, mat2str(moments(usage_tail_idx)), ' Target Across Tenants')); 
ylabel('Histogram');
set(gca, 'xticklabel', bins);
l = legend('Testing Set Before', 'Testing Set After');
set(l, 'location', 'northeast');
set(gca, 'ylim', [0 0.4]); set(gca, 'ytick', [0: 0.05 : 0.4]); 
set(gca, 'fontsize', font_size);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'tail_compare_test_predict_', mat2str(moments(usage_tail_idx))));

% check percentage of violation
fig = figure;
set(fig, 'Position', [200 200 450 300]);
num_case = numel(tenant_size);
[BOX_VIOLATION_RATIO_TRAIN_new, idx] = sort(BOX_VIOLATION_RATIO_TRAIN);
plot(1 : num_case, tail_violation(idx, 1) * 100, 'r*:', 'linewidth', 1.5);
hold on
plot(1 : num_case, tail_violation(idx, 2) * 100, 'ko--', 'linewidth', 1.5);
xlabel('Percent of Box w/ Tail Violation in Training Satge (%)'); 
ylabel('Ratio of Boxes w/ Tail Violation (%)');
title1 = strcat('Before: Mean Violation = ', mat2str(round(mean(tail_violation(:,1)) * 100, 2)), '%');
title2 = strcat('After: Mean Violation = ', mat2str(round(mean(tail_violation(:,2)) * 100, 2)), '%');
% set(gca, 'xlim', [0 100]); set(gca, 'xtick', [0 : 10 : 100]); 
set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:1:num_case]); 
set(gca, 'xticklabel', round(BOX_VIOLATION_RATIO_TRAIN_new * 100));
% xticklabel_rotate([], 45);
title(strcat(title1, {', '}, title2), 'fontweight', 'normal');
set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0:10:100]); 
h = legend('Original', 'Migration + Cloning');
set(h, 'Location', 'northwest');
set(gca, 'fontsize', font_size);

set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'tail_violation_compare_', mat2str(moments(usage_tail_idx))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check 2: Finished %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Check 3: # of Migrate + Clone %%%%%%%%%%%%%%%%%%%%%
% fig = figure;
% set(fig, 'Position', [200 200 450 300]);
% [ticket_tenant, idx] = sortrows([str2double(tenant_size)', round(ticket_compare(:,1))], [-1, -2]);
% stem3(1:num_case, ticket_tenant(:,2), num_of_clone(idx), 'r*:', 'linewidth', 1.5);
% hold on
% stem3(1:num_case, ticket_tenant(:,2), num_of_migrate(idx), 'ko--', 'linewidth', 1.5);
% xlabel('Tenant Size'); ylabel('Original Ticket Number')
% zlabel('# of Migration/Cloning');
% set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1: 2 :num_case]);
% set(gca, 'xticklabel', tenant_size(idx(1:2:num_case)));
% h = legend('Cloning', 'Migration');
% set(h, 'Location', 'northoutside');
% set(gca, 'fontsize', font_size -3); 
% % upper = ceil(max(max(num_of_clone), max(num_of_migrate)) / 10) * 10;
% % set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0: upper/10: upper], 'fontsize', font_size - 3); 
% set(gca, 'zscale', 'log');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(fig_path, 'num_migrate_clone_based_on_ticket', mat2str(moments(usage_tail_idx))));

disp(strcat('Tail target is ', mat2str(moments(usage_tail_idx))));
disp(strcat('Mean # of Migration is ', mat2str(round(mean(num_of_migrate),2))));
disp(strcat('Mean # of Clone is ', mat2str(round(mean(num_of_clone), 2))));

fig = figure;
set(fig, 'Position', [200 200 450 300]);
[f_migrate, x_migrate] = ecdf(num_of_migrate);
[f_clone, x_clone] = ecdf(num_of_clone);
plot(x_migrate, f_migrate, 'r-', 'linewidth', 1.5);
hold on
plot(x_clone, f_clone, 'k--', 'linewidth', 1.5);
xlabel('Number of Movement'); 
ylabel('CDF');
%set(gca, 'xlim', [1 30]); set(gca, 'xtick', [1:3:30]); 
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0:0.1:1]); 
h = legend('# of Migration', '# of Clones');
set(h, 'Location', 'southeast');
set(gca, 'fontsize', font_size);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'cdf_movement_', mat2str(moments(usage_tail_idx))));

fig = figure;
set(fig, 'Position', [200 200 450 300]);
h = bar([num_of_migrate', num_of_clone'], 'grouped');
set(h(1), 'facecolor', 'r');
set(h(2), 'facecolor', 'k');
xlabel('Tenant'); 
ylabel('# of Movement');
%set(gca, 'xlim', [1 30]); set(gca, 'xtick', [1:3:30]); 
set(gca, 'ylim', [0.1 1000]); 
set(gca, 'ytick', [0.1, 1, 10, 100,1000]);
set(gca, 'yscale', 'log');
h = legend('Migration', 'Clones');
set(h, 'Location', 'northeast');
set(gca, 'fontsize', font_size + 3);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'num_movement_', mat2str(moments(usage_tail_idx))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check 3: Finished %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Check 4: Memory Vilation Compare %%%%%%%%%%%%%%%%%%%%
num_of_move = num_of_migrate + num_of_clone;
memory_violation_compare(:,1) = memory_violation_compare(:,1) ./ str2double(tenant_size)';
memory_violation_compare(:,2) = memory_violation_compare(:,2) ./ str2double(tenant_size)';

[num_of_move, idx] = sort(num_of_move);

num_case = numel(tenant_size);

fig = figure;
set(fig, 'Position', [200 200 450 300]);
plot(1:num_case, memory_violation_compare(idx, 1), 'r*:', 'linewidth', 1.5);
hold on
plot(1:num_case, memory_violation_compare(idx, 2), 'ko--', 'linewidth', 1.5);
xlabel('Number of Migration + Clones'); 
ylabel('Memory Violation Times per Box');
title1 = strcat('Mean RAM Violation Before = ', mat2str(round(mean(memory_violation_compare(:,1)), 2)));
title2 = strcat('Mean RAM Violation After = ', mat2str(round(mean(memory_violation_compare(:,2)), 2)));
set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:1:num_case]); 
set(gca, 'xticklabel', num2cell(num_of_move));
xticklabel_rotate([], 45);
title(strcat(title1, {', '}, title2), 'fontweight', 'normal');
h = legend('Original', 'Migration + Cloning');
set(h, 'Location', 'northwest');
set(gca, 'fontsize', font_size);
% set(gca, 'ylim', [0 20]); set(gca, 'ytick', [0:2:20], 'fontsize', font_size - 3); 
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'ram_violation_compare_', mat2str(moments(usage_tail_idx))));

disp(strcat('Mean # of Increased Memory Violation is ', mat2str(round(mean(memory_violation_compare(:, 2) - memory_violation_compare(:,1)),2))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check 4: Finished %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Check 5: Ticket Reduction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% (Quantification of How many) %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% (Make up for the tail) %%%%%%%%%%%%%%%%%%%%%%
non_zero_idx = find(ticket_compare(:,1) ~= 0);
ticket_compare = ticket_compare(non_zero_idx, :);
BOX_VIOLATION_RATIO_TRAIN = BOX_VIOLATION_RATIO_TRAIN(non_zero_idx);

fig = figure;
set(fig, 'Position', [200 200 450 300]);
[ticket_tenant, idx] = sortrows([BOX_VIOLATION_RATIO_TRAIN', round(ticket_compare(:,1))], [-1, -2]);
num_case = numel(idx);
stem3(1:num_case, ticket_tenant(:,2), (ticket_compare(idx,1) - ticket_compare(idx,2)) ./ ticket_compare(idx,1) * 100, 'r*:', 'linewidth', 1.5);
xlabel('PCT of Box w/ Tail Violation(%)'); ylabel('Original Ticket Number')
zlabel('Ticket Reduction (%) per Box');
set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:2:num_case]); 
set(gca, 'xticklabel', round(ticket_tenant(1:2:end, 1) * 100));
upper = ceil(max(ticket_compare(:,1)) / 10) * 10;
set(gca, 'ylim', [0 upper]); set(gca, 'ytick', [0: upper/5: upper]); 
set(gca, 'zlim', [0 100]); set(gca, 'ztick', [0 : 20 : 100]);
title(strcat('Mean Ticket Reduction = ', mat2str(round(...
             nanmean((ticket_compare(:,1) - ticket_compare(:,2)) ./ ticket_compare(:,1) * 100), ...
             2)), '%'))
set(gca, 'fontsize', font_size); 
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'ticket_reduction_', mat2str(moments(usage_tail_idx))));

% Plot Ticket Reduction
font_size = 10;
fig = figure;
set(fig, 'Position', [200 200 450 300]);
[BOX_VIOLATION_RATIO_TRAIN_new, idx] = sort(BOX_VIOLATION_RATIO_TRAIN);
num_case = numel(idx);
plot(1:num_case, (ticket_compare(idx,1) - ticket_compare(idx,2)) ./ ticket_compare(idx,1) * 100, 'b*', 'linewidth', 1.5);
xlabel('PCT of Box w/ Tail Violation(%)'); 
ylabel('Ticket Reduction (%) per Box');
title(strcat('Mean Ticket Reduction = ', mat2str(round(...
             nanmean((ticket_compare(:,1) - ticket_compare(:,2)) ./ ticket_compare(:,1) * 100), ...
             2)), '%'))
set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:1:num_case]); 
set(gca, 'xticklabel', round(BOX_VIOLATION_RATIO_TRAIN_new * 100));
xticklabel_rotate([], 45);
set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0:10:100]); 
% h = legend('Original Box Ticket', 'Box Ticket w/ cloning');
% set(h, 'Location', 'northwest');
set(gca, 'fontsize', font_size);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'ticket_compare_', mat2str(moments(usage_tail_idx))));

font_size = 10;
fig = figure;
set(fig, 'Position', [200 200 450 300]);
[~, idx] = sort(ticket_compare(:,1));
num_case = numel(idx);
plot(1:num_case, (ticket_compare(idx,1) - ticket_compare(idx,2)) ./ ticket_compare(idx,1) * 100, 'ro', 'linewidth', 1.5);
xlabel('Original Ticket Number'); 
ylabel('Ticket Reduction (%) per Box');
title(strcat('Mean Ticket Reduction = ', mat2str(round(...
             nanmean((ticket_compare(:,1) - ticket_compare(:,2)) ./ ticket_compare(:,1) * 100), ...
             2)), '%'))
set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:1:num_case]); 
set(gca, 'xticklabel', ceil(ticket_compare(idx,1)));
xticklabel_rotate([], 45);
set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0:10:100]); 
% h = legend('Original Box Ticket', 'Box Ticket w/ cloning');
% set(h, 'Location', 'northwest');
set(gca, 'fontsize', font_size);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(fig_path, 'ticket_compare_based_on_ticket_', mat2str(moments(usage_tail_idx))));
