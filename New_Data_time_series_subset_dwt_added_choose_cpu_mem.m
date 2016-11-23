% This script is used to contains *2* results
% 1. No reduction for VM at all
% 2. Reduce VM number using: CORRELATION + VIF + STEPWISE

close all; clear; clc

load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat

size_box_vm = size(box_vm_time_series_summary);

% Determine the maximum time length
max_time = 96;
grat_small = 900;
max_lag = 1;
xcf_thres = 0.7;
vif_thres = 4;

percentile_to_check = 90;
grat = 40;

SEPERATE_VALUE = false;

mkdir('../New_Data_time_series_linear_subset_vif_cpu+mem_figure');
path = '../New_Data_time_series_linear_subset_vif_cpu+mem_figure/';

BOX_ID = []; ORIGINAL_VM_NUM = []; REDUCED_VM_NO = [];

ALL_ABS_ERROR = []; ALL_APE = []; ALL_R_SQUARE = [];
ALL_ABS_ERROR_MEM = []; ALL_APE_MEM = []; ALL_R_SQUARE_MEM = [];

BOX_APE_ALL_VM = []; BOX_R_SQUARE_ALL_VM = [];
BOX_APE_ALL_VM_MEM = []; BOX_R_SQUARE_ALL_VM_MEM = [];

ORIGINAL_TICKET = []; 
ORIGINAL_TICKET_MEM = [];

PRIO_RESIZE_TICKET = [];
PRIO_RESIZE_TICKET_MEM = [];

RESIZE_TICKET = [];
RESIZE_TICKET_MEM = [];

ticket_thres = 60;
ticket_thres_mem = 40;

for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary{1, box_id});
    
    % If we don't have time series
    if size_box <= 2
        continue;
    end
%     
    if box_id > 100
        break;
    end
%     
    if numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1)) <= 70
        continue;
    end
    
    %disp(strcat('sizebox', mat2str(size_box-1)));
    
    % for the mem, we change it as the format of cpu
    for vm_id = 1 : size_box
        box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, end-1) = ...
            box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, end-1) .* ...
            box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, end) / 100;
        nan_idx = isnan(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, end-1));
        box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(nan_idx, end-1) = 0;
    end
    
    % First, extract all of the time series
    time = box_vm_time_series_summary{1, box_id}{1, 1}(:,1) - box_vm_time_series_summary{1, box_id}{1, 1}(1,1);
    pm_id = box_vm_time_series_summary{1, box_id}{1, 1}(1,2);
    
    all_time_series = []; %box_vm_time_series_summary{1, box_id}{1, 1}(:,4);
    all_time_series_used_cap = [];
    all_time_series_mem = [];
    all_time_series_used_cap_mem = [];
    
    original_allocation_each_vm = [];
    original_allocation_each_vm_mem = [];
    
    box_time_series = box_vm_time_series_summary{1, box_id}{1, 1}(:,4);
    box_time_series_mem = box_vm_time_series_summary_mem{1, box_id}{1, 1}(:,4);
    
    total_available_capacity = 0;
    total_available_capacity_mem = 0;
    
    true_sample = {}; true_counts = {};
    original_ticket = 0;
    
    true_sample_mem = {}; true_counts_mem = {};
    original_ticket_mem = 0;
    
    for vm_id = 2 : size_box
        all_time_series = [all_time_series, ...
                    box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5)];
        all_time_series_used_cap = [all_time_series_used_cap, ...
                    box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4)];
                
        all_time_series_mem = [all_time_series_mem, ...
                    box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5)];
        all_time_series_used_cap_mem = [all_time_series_used_cap_mem, ...
                    box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 4)];
                
        original_allocation = mean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4) ./ ...
                                   box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5) * 100);                              
        original_allocation_mem = mean(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 4) ./ ...
                                   box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5) * 100);  
        original_allocation_mem(isnan(original_allocation_mem)) = 0;
        
        original_allocation_each_vm = [original_allocation_each_vm, original_allocation]; 
        original_allocation_each_vm_mem = [original_allocation_each_vm_mem, original_allocation_mem]; 
        
        total_available_capacity = total_available_capacity + original_allocation;
        total_available_capacity_mem = total_available_capacity_mem + original_allocation_mem;
        
        %%%%%%%%%%%%%%%%%%Find the ticket of CPU%%%%%%%%%%%%%%%%%%%%%%%%%%
        box_vm_time_series_demands = box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4) / ticket_thres * 100;
        true_sample{end+1} = unique(box_vm_time_series_demands);
        if true_sample{vm_id - 1}(1) ~= 0
            true_sample{vm_id -1} = [0; true_sample{vm_id - 1}];
        end
        true_counts{end+1} = hist(box_vm_time_series_demands, true_sample{vm_id-1});
        
        true_sample{vm_id-1} = flipud(true_sample{vm_id-1});
        true_counts{vm_id-1} = cumsum(fliplr(true_counts{vm_id-1})) - fliplr(true_counts{vm_id-1});
        
        % Find the first one less than the original allocation
        [~, idx] = min(true_sample{vm_id - 1} >= original_allocation);
        original_ticket = original_ticket + true_counts{vm_id -1}(idx);
        
        %%%%%%%%%%%%%%%%%%Find the ticket of MEM%%%%%%%%%%%%%%%%%%%%%%%%%%
        box_vm_time_series_demands = box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 4) / ticket_thres_mem * 100;
        true_sample_mem{end+1} = unique(box_vm_time_series_demands);
        if true_sample_mem{vm_id - 1}(1) ~= 0
            true_sample_mem{vm_id -1} = [0; true_sample_mem{vm_id - 1}];
        end
        true_counts_mem{end+1} = hist(box_vm_time_series_demands, true_sample_mem{vm_id-1});
        
        true_sample_mem{vm_id-1} = flipud(true_sample_mem{vm_id-1});
        true_counts_mem{vm_id-1} = cumsum(fliplr(true_counts_mem{vm_id-1})) - fliplr(true_counts_mem{vm_id-1});
        
        % Find the first one less than the original allocation
        [~, idx] = min(true_sample_mem{vm_id - 1} >= original_allocation_mem);
        original_ticket_mem = original_ticket_mem + true_counts{vm_id -1}(idx);
    end
    
    % Ticket summary
    ORIGINAL_TICKET(end+1) = original_ticket;
    ORIGINAL_TICKET_MEM(end+1) = original_ticket_mem;
    
    %%%%%%%%%%%%%CPU: PRIOR-KNOWLEDGE TICKET REDUCTION RESULTS%%%%%%%%%%%%%
    [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(true_sample, true_counts, total_available_capacity);
    sum_ticket = 0; each_left_for_vm = left_capacity / (size_box-1);
    for vm_id = 1 : size_box - 1
        allocation = true_sample{vm_id}(candidate_greedy(vm_id)) + each_left_for_vm;
        [~, idx] = min(true_sample{vm_id} >= allocation);    
        sum_ticket = sum_ticket + true_counts{vm_id}(candidate_greedy(vm_id));     
    end
    PRIO_RESIZE_TICKET(end+1) = sum_ticket;
    
    %%%%%%%%%%%%%MEM: PRIOR-KNOWLEDGE TICKET REDUCTION RESULTS%%%%%%%%%%%%%
    [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(true_sample_mem, true_counts_mem, total_available_capacity_mem);
    sum_ticket = 0; each_left_for_vm = left_capacity / (size_box-1);
    for vm_id = 1 : size_box - 1
        allocation = true_sample_mem{vm_id}(candidate_greedy(vm_id)) + each_left_for_vm;
        [~, idx] = min(true_sample_mem{vm_id} >= allocation);  
        if numel(true_counts_mem{vm_id}) ~= 0
            sum_ticket = sum_ticket + true_counts_mem{vm_id}(candidate_greedy(vm_id));    
        end
    end
    PRIO_RESIZE_TICKET_MEM(end+1) = sum_ticket;
    
    % Second, calculate the dissimilarities among these time series using
    % cross-correlation. Need to notice that, here we change the xcf into
    % the range of [0,1], which makes it more sense to measure the
    % dissimilarity

    % Merge CPU and MEM time series trace
    all_time_series = [all_time_series, all_time_series_mem];
    all_time_series_used_cap = [all_time_series_used_cap, all_time_series_used_cap_mem];
    original_allocation_each_vm = [original_allocation_each_vm, original_allocation_each_vm_mem];
    
    total_metric_num = size(all_time_series, 2);
    d = ones(total_metric_num); d_half = zeros(total_metric_num);
    cross_all = ones(total_metric_num);
    for row = 2 : total_metric_num
        for col = 1 : row -1    
            [xcf, lags]= crosscorr(all_time_series(:, row), all_time_series(:, col), max_lag);
            d(row, col) = abs(xcf(2));
            d(col, row) = d(row, col);
            d_half(row, col) = d(row, col);                     
        end
    end    
    
    % Step 1: Do clustering and pick up the representative for each cluster 
    great_than_thres = d >= xcf_thres;
    overall_xcf = [];
    for vm_no = 1 : total_metric_num
        chosen_xcf = d(vm_no, great_than_thres(vm_no, :));
        overall_xcf(end+1, 1:2) = [sum(great_than_thres(vm_no,:)), mean(chosen_xcf)];
    end
    
    [overall_xcf, idx] = sortrows(overall_xcf, [-1,-2]);
    
    % Test
    original_overall_xcf = overall_xcf; original_idx = idx;
    
    xcf_cluster = {};
    representative_cluster = [];
    representative_cluster_time_series = [];
    vm_no = 1;
    while vm_no <= total_metric_num
        if overall_xcf(vm_no, 1) == 0 
            vm_no = vm_no + 1; continue;
        end
        
        xcf_cluster{end+1} = [idx(vm_no)];
        representative_cluster(end+1) = idx(vm_no);
        representative_cluster_time_series(:,end+1) = all_time_series(:, idx(vm_no));
        % exclude itself from the candidates in the cluster
        d(idx(vm_no),idx(vm_no)) = 0;
        cand_in_cluster = find(d(idx(vm_no),:) >= xcf_thres);
        xcf_cluster{end} = [xcf_cluster{end}, cand_in_cluster];
        
        % change all the related column to '0' in xcf matrix
        d(idx(vm_no), :) = 0; d(:, idx(vm_no)) = 0;
        
        % change all the chosen candidate as '0'
        for cand_id = 1 : numel(cand_in_cluster)
            idx_in_overall_xcf = find(idx == cand_in_cluster(cand_id));
            overall_xcf(idx_in_overall_xcf, 1:2) = 0;
            d(cand_in_cluster(cand_id), :) = 0; 
            d(:, cand_in_cluster(cand_id)) = 0;
        end
        
        vm_no = vm_no + 1;
    end
    
    % Step 2: Calculate VIF for each representative in all clusters, if we
    % find significant VIF values, we regard them as the sign of
    % *multi-collinearity*, then we pick these *risky collinear
    % representatives* up to do stepwise fit in Step 3
    risky_representative_cluster = [];
    risky_representative_cluster_time_series = [];
    final_representative_cluster = [];
    VIF = [];
    cand_id = 1;
    time_len = numel(all_time_series(:,1));
    while cand_id <= numel(representative_cluster) 
        Y = representative_cluster_time_series(:, cand_id);
        X = [ones(time_len,1), representative_cluster_time_series(:, 1: cand_id -1), ...
             representative_cluster_time_series(:, cand_id+1:end)];
        [b, bint, r, rint, stats] = regress(Y,X);
        
        VIF(end+1) = 1/(1-stats(1));
        
        if VIF(end) > vif_thres
            risky_representative_cluster(end+1) = representative_cluster(cand_id);
            risky_representative_cluster_time_series(:,end+1) = Y;
        else
            final_representative_cluster(end+1) = representative_cluster(cand_id);
        end
                  
        cand_id = cand_id + 1;
    end
    
    % Step 3: Do stepwise choose, the representative in each cluster are
    % grouped in a big representative cluster. Then get the best subset of
    % the representative cluster
    
    cand_id = numel(risky_representative_cluster);
    while cand_id >= 1 
        Y = risky_representative_cluster_time_series(:, cand_id);
        X = [risky_representative_cluster_time_series(:, 1: cand_id -1), ...
             risky_representative_cluster_time_series(:, cand_id+1:end)];
        [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X, Y, 'display', 'off');
        
        if sum(inmodel) >= 2
            risky_representative_cluster_time_series(:, cand_id) = 0;
        else
            final_representative_cluster(end+1) = risky_representative_cluster(cand_id);
        end    
        
        cand_id = cand_id - 1;
    end
    % final_representative_cluster = fliplr(final_representative_cluster);
       
    % Step 4: fitting all the time series using the representatives
    X = [ones(time_len, 1), all_time_series(:, final_representative_cluster)];
    
    abs_error = []; ape = []; r_square = [];
    abs_error_mem = []; ape_mem = []; r_square_mem = [];
   
    sample = {}; counts = {};
    %sample_mem = {}; counts_mem = {};
    for cand_id = 1 : numel(representative_cluster)
        X_fit = X;
        % if the time series is in the representative
        vm_no = representative_cluster(cand_id);
        if vm_no <= total_metric_num/2
            used_thres = ticket_thres;
        else
            used_thres = ticket_thres_mem;
        end
            
        if numel(find(final_representative_cluster == vm_no)) == 1
            
            box_vm_time_series_demands = all_time_series_used_cap(:, vm_no) / used_thres * 100;   
            sample{vm_no} = unique(box_vm_time_series_demands);
            if sample{vm_no}(1) ~= 0
                sample{vm_no} = [0; sample{vm_no}];
            end
            counts{vm_no} = hist(box_vm_time_series_demands, sample{vm_no});
            sample{vm_no} = flipud(sample{vm_no});
            counts{vm_no} = cumsum(fliplr(counts{vm_no})) - fliplr(counts{vm_no});
            
            xcf_cluster{cand_id} = xcf_cluster{cand_id}(2:end);
            if SEPERATE_VALUE
                X_fit = [ones(time_len,1), all_time_series(:, vm_no)];
            end
        end       
        
        for cand_in_cluster = 1 : numel(xcf_cluster{cand_id})
            vm_no = xcf_cluster{cand_id}(cand_in_cluster);
            Y = all_time_series(:, vm_no);
            [b,bint,r,rint,stats] = regress(Y,X_fit);

            y_fit = max(0, Y - r);
            box_vm_time_series_demands = original_allocation_each_vm(vm_no)/used_thres * y_fit;
            sample{vm_no} = unique(box_vm_time_series_demands);
            if sample{vm_no}(1) ~= 0
                sample{vm_no} = [0; sample{vm_no}];
            end
            counts{vm_no} = hist(box_vm_time_series_demands, sample{vm_no});
            sample{vm_no} = flipud(sample{vm_no});
            counts{vm_no} = cumsum(fliplr(counts{vm_no})) - fliplr(counts{vm_no});
            
            if vm_no <= total_metric_num / 2
                abs_error(end+1) = mean(abs(r));
                ape(end+1) = nanmean(abs(r) ./ Y);
                r_square(end+1) = stats(1); 
            else
                abs_error_mem(end+1) = mean(abs(r));
                ape_mem(end+1) = nanmean(abs(r) ./ Y);
                r_square_mem(end+1) = stats(1); 
            end
        end
  
    end
        
    % Split MEM and CPU
    sample_mem = {sample{1, total_metric_num/2+1:end}}; 
    counts_mem = {counts{1, total_metric_num/2+1:end}};
    
    sample = {sample{1, 1:total_metric_num/2}};
    counts = {counts{1, 1:total_metric_num/2}};
    
    %disp(strcat('totalmetric', mat2str(total_metric_num)));
    
    % Step 5: fit box time series using the representative
    % For CPU
    [b, bint, r, rint, stats] = regress(box_time_series, X);
    box_abs_error = mean(abs(r));
    box_ape = nanmean(abs(r) ./ box_time_series);
    box_r_square = stats(1);
    
    % For MEM
    [b, bint, r, rint, stats] = regress(box_time_series_mem, X);
    box_abs_error_mem = mean(abs(r));
    box_ape_mem = nanmean(abs(r) ./ box_time_series);
    box_r_square_mem = stats(1);
    
    %%%%%%%%%%%CPU: Step 6: Ticket reduction results with fitting%%%%%%%%%%
    [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(sample, counts, total_available_capacity);
    sum_ticket = 0; each_left_for_vm = left_capacity / (size_box-1);
    for vm_id = 1 : size_box - 1
        allocation = sample{vm_id}(candidate_greedy(vm_id)) + each_left_for_vm;
        [~, idx] = min(true_sample{vm_id} >= allocation);     
        sum_ticket = sum_ticket + true_counts{vm_id}(idx);     
    end
    RESIZE_TICKET(end+1) = sum_ticket;
    
    %%%%%%%%%%%MEM: Step 6: Ticket reduction results with fitting%%%%%%%%%%
    [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(sample_mem, counts_mem, total_available_capacity_mem);
    sum_ticket = 0; each_left_for_vm = left_capacity / (size_box-1);
    for vm_id = 1 : size_box - 1
        allocation = sample_mem{vm_id}(candidate_greedy(vm_id)) + each_left_for_vm;
        [~, idx] = min(true_sample_mem{vm_id} >= allocation); 
        if numel(true_counts_mem{vm_id}) ~= 0
            sum_ticket = sum_ticket + true_counts_mem{vm_id}(idx);  
        end
    end
    RESIZE_TICKET_MEM(end+1) = sum_ticket;
    
    % write the summary of this box information after linear fitting
    BOX_ID(end+1) = pm_id; 
    ORIGINAL_VM_NUM(end+1) = total_metric_num; 
    REDUCED_VM_NO(end+1) = numel(final_representative_cluster);
    
    %%%%%%%%%%%%%%%%%%%% CPU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALL_ABS_ERROR(:, end+1) = [box_abs_error; mean(abs_error); median(abs_error); prctile(abs_error,percentile_to_check)]; 
    ALL_APE(:,end+1) = [box_ape; mean(ape); median(ape); prctile(ape,percentile_to_check)]; 
    ALL_R_SQUARE(:,end+1) = [box_r_square; mean(r_square); median(r_square); prctile(r_square,percentile_to_check)];
    
    %%%%%%%%%%%%%%%%%%%% MEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALL_ABS_ERROR_MEM(:, end+1) = [box_abs_error_mem; mean(abs_error_mem); median(abs_error_mem); prctile(abs_error_mem,percentile_to_check)]; 
    ALL_APE_MEM(:,end+1) = [box_ape_mem; mean(ape_mem); median(ape_mem); prctile(ape_mem,percentile_to_check)]; 
    ALL_R_SQUARE_MEM(:,end+1) = [box_r_square_mem; mean(r_square_mem); median(r_square_mem); prctile(r_square_mem,percentile_to_check)];
    
    %%%%%%%%%%%%%CPU: Step 7: fit box time series using all the VMs%%%%%%%%
    X = [ones(time_len, 1), all_time_series(:, 1:total_metric_num/2)];
    [b, bint, r, rint, stats] = regress(box_time_series, X);
    box_abs_error = mean(abs(r));
    box_ape = nanmean(abs(r) ./ box_time_series);
    box_r_square = stats(1);
    BOX_APE_ALL_VM(end+1) = box_ape;
    BOX_R_SQUARE_ALL_VM(end+1) = box_r_square;
    
    %%%%%%%%%%%%%MEM: Step 7: fit box time series using all the VMs%%%%%%%%
    X = [ones(time_len, 1), all_time_series(:, total_metric_num/2+1:end)];
    [b, bint, r, rint, stats] = regress(box_time_series_mem, X);
    box_abs_error = mean(abs(r));
    box_ape = nanmean(abs(r) ./ box_time_series);
    box_r_square = stats(1);
    BOX_APE_ALL_VM_MEM(end+1) = box_ape;
    BOX_R_SQUARE_ALL_VM_MEM(end+1) = box_r_square;
end
% 
% SUMMARY_METRICS = [BOX_ID', ORIGINAL_VM_NUM', ORIGINAL_TICKET',...
%                    ALL_APE(1:2,:)', ALL_R_SQUARE(1:2,:)',...
%                    ((ORIGINAL_VM_NUM - REDUCED_VM_NO) ./ ORIGINAL_VM_NUM)', ...
%                    ((ORIGINAL_TICKET - PRIO_RESIZE_TICKET) ./ ORIGINAL_TICKET)', ...
%                    ((ORIGINAL_TICKET - RESIZE_TICKET) ./ ORIGINAL_TICKET)'];
% SUMMARY_METRICS = sortrows(SUMMARY_METRICS, [-2, -3]);
% metric_name = {{'BOX ID'}, {'Original # of VMs'}, {'Original # of tickets'}, ...
%                {'APE (BOX/VM)'}, {'R2 (BOX/VM)'}, {'VM REDUCTION PERCENTAGE'}, ...
%                {'PRIO TICKET REDUCTION PERCENTAGE'}, {'FITTING TICKET REDUCTION PERCENTAGE'}};
% dlmwrite('All_metrics_vif.txt', metric_name);
% dlmwrite('All_metrics_vif.txt', SUMMARY_METRICS, '-append');

% CPU
save(strcat(path, 'BOX_APE_BASELINE'), 'BOX_APE_ALL_VM');
save(strcat(path, 'BOX_R_SQUARE_BASELINE'), 'BOX_R_SQUARE_ALL_VM');
save(strcat(path, 'BOX_VM_APE_COX_VIF_STEPWISE'), 'ALL_APE');
save(strcat(path, 'BOX_VM_R_SQUARE_COX_VIF_STEPWISE'), 'ALL_R_SQUARE');

% MEM 
save(strcat(path, 'BOX_APE_BASELINE_MEM'), 'BOX_APE_ALL_VM_MEM');
save(strcat(path, 'BOX_R_SQUARE_BASELINE_MEM'), 'BOX_R_SQUARE_ALL_VM_MEM');
save(strcat(path, 'BOX_VM_APE_COX_VIF_STEPWISE_MEM'), 'ALL_APE_MEM');
save(strcat(path, 'BOX_VM_R_SQUARE_COX_VIF_STEPWISE_MEM'), 'ALL_R_SQUARE_MEM');

% Number of Predictors
VM_REDUCED_PCT = (ORIGINAL_VM_NUM-REDUCED_VM_NO) ./ ORIGINAL_VM_NUM;
save(strcat(path, 'VM_REDUCED_PCT_COX_VIF_STEPWISE'), 'VM_REDUCED_PCT');
save(strcat(path, 'ORIGINAL_VM_NUM'), 'ORIGINAL_VM_NUM');

% CPU
not_zero_idx = ORIGINAL_TICKET > 0;
PRIO_REDUCED_PCT = (ORIGINAL_TICKET(not_zero_idx) - PRIO_RESIZE_TICKET(not_zero_idx))...
                   ./ ORIGINAL_TICKET(not_zero_idx);
REDUCED_PCT = (ORIGINAL_TICKET(not_zero_idx) - RESIZE_TICKET(not_zero_idx))...
               ./ ORIGINAL_TICKET(not_zero_idx);
save(strcat(path, 'PRIO_REDUCED_TICKET_PCT'), 'PRIO_REDUCED_PCT');
save(strcat(path, 'REDUCED_TICKET_PCT'), 'REDUCED_PCT');

% MEM
not_zero_idx_mem = ORIGINAL_TICKET_MEM > 0;
PRIO_REDUCED_PCT_MEM = (ORIGINAL_TICKET_MEM(not_zero_idx_mem) - PRIO_RESIZE_TICKET_MEM(not_zero_idx_mem))...
                   ./ ORIGINAL_TICKET_MEM(not_zero_idx_mem);
REDUCED_PCT_MEM = (ORIGINAL_TICKET_MEM(not_zero_idx_mem) - RESIZE_TICKET_MEM(not_zero_idx_mem))...
               ./ ORIGINAL_TICKET_MEM(not_zero_idx_mem);
save(strcat(path, 'PRIO_REDUCED_TICKET_PCT_MEM'), 'PRIO_REDUCED_PCT_MEM');
save(strcat(path, 'REDUCED_TICKET_PCT_MEM'), 'REDUCED_PCT_MEM');

% Box plot: vm reduction
box_size = floor(ORIGINAL_VM_NUM / grat + 1) * grat;
vm_reduced_pct = VM_REDUCED_PCT * 100;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(vm_reduced_pct', box_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
%set(gca, 'xlim', [0 120]); set(gca, 'xtick', [0 : 20 : 120], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('Reduced Percentage of Used Predictors (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_vm_percent'));

% Box plot: ticket Reduction w/ Original Tickets
% VM + Ticket
grat = 30;
ticket_size = floor(ORIGINAL_TICKET(not_zero_idx) / grat + 1) * grat;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(REDUCED_PCT'*100, ticket_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Predictor Reduction + Ticket Reduction')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of Tickets', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_both_cpu'));

% TICKEY ONLY
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(PRIO_REDUCED_PCT'*100, ticket_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Ticket Reduction Only')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of Tickets', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_only_cpu'));

% Box plot: ticket Reduction w/ BOX size
% VM + Ticket
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(REDUCED_PCT'*100, box_size(not_zero_idx)'/2);
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Predictor Reduction + Ticket Reduction')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_both_with_box_size_cpu'));

% TICKEY ONLY
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(PRIO_REDUCED_PCT'*100, box_size(not_zero_idx)'/2);
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Ticket Reduction Only')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_only_with_box_size_cpu'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Box plot: ticket Reduction w/ Original Tickets
% VM + Ticket
grat = 30;
ticket_size = floor(ORIGINAL_TICKET_MEM(not_zero_idx_mem) / grat + 1) * grat;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(REDUCED_PCT_MEM'*100, ticket_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Predictor Reduction + Ticket Reduction')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of Tickets', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_both_mem'));

% TICKEY ONLY
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(PRIO_REDUCED_PCT_MEM'*100, ticket_size');
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Ticket Reduction Only')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of Tickets', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_only_mem'));

% Box plot: ticket Reduction w/ BOX size
% VM + Ticket
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(REDUCED_PCT_MEM'*100, box_size(not_zero_idx_mem)'/2);
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Predictor Reduction + Ticket Reduction')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_both_with_box_size_mem'));

% TICKEY ONLY
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
boxplot(PRIO_REDUCED_PCT_MEM'*100, box_size(not_zero_idx_mem)'/2);
set(gca,'ylim', [0 100]); set(gca, 'ytick', [0 : 20 : 100], 'fontsize', 15);
title('Ticket Reduction Only')
ylabel('Reduced Percentage of Tickets on VMs (%)', 'fontsize', 15); 
xlabel('Original Number of VMs', 'fontsize', 15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'box_plot_reduce_ticket_pct_only_with_box_size_mem'));
