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
grat = 10;

SEPERATE_VALUE = false;

mkdir('../New_Data_time_series_linear_subset_dwt_figure');
path = '../New_Data_time_series_linear_subset_dwt_figure/';

BOX_ID = []; ORIGINAL_VM_NUM = []; REDUCED_VM_NO = [];
ALL_ABS_ERROR = []; ALL_APE = []; ALL_R_SQUARE = [];

BOX_APE_ALL_VM = []; BOX_R_SQUARE_ALL_VM = [];

ORIGINAL_TICKET = []; 
PRIO_RESIZE_TICKET = [];
RESIZE_TICKET = [];

UNDER_PROVISION = [];
PRIO_UNDER_PROVISION = [];
FIT_UNDER_PROVISION = [];

ticket_thres = 60;

for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary{1, box_id});
    
    % If we don't have time series
    if size_box <= 2
        continue;
    end
%     
%     if box_id > 1000
%         break;
%     end
%     
    if numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1)) <= 70
        continue;
    end
    
    % First, extract all of the time series
    time = box_vm_time_series_summary{1, box_id}{1, 1}(:,1) - box_vm_time_series_summary{1, box_id}{1, 1}(1,1);
    pm_id = box_vm_time_series_summary{1, box_id}{1, 1}(1,2);
    all_time_series = []; %box_vm_time_series_summary{1, box_id}{1, 1}(:,4);
    all_time_series_used_cap = [];
    original_allocation_each_vm = [];
    box_time_series = box_vm_time_series_summary{1, box_id}{1, 1}(:,4);
    
    total_available_capacity = 0; 
    true_sample = {}; true_counts = {}; true_lower_bound = [];
    original_ticket = 0;
    for vm_id = 2 : size_box
        all_time_series = [all_time_series, ...
                    box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5)];
        all_time_series_used_cap = [all_time_series_used_cap, ...
                    box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4)];
                
        original_allocation = mean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4) ./ ...
                                   box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5) * 100);
        
        original_allocation_each_vm = [original_allocation_each_vm, original_allocation];                      
        
        total_available_capacity = total_available_capacity + original_allocation;
        
        box_vm_time_series_demands = box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4) / ticket_thres * 100;
        true_sample{end+1} = unique(box_vm_time_series_demands);
        if true_sample{vm_id - 1}(1) ~= 0
            true_sample{vm_id -1} = [0; true_sample{vm_id - 1}];
        end
        true_counts{end+1} = hist(box_vm_time_series_demands, true_sample{vm_id-1});
        
        true_sample{vm_id-1} = flipud(true_sample{vm_id-1});
        true_counts{vm_id-1} = cumsum(fliplr(true_counts{vm_id-1})) - fliplr(true_counts{vm_id-1});
        
        % Add the lower bound to the case, which in case that the new
        % allocation will not furthur degrade the performance, especially
        % for the *truncated* case.
        lower_bound = max(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4));
        [~, idx_lower] = min(true_sample{vm_id-1} >= lower_bound);
        
        true_lower_bound(end+1) = lower_bound;
        true_sample{vm_id - 1} = [true_sample{vm_id - 1}(1 : idx_lower-1); lower_bound];
        true_counts{vm_id - 1} = true_counts{vm_id - 1}(1 : idx_lower);      
        
        % Due to allocation to each VM cannot exceed the total capacity in
        % BOX, so the first step is to change *true_sample* and
        % *true_counts*
        [~, idx] = min(true_sample{vm_id - 1} > box_available_capacity);
        if idx ~= 1
            true_sample{vm_id - 1} = true_sample{vm_id - 1}(idx:end);
            true_counts{vm_id - 1} = true_counts{vm_id - 1}(idx:end);
        end
            
        % Find the first one less than the original allocation
        [~, idx] = min(true_sample{vm_id - 1} > original_allocation);
        original_ticket = original_ticket + true_counts{vm_id -1}(idx);
    end
    
    % Ticket summary
    ORIGINAL_TICKET(end+1) = original_ticket;
    
    % PRIOR-KNOWLEDGE TICKET REDUCTION RESULTS
    [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(true_sample, true_counts, total_available_capacity);
    if ~solution_or_not
        disp('It can never happen');
    end
    sum_ticket = 0; each_left_for_vm = left_capacity / (size_box-1);   
    exceed_allocation = [];
    allocation = [];
    for vm_id = 1 : size_box - 1
        allocation(vm_id) = true_sample{vm_id}(candidate_greedy(vm_id)) + each_left_for_vm;
        exceed_allocation(vm_id) = allocation(vm_id) - box_available_capacity;
    end
    
    if sum(exceed_allocation) >= 0 
        % This means each should at least have the *box available capacity*
        % however, this case should be very few
        when_all_resize_reach_box_cap = when_all_resize_reach_box_cap + 1;
        allocation = ones(size_box - 1) * box_available_capacity;
    else
        pos_idx = exceed_allocation >= 0;
        sum_pos = sum(exceed_allocation(pos_idx));
        allocation(pos_idx) = box_available_capacity;
        
        neg_idx = find(exceed_allocation < 0);
        neg_val = exceed_allocation(neg_idx);
        
        re_size_neg = [neg_idx', abs(neg_val)'];
        re_size_neg = sortrows(re_size_neg, 2);
        
        no_neg = numel(neg_idx);
        
        for vm_no = 1 : no_neg
            if (no_neg - vm_no + 1) * re_size_neg(vm_no,2) > sum_pos
                allocation(re_size_neg(vm_no:end,1)) = allocation(re_size_neg(vm_no:end,1)) + sum_pos / (no_neg - vm_no + 1);
                sum_pos = 0;
                break;
            else
                allocation(re_size_neg(vm_no:end,1)) = allocation(re_size_neg(vm_no:end,1)) + re_size_neg(vm_no, 2);
                re_size_neg(vm_no:end,2) = re_size_neg(vm_no:end,2) - re_size_neg(vm_no, 2);
                sum_pos = sum_pos - (no_neg - vm_no + 1) * re_size_neg(vm_no,2);
            end
        end      
    end
    
    under_provision = 0; prio_under_provision = 0;
    for vm_id = 1 : size_box - 1
        [~, idx] = min(true_sample{vm_id} > allocation(vm_id));    
        sum_ticket = sum_ticket + true_counts{vm_id}(idx); 
        under_provision = under_provision + sum(all_time_series(:, vm_id) == 100);
        prio_under_provision = prio_under_provision + sum(all_time_series_used_cap(:, vm_id) > allocation(vm_id));
    end
    PRIO_RESIZE_TICKET(end+1) = sum_ticket;
    UNDER_PROVISION(end+1, 1:2) = [under_provision, vm_id * numel(all_time_series_used_cap(:, vm_id))];
    PRIO_UNDER_PROVISION(end+1, 1:2) = [prio_under_provision, vm_id * numel(all_time_series_used_cap(:, vm_id))];
    
    % Second, calculate the dissimilarities among these time series using
    % cross-correlation. Need to notice that, here we change the xcf into
    % the range of [0,1], which makes it more sense to measure the
    % dissimilarity

    total_metric_num = size(all_time_series, 2);
    d = zeros(total_metric_num); d_half = zeros(total_metric_num);
    for row = 2 : total_metric_num
        for col = 1 : row -1    
            d(row, col) = dtw(all_time_series(:, row), ...
                              all_time_series(:, col), max_lag);
            
            d(col, row) = d(row, col);
            
            d_half(row, col) = d(row, col);  
            % in case of totally same case
            if d_half(row, col) == 0
                d_half(row, col) = -1;
            end
        end
    end    
           
    % Step 1: Do clustering and pick up the representative for each cluster 
    d_half_vec = reshape(d_half, 1, total_metric_num^2);
    d_vec = nonzeros(d_half_vec)';
    % revert to '0' if we have '-1'
    d_vec(find(d_vec == -1)) = 0;
    z = linkage(d_vec, 'single'); % hierarchical clustering
    
    % Create clusters and pick up the optimal number of clusters
    s_diff_cluster = [];
    for num_cluster = 2 : max(2, ceil(sqrt(total_metric_num)))
        cidx = cluster(z, 'maxclust', num_cluster);
        s_diff_cluster(num_cluster-1) = mean(silhouette([], cidx, d_vec));
    end
    
    [max_s, idx] = max(s_diff_cluster);
    
    % Do the optimal clustering 
    cidx = cluster(z, 'maxclust', idx + 1);
    
    representative_cluster = [];
    representative_cluster_time_series = [];
    cluster_time_series_idx = {};
    cluster_no = 1;
    while cluster_no <= idx + 1
        cluster_time_series_idx{cluster_no} = [];
        for vm_no = 1 : total_metric_num
            if cidx(vm_no) ~= cluster_no
                continue;
            end
            cluster_time_series_idx{cluster_no}(end+1) = vm_no;
        end
        
        % find the representative, who has the minimum distance with others
        [val, min_idx] = min(sum(d(:, cluster_time_series_idx{cluster_no})));
        
        vm_no = cluster_time_series_idx{cluster_no}(min_idx);
        
        representative_cluster(end+1) = vm_no;
        representative_cluster_time_series(:,end+1) = all_time_series(:, vm_no);

        cluster_no = cluster_no + 1;
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
    sample = {}; counts = {};
    for cand_id = 1 : numel(representative_cluster)
        X_fit = X;
        % if the time series is in the representative
        vm_no = representative_cluster(cand_id);
        if numel(find(final_representative_cluster == vm_no)) == 1
            box_vm_time_series_demands = all_time_series_used_cap(:, vm_no) / ticket_thres * 100;   
            sample{vm_no} = unique(box_vm_time_series_demands);
            if sample{vm_no}(1) ~= 0
                sample{vm_no} = [0; sample{vm_no}];
            end
            counts{vm_no} = hist(box_vm_time_series_demands, sample{vm_no});
            sample{vm_no} = flipud(sample{vm_no});
            counts{vm_no} = cumsum(fliplr(counts{vm_no})) - fliplr(counts{vm_no});
            
            inter_idx = find(cluster_time_series_idx{cand_id} == vm_no);
            cluster_time_series_idx{cand_id} = [cluster_time_series_idx{cand_id}(1:inter_idx-1), ...
                                                cluster_time_series_idx{cand_id}(inter_idx+1:end)];
            if SEPERATE_VALUE
                X_fit = [ones(time_len,1), all_time_series(:, vm_no)];
            end
            
            % we need to consider upper bound and lower bound
            % Add the lower bound to the case, which in case that the new
            % allocation will not furthur degrade the performance, especially
            % for the *truncated* case.
            lower_bound = max(all_time_series_used_cap(:, vm_no));
            [~, idx_lower] = min(sample{vm_no} >= lower_bound);

            sample{vm_no} = [sample{vm_no}(1 : idx_lower-1); lower_bound];
            counts{vm_no} = counts{vm_no}(1 : idx_lower);      

            % Due to allocation to each VM cannot exceed the total capacity in
            % BOX, so the first step is to change *true_sample* and
            % *true_counts*
            % here we could assume *box_available_capacity* stays same due to
            % the really small error in box capacity prediction.
            [~, idx] = min(sample{vm_no} > box_available_capacity);
            if idx ~= 1
                sample{vm_no} = sample{vm_no}(idx:end);
                counts{vm_no} = counts{vm_no}(idx:end);
            end
        end       
        
        for cand_in_cluster = 1 : numel(cluster_time_series_idx{cand_id})
            vm_no = cluster_time_series_idx{cand_id}(cand_in_cluster);
            Y = all_time_series(:, vm_no);
            [b,bint,r,rint,stats] = regress(Y,X_fit);

            y_fit = max(0, Y - r);
            box_vm_time_series_demands = original_allocation_each_vm(vm_no)/ticket_thres * y_fit;
            sample{vm_no} = unique(box_vm_time_series_demands);
            if sample{vm_no}(1) ~= 0
                sample{vm_no} = [0; sample{vm_no}];
            end
            counts{vm_no} = hist(box_vm_time_series_demands, sample{vm_no});
            sample{vm_no} = flipud(sample{vm_no});
            counts{vm_no} = cumsum(fliplr(counts{vm_no})) - fliplr(counts{vm_no});

            % we need to consider upper bound and lower bound
            % Add the lower bound to the case, which in case that the new
            % allocation will not furthur degrade the performance, especially
            % for the *truncated* case.
            lower_bound = max(original_allocation_each_vm(vm_no) * y_fit / 100);
            [~, idx_lower] = min(sample{vm_no} >= lower_bound);

            sample{vm_no} = [sample{vm_no}(1 : idx_lower-1); lower_bound];
            counts{vm_no} = counts{vm_no}(1 : idx_lower);      

            % Due to allocation to each VM cannot exceed the total capacity in
            % BOX, so the first step is to change *true_sample* and
            % *true_counts*
            % here we could assume *box_available_capacity* stays same due to
            % the really small error in box capacity prediction.
            [~, idx] = min(sample{vm_no} > box_available_capacity);
            if idx ~= 1
                sample{vm_no} = sample{vm_no}(idx:end);
                counts{vm_no} = counts{vm_no}(idx:end);
            end
            
            abs_error(end+1) = mean(abs(r));
            ape(end+1) = nanmean(abs(r) ./ Y);
            r_square(end+1) = stats(1); 
        end
  
    end
    
    % Step 5: fit box time series using the representative
    [b, bint, r, rint, stats] = regress(box_time_series, X);
    box_abs_error = mean(abs(r));
    box_ape = nanmean(abs(r) ./ box_time_series);
    box_r_square = stats(1);
    
    % Step 6: Ticket reduction results with fitting
    [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(sample, counts, total_available_capacity);
    if ~solution_or_not
        disp('It may happen');
    end
    sum_ticket = 0; each_left_for_vm = left_capacity / (size_box-1);
    exceed_allocation = [];
    allocation = [];
    for vm_id = 1 : size_box - 1
        allocation(vm_id) = sample{vm_id}(candidate_greedy(vm_id)) + each_left_for_vm;
        exceed_allocation(vm_id) = allocation(vm_id) - box_available_capacity;
    end
    
    if sum(exceed_allocation) >= 0 
        % This means each should at least have the *box available capacity*
        % however, this case should be very few
        when_all_resize_reach_box_cap = when_all_resize_reach_box_cap + 1;
        allocation = ones(size_box - 1) * box_available_capacity;
    else
        pos_idx = exceed_allocation >= 0;
        sum_pos = sum(exceed_allocation(pos_idx));
        allocation(pos_idx) = box_available_capacity;
        
        neg_idx = find(exceed_allocation < 0);
        neg_val = exceed_allocation(neg_idx);
        
        re_size_neg = [neg_idx', abs(neg_val)'];
        re_size_neg = sortrows(re_size_neg, 2);
        
        no_neg = numel(neg_idx);
        
        for vm_no = 1 : no_neg
            if (no_neg - vm_no + 1) * re_size_neg(vm_no,2) > sum_pos
                allocation(re_size_neg(vm_no:end,1)) = allocation(re_size_neg(vm_no:end,1)) + sum_pos / (no_neg - vm_no + 1);
                sum_pos = 0;
                break;
            else
                allocation(re_size_neg(vm_no:end,1)) = allocation(re_size_neg(vm_no:end,1)) + re_size_neg(vm_no, 2);
                re_size_neg(vm_no:end,2) = re_size_neg(vm_no:end,2) - re_size_neg(vm_no, 2);
                sum_pos = sum_pos - (no_neg - vm_no + 1) * re_size_neg(vm_no,2);
            end
        end      
    end
    fit_under_provision = 0;
    for vm_id = 1 : size_box - 1
        [~, idx] = min(true_sample{vm_id} >= allocation(vm_id));     
        sum_ticket = sum_ticket + true_counts{vm_id}(idx);   
        fit_under_provision = fit_under_provision + sum(all_time_series_used_cap(:, vm_id) > allocation(vm_id));
        if sum(all_time_series_used_cap(:, vm_id) > allocation(vm_id)) > 0
            FIT_UNDER_PROVISION_NO_VM = FIT_UNDER_PROVISION_NO_VM + 1;
        end
    end
    TOTAL_VM_NO = TOTAL_VM_NO + size_box -1;
    RESIZE_TICKET(end+1) = sum_ticket;
    FIT_UNDER_PROVISION(end+1, 1:2) = [fit_under_provision, vm_id * numel(all_time_series_used_cap(:, vm_id))];
     
    % write the summary of this box information after linear fitting
    BOX_ID(end+1) = pm_id; 
    ORIGINAL_VM_NUM(end+1) = total_metric_num; 
    REDUCED_VM_NO(end+1) = numel(final_representative_cluster);
    
    ALL_ABS_ERROR(:, end+1) = [box_abs_error; mean(abs_error); median(abs_error); prctile(abs_error,percentile_to_check)]; 
    ALL_APE(:,end+1) = [box_ape; mean(ape); median(ape); prctile(ape,percentile_to_check)]; 
    ALL_R_SQUARE(:,end+1) = [box_r_square; mean(r_square); median(r_square); prctile(r_square,percentile_to_check)];
    
end

% save(strcat(path, 'BOX_VM_APE_DWT_STEPWISE'), 'ALL_APE');
% save(strcat(path, 'BOX_VM_R_SQUARE_DWT_STEPWISE'), 'ALL_R_SQUARE');

VM_REDUCED_PCT = (ORIGINAL_VM_NUM-REDUCED_VM_NO) ./ ORIGINAL_VM_NUM;
% save(strcat(path, 'VM_REDUCED_PCT_DWT_STEPWISE'), 'VM_REDUCED_PCT');

not_zero_idx = ORIGINAL_TICKET > 0;
REDUCED_PCT = (ORIGINAL_TICKET(not_zero_idx) - RESIZE_TICKET(not_zero_idx))...
               ./ ORIGINAL_TICKET(not_zero_idx);
% save(strcat(path, 'REDUCED_TICKET_PCT'), 'REDUCED_PCT');

% Box plot: vm reduction

vm_reduced_pct = VM_REDUCED_PCT * 100;
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
ecdf(vm_reduced_pct);
xlabel('Number of VM Predictor Reduction (%)');
title(strcat('DTW Cluster, mean = ', {}, mat2str(int32(nanmean(vm_reduced_pct))), '%'));
set(gca, 'fontsize',15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'new_vm_reduction_cdf'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
ecdf(REDUCED_PCT*100);
xlabel('Number of Ticket Reduction (%)');
title(strcat('DTW Cluster, mean = ', {}, mat2str(int32(nanmean(REDUCED_PCT))), '%'));
set(gca, 'fontsize',15);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'new_ticket_reduction_cdf'));


