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

mkdir('../New_Data_time_series_linear_subset_vif_different_cost_budget_figure');
path = '../New_Data_time_series_linear_subset_vif_different_cost_budget_figure/';

ticket_thres_cand = [40, 50, 60];
cost_reduction_cand = [0, 10, 20, 30, 40, 50];
% 
% when_all_resize_reach_box_cap = 0;
% total_box_no = 0;

ORIGINAL_TICKET = {};

PRIO_RESIZE_TICKET = {}; 

FIT_RESIZE_TICKET = {};

BOX_SIZE = [];

for ticket_id = 1 : numel(ticket_thres_cand)
    ticket_thres = ticket_thres_cand(ticket_id);
    
    ORIGINAL_TICKET_SUMMARY = []; ORIGINAL_TICKET_SIG = [];
    PRIO_RESIZE_TICKET_SUMMARY = []; PRIO_RESIZE_TICKET_SIG = [];
    FIT_RESIZE_TICKET_SUMMARY = []; FIT_RESIZE_TICKET_SIG = [];
    
    for cost_id = 1 : numel(cost_reduction_cand)
        cost_thres = 100 - cost_reduction_cand(cost_id);
        
        ORIGINAL_TICKET{ticket_id, cost_id} = [];
        PRIO_RESIZE_TICKET{ticket_id, cost_id} = [];
        FIT_RESIZE_TICKET{ticket_id, cost_id} = [];
        
        for box_id = 1 : size_box_vm(2)

            size_box = numel(box_vm_time_series_summary{1, box_id});

            % If we don't have time series
            if size_box <= 2
                continue;
            end
        %     
%             if box_id > 100
%                 break;
%             end
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

            box_available_capacity = mean(box_vm_time_series_summary{1, box_id}{1,1}(:,3) ./ ...
                                          box_vm_time_series_summary{1, box_id}{1,1}(:,4) * 100);
            total_available_capacity = 0; 
            true_sample = {}; true_counts = {};
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
            ORIGINAL_TICKET{ticket_id, cost_id}(end+1) = original_ticket;
            ORIGINAL_TICKET_SUMMARY(end+1, 1:2) = [original_ticket, 100 - cost_thres];

            % PRIOR-KNOWLEDGE TICKET REDUCTION RESULTS
            [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(true_sample, true_counts, total_available_capacity * cost_thres / 100);
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
%                 when_all_resize_reach_box_cap = when_all_resize_reach_box_cap + 1;
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

            for vm_id = 1 : size_box - 1
                [~, idx] = min(true_sample{vm_id} > allocation(vm_id));    
                sum_ticket = sum_ticket + true_counts{vm_id}(candidate_greedy(vm_id));     
            end
            PRIO_RESIZE_TICKET{ticket_id, cost_id}(end+1) = sum_ticket;
            PRIO_RESIZE_TICKET_SUMMARY(end+1) = sum_ticket;

            % Second, calculate the dissimilarities among these time series using
            % cross-correlation. Need to notice that, here we change the xcf into
            % the range of [0,1], which makes it more sense to measure the
            % dissimilarity

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

            labels = {};
            for vm_id = 1 : size_box
                vm_label_no = {strcat(mat2str(vm_id))};
                labels = {labels{:}, vm_label_no{:}};
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
                    box_vm_time_series_demands = original_allocation_each_vm(vm_no)/ticket_thres * y_fit;
                    sample{vm_no} = unique(box_vm_time_series_demands);
                    if sample{vm_no}(1) ~= 0
                        sample{vm_no} = [0; sample{vm_no}];
                    end
                    counts{vm_no} = hist(box_vm_time_series_demands, sample{vm_no});
                    sample{vm_no} = flipud(sample{vm_no});
                    counts{vm_no} = cumsum(fliplr(counts{vm_no})) - fliplr(counts{vm_no});

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
            [candidate_greedy, solution_or_not, left_capacity] = greedy_find_approximate(sample, counts, total_available_capacity * cost_thres / 100);
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
            for vm_id = 1 : size_box - 1
                [~, idx] = min(true_sample{vm_id} >= allocation(vm_id));     
                sum_ticket = sum_ticket + true_counts{vm_id}(idx);     
            end
            FIT_RESIZE_TICKET{ticket_id, cost_id}(end+1) = sum_ticket;
            FIT_RESIZE_TICKET_SUMMARY(end+1) = sum_ticket;

            % write the summary of this box information after linear fitting
            if ticket_id * cost_id == 1
                BOX_SIZE(end+1) = total_metric_num; 
            end
            
        end
        
        % For a certain COST REDUCTION CASE
        PRIO_REDUCTION = (ORIGINAL_TICKET{ticket_id,cost_id} - PRIO_RESIZE_TICKET{ticket_id,cost_id}) ./ ORIGINAL_TICKET{ticket_id,cost_id};
        PRIO_RESIZE_TICKET_SIG(end+1, 1:3) = [nanmean(PRIO_REDUCTION), ...
                                              nanmedian(PRIO_REDUCTION),...
                                              prctile(PRIO_REDUCTION,90)];
                                          
        FIT_REDUCTION = (ORIGINAL_TICKET{ticket_id,cost_id} - FIT_RESIZE_TICKET{ticket_id,cost_id}) ./ ORIGINAL_TICKET{ticket_id,cost_id};                                  
        FIT_RESIZE_TICKET_SIG(end+1, 1:3) = [nanmean(FIT_REDUCTION), ...
                                              nanmedian(FIT_REDUCTION),...
                                              prctile(FIT_REDUCTION,90)];
    end
    
    save(strcat(path, 'PRIO_RESIZE_TICKET_SUMMARY_', mat2str(ticket_thres)), 'PRIO_RESIZE_TICKET_SUMMARY');
    save(strcat(path, 'FIT_REDUCED_TICKET_SUMMARY_', mat2str(ticket_thres)), 'FIT_RESIZE_TICKET_SUMMARY');
    save(strcat(path, 'ORIGINAL_TICKET_SUMMARY_', mat2str(ticket_thres)), 'ORIGINAL_TICKET_SUMMARY');

    save(strcat(path, 'PRIO_RESIZE_TICKET_SIG_', mat2str(ticket_thres)), 'PRIO_RESIZE_TICKET_SIG');
    save(strcat(path, 'FIT_REDUCED_TICKET_SIG_', mat2str(ticket_thres)), 'FIT_RESIZE_TICKET_SIG');
    save(strcat(path, 'ORIGINAL_TICKET_SIG_', mat2str(ticket_thres)), 'ORIGINAL_TICKET_SIG');

    PRIO_REDUCTION = [(ORIGINAL_TICKET_SUMMARY(:,1) - PRIO_RESIZE_TICKET_SUMMARY') ./ ...
                       ORIGINAL_TICKET_SUMMARY(:,1), ORIGINAL_TICKET_SUMMARY(:,2)];

    FIT_REDUCTION = [(ORIGINAL_TICKET_SUMMARY(:,1) - FIT_RESIZE_TICKET_SUMMARY') ./ ...
                       ORIGINAL_TICKET_SUMMARY(:,1), ORIGINAL_TICKET_SUMMARY(:,2)];

    fig = figure;
    set(fig, 'Position', [200, 200, 600, 400]);
    boxplot(PRIO_REDUCTION(:,1) * 100, PRIO_REDUCTION(:,2));
    xlabel('Cost Reduction (%)');
    ylabel('Ticket Reduction (%)');
    set(gca, 'ylim', [-100 100]); set(gca, 'ytick', [-100:20:100]);
    title(strcat('Alg 1: Ticket Threshold(%) is', {' '}, mat2str(ticket_thres)));
    set(gca, 'fontsize', 15);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path,'box_prio_resize_', mat2str(ticket_thres)));
    
    fig = figure;
    set(fig, 'Position', [200, 200, 600, 400]);
    boxplot(FIT_REDUCTION(:,1) * 100, FIT_REDUCTION(:,2));
    ylabel('Ticket Reduction (%)');
    xlabel('Cost Reduction (%)');
    set(gca, 'ylim', [-100 100]); set(gca, 'ytick', [-100:20:100]);
    title(strcat('Alg 2: Ticket Threshold(%) is', {' '},  mat2str(ticket_thres)));
    set(gca, 'fontsize', 15);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path,'box_fit_resize_', mat2str(ticket_thres)));
end

save(strcat(path, 'PRIO_RESIZE_TICKET'), 'PRIO_RESIZE_TICKET');
save(strcat(path, 'FIT_REDUCED_TICKET'), 'FIT_RESIZE_TICKET');
save(strcat(path, 'ORIGINAL_TICKET'), 'ORIGINAL_TICKET');
save(strcat(path, 'BOX_SIZE'), 'BOX_SIZE');



            

