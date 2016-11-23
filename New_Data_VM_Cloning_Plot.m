% Compare the VM cloning simulation for w/ and w/o tenant clustering

close all; clear; clc

moments = 5 : 5 : 95;
moments = [moments, 98];

server_num_thres = 10;
ticket_thres = 60;

check_usage_tail = [numel(moments) - 2, numel(moments) - 1];



fig_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Clone_with_or_without_Clustering_Fig/'); 
mkdir(fig_path);

font_size = 13;

first_try = false;
if first_try
    for tail_idx = 1 : numel(check_usage_tail)
        usage_tail_idx = check_usage_tail(tail_idx);

        result_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Target_Results/'); 
        load (strcat(result_path, 'USED_TENANT_CLUSTER_', moments(usage_tail_idx)));
        USED_TENANT_CLUSTER_ORIGINAL = USED_TENANT_CLUSTER;
        load (strcat(result_path, 'BOX_VIOLATION_RATIO_TRAIN_', moments(usage_tail_idx)));
        BOX_VIOLATION_RATIO_TRAIN_ORIGINAL = BOX_VIOLATION_RATIO_TRAIN;
        load (strcat(result_path, 'tail_violation_', moments(usage_tail_idx)));
        tail_violation_original = tail_violation;
        load (strcat(result_path, 'ticket_compare_', moments(usage_tail_idx)));
        ticket_compare_original = ticket_compare;

        result_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Target_Cluster_Added_Results/'); 
        load (strcat(result_path, 'USED_TENANT_CLUSTER_', moments(usage_tail_idx)));
        load (strcat(result_path, 'BOX_VIOLATION_RATIO_TRAIN_', moments(usage_tail_idx)));
        load (strcat(result_path, 'tail_violation_', moments(usage_tail_idx)));
        load (strcat(result_path, 'ticket_compare_', moments(usage_tail_idx)));

        % first unite two cases
        [all_tenant_id, ia, ib] = intersect(USED_TENANT_CLUSTER_ORIGINAL, USED_TENANT_CLUSTER(:,1)); 
        box_violation_original_train = BOX_VIOLATION_RATIO_TRAIN_ORIGINAL(ia);    
        tail_violation_original_new = tail_violation_original(ia, :);
        ticket_compare_original_new = ticket_compare_original(ia, :);

        cluster_new = USED_TENANT_CLUSTER(ib,2);

        box_violation_train = BOX_VIOLATION_RATIO_TRAIN(ib);
        tail_violation_new = tail_violation(ib, :);
        ticket_compare_new = ticket_compare(ib, :);

        % Different Cluster, Different Results
        disp('*********************************');
        disp('*********************************');
        disp('*********************************');
        disp('*********************************');
        disp('*********************************');

        ticket_reduction_original_new = (ticket_compare_original_new(:,1) - ticket_compare_original_new(:,2)) ...
                                         ./ ticket_compare_original_new(:,1) * 100;
        ticket_reduction_new = (ticket_compare_new(:,1) - ticket_compare_new(:,2)) ./ ticket_compare_new(:,1) * 100;

        for cluster_id = 1 : 3
            num_cluster = find(cluster_new == cluster_id);
            disp('*********************************');
            disp(strcat('num of cluster is ', {' '}, mat2str(numel(num_cluster))));
            disp(strcat('original - mean tail violation is ', {' '}, mat2str(mean(tail_violation_original_new(num_cluster, :)))));
            disp(strcat('original - mean ticket reduction is ', {' '}, mat2str(nanmean(ticket_reduction_original_new(num_cluster)))));

            disp(strcat('new - mean tail violation is ', {' '}, mat2str(mean(tail_violation_new(num_cluster, :)))));
            disp(strcat('new - mean ticket reduction is ', {' '}, mat2str(nanmean(ticket_reduction_new(num_cluster)))));
        end

        % Figure 1: Tail Target Violation
        fig = figure;
        set(fig, 'Position', [200 200 600 300]);
        num_case = numel(all_tenant_id);
        [train_violation, idx] = sort(box_violation_original_train);
        plot(1 : num_case, tail_violation_original_new(idx, 1) * 100, 'b*:', 'linewidth', 1.5);
        hold on
        plot(1 : num_case, tail_violation_original_new(idx, 2) * 100, 'r^--', 'linewidth', 1.5);
        hold on
        plot(1 : num_case, tail_violation_new(idx, 2) * 100, 'ko-.', 'linewidth', 1.5);
        xlabel('Ratio of Box w/ Tail Violation in Training Stage (%)'); 
        ylabel('Ratio of Boxes w/ Tail Violation(%)');
        title1 = strcat('Before: Mean Violation = ', mat2str(round(mean(tail_violation_original_new(:,1)) * 100, 2)), '%');
        title2 = strcat('After w/o Clutering: Mean Violation = ', mat2str(round(mean(tail_violation_original_new(:,2)) * 100, 2)), '%');
        title3 = strcat('After w/ Clutering: Mean Violation = ', mat2str(round(mean(tail_violation_new(:,2)) * 100, 2)), '%');
        set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:1:num_case]); 
        set(gca, 'xticklabel', round(train_violation * 100));
        xticklabel_rotate([], 45);
        set(gca, 'ylim', [0 80]); set(gca, 'ytick', [0:10:80]); 
        h = legend(title1, title2, title3);
        set(h, 'Location', 'northwest');
        set(gca, 'fontsize', font_size);

        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(fig_path, 'tail_violation_compare_', mat2str(moments(usage_tail_idx))));

        % Figure 2: Ticket Violation
        non_zero_idx = find(ticket_compare_new(:,1) ~= 0);
        ticket_compare_new = ticket_compare_new(non_zero_idx, :);
        box_violation_train = box_violation_train(non_zero_idx);

        ticket_compare_original_new = ticket_compare_original_new(non_zero_idx, :);
        box_violation_original_train = box_violation_original_train(non_zero_idx);

        fig = figure;
        set(fig, 'Position', [200 200 600 300]);
        [ticket_tenant, idx] = sortrows([box_violation_train', round(ticket_compare_new(:,1))], [-1, -2]);
        num_case = numel(idx);
        stem3(1:num_case, ticket_tenant(:,2), ...
              (ticket_compare_original_new(idx,1) - ticket_compare_original_new(idx,2)) ./ ticket_compare_original_new(idx,1) * 100, ...
              'r^');
        hold on
        stem3(1:num_case, ticket_tenant(:,2), ...
              (ticket_compare_new(idx,1) - ticket_compare_new(idx,2)) ./ ticket_compare_new(idx,1) * 100, ...
              'ko');

        xlabel('Ratio of Box w/ Tail Violation(%)'); ylabel('Original Ticket Number')
        zlabel('Ticket Reduction (%) per Box');
        set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:2:num_case]); 
        set(gca, 'xticklabel', round(ticket_tenant(1:2:end, 1) * 100));
        upper = ceil(max(ticket_compare(:,1)) / 10) * 10;
        set(gca, 'ylim', [0 upper]); set(gca, 'ytick', [0: upper/5: upper]); 
        set(gca, 'zlim', [0 100]); set(gca, 'ztick', [0 : 20 : 100]);
        title1 = strcat('After w/o Clutering: Mean Ticket Reduction = ', mat2str(round(...
                     nanmean((ticket_compare_original_new(idx,1) - ticket_compare_original_new(idx,2)) ./ ticket_compare_original_new(idx,1) * 100), ...
                     2)), '%');
        title2 = strcat('After w/ Clutering: Mean Ticket Reduction = ', mat2str(round(...
                     nanmean((ticket_compare_new(idx,1) - ticket_compare_new(idx,2)) ./ ticket_compare_new(idx,1) * 100), ...
                     2)), '%');
        h = legend(title1, title2);
        set(h, 'location', 'northoutside');
        set(gca, 'fontsize', font_size); 
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(fig_path, 'ticket_reduction_', mat2str(moments(usage_tail_idx))));

    end   
end


fig_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Clone_Optimal_Fig/'); 
mkdir(fig_path);

second_try = true;
if second_try
    for tail_idx = 1 : 1%numel(check_usage_tail)
        usage_tail_idx = check_usage_tail(tail_idx);

        result_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Target_OptimalPlacement_Results/'); 
        load (strcat(result_path, 'USED_TENANT_CLUSTER_', moments(usage_tail_idx)));
        USED_TENANT_CLUSTER_ORIGINAL = USED_TENANT_CLUSTER;
        load (strcat(result_path, 'BOX_VIOLATION_RATIO_TRAIN_', moments(usage_tail_idx)));
        BOX_VIOLATION_RATIO_TRAIN_ORIGINAL = BOX_VIOLATION_RATIO_TRAIN;
        load (strcat(result_path, 'tail_violation_', moments(usage_tail_idx)));
        tail_violation_original = tail_violation;
        load (strcat(result_path, 'ticket_compare_', moments(usage_tail_idx)));
        ticket_compare_original = ticket_compare;

        result_path = strcat('../New_Data_box_moments_fit/', mat2str(ticket_thres),'_Target_Cluster_Added_Results/'); 
        load (strcat(result_path, 'USED_TENANT_CLUSTER_', moments(usage_tail_idx)));
        load (strcat(result_path, 'BOX_VIOLATION_RATIO_TRAIN_', moments(usage_tail_idx)));
        load (strcat(result_path, 'tail_violation_', moments(usage_tail_idx)));
        load (strcat(result_path, 'ticket_compare_', moments(usage_tail_idx)));

        % first unite two cases
        [all_tenant_id, ia, ib] = intersect(USED_TENANT_CLUSTER_ORIGINAL, USED_TENANT_CLUSTER(:,1)); 
        box_violation_original_train = BOX_VIOLATION_RATIO_TRAIN_ORIGINAL(ia);    
        tail_violation_original_new = tail_violation_original(ia, :);
        ticket_compare_original_new = ticket_compare_original(ia, :);

        cluster_new = USED_TENANT_CLUSTER(ib,2);

        box_violation_train = BOX_VIOLATION_RATIO_TRAIN(ib);
        tail_violation_new = tail_violation(ib, :);
        ticket_compare_new = ticket_compare(ib, :);

        % Different Cluster, Different Results
        disp('*********************************');
        disp('*********************************');
        disp('*********************************');
        disp('*********************************');
        disp('*********************************');

        ticket_reduction_original_new = (ticket_compare_original_new(:,1) - ticket_compare_original_new(:,2)) ...
                                         ./ ticket_compare_original_new(:,1) * 100;
        ticket_reduction_new = (ticket_compare_new(:,1) - ticket_compare_new(:,2)) ./ ticket_compare_new(:,1) * 100;

        for cluster_id = 1 : 3
            num_cluster = find(cluster_new == cluster_id);
            if numel(num_cluster) > 1
                disp('*********************************');
                disp(strcat('num of cluster is ', {' '}, mat2str(numel(num_cluster))));
                disp(strcat('optimal - mean tail violation is ', {' '}, mat2str(mean(tail_violation_original_new(num_cluster, :)))));
                disp(strcat('optimal - mean ticket reduction is ', {' '}, mat2str(nanmean(ticket_reduction_original_new(num_cluster)))));

                disp(strcat('new - mean tail violation is ', {' '}, mat2str(mean(tail_violation_new(num_cluster, :)))));
                disp(strcat('new - mean ticket reduction is ', {' '}, mat2str(nanmean(ticket_reduction_new(num_cluster)))));
            else
                disp('*********************************');
                disp(strcat('num of cluster is ', {' '}, mat2str(numel(num_cluster))));
                disp(strcat('optimal - mean tail violation is ', {' '}, mat2str((tail_violation_original_new(num_cluster, :)))));
                disp(strcat('optimal - mean ticket reduction is ', {' '}, mat2str((ticket_reduction_original_new(num_cluster)))));

                disp(strcat('new - mean tail violation is ', {' '}, mat2str((tail_violation_new(num_cluster, :)))));
                disp(strcat('new - mean ticket reduction is ', {' '}, mat2str((ticket_reduction_new(num_cluster)))));
            end
        end

        % Figure 1: Tail Target Violation
        fig = figure;
        set(fig, 'Position', [200 200 600 300]);
        num_case = numel(all_tenant_id);
        [train_violation, idx] = sort(box_violation_original_train);
        plot(1 : num_case, tail_violation_original_new(idx, 1) * 100, 'b*:', 'linewidth', 1.5);
        hold on
        plot(1 : num_case, tail_violation_original_new(idx, 2) * 100, 'r^--', 'linewidth', 1.5);
        hold on
        plot(1 : num_case, tail_violation_new(idx, 2) * 100, 'ko-.', 'linewidth', 1.5);
        xlabel('Ratio of Box w/ Tail Violation in Training Stage (%)'); 
        ylabel('Ratio of Boxes w/ Tail Violation(%)');
        title1 = strcat('Before: Mean Violation = ', mat2str(round(mean(tail_violation_original_new(:,1)) * 100, 2)), '%');
        title2 = strcat('Optimal: Mean Violation = ', mat2str(round(mean(tail_violation_original_new(:,2)) * 100, 2)), '%');
        title3 = strcat('After w/ Clutering: Mean Violation = ', mat2str(round(mean(tail_violation_new(:,2)) * 100, 2)), '%');
        set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:1:num_case]); 
        set(gca, 'xticklabel', round(train_violation * 100));
        xticklabel_rotate([], 45);
        set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0:10:100]); 
        h = legend(title1, title2, title3);
        set(h, 'Location', 'northwest');
        set(gca, 'fontsize', font_size);

        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(fig_path, 'tail_violation_compare_', mat2str(moments(usage_tail_idx))));

        % Figure 2: Ticket Violation
        non_zero_idx = find(ticket_compare_new(:,1) ~= 0);
        ticket_compare_new = ticket_compare_new(non_zero_idx, :);
        box_violation_train = box_violation_train(non_zero_idx);

        ticket_compare_original_new = ticket_compare_original_new(non_zero_idx, :);
        box_violation_original_train = box_violation_original_train(non_zero_idx);

        fig = figure;
        set(fig, 'Position', [200 200 600 300]);
        [ticket_tenant, idx] = sortrows([box_violation_train', round(ticket_compare_new(:,1))], [-1, -2]);
        num_case = numel(idx);
        stem3(1:num_case, ticket_tenant(:,2), ...
              (ticket_compare_original_new(idx,1) - ticket_compare_original_new(idx,2)) ./ ticket_compare_original_new(idx,1) * 100, ...
              'r^');
        hold on
        stem3(1:num_case, ticket_tenant(:,2), ...
              (ticket_compare_new(idx,1) - ticket_compare_new(idx,2)) ./ ticket_compare_new(idx,1) * 100, ...
              'ko');

        xlabel('Ratio of Box w/ Tail Violation(%)'); ylabel('Original Ticket Number')
        zlabel('Ticket Reduction (%) per Box');
        set(gca, 'xlim', [1 num_case]); set(gca, 'xtick', [1:2:num_case]); 
        set(gca, 'xticklabel', round(ticket_tenant(1:2:end, 1) * 100));
        upper = ceil(max(ticket_compare(:,1)) / 10) * 10;
        set(gca, 'ylim', [0 upper]); set(gca, 'ytick', [0: upper/5: upper]); 
        set(gca, 'zlim', [0 100]); set(gca, 'ztick', [0 : 20 : 100]);
        title1 = strcat('Optimal: Mean Ticket Reduction = ', mat2str(round(...
                     nanmean((ticket_compare_original_new(idx,1) - ticket_compare_original_new(idx,2)) ./ ticket_compare_original_new(idx,1) * 100), ...
                     2)), '%');
        title2 = strcat('After w/ Clutering: Mean Ticket Reduction = ', mat2str(round(...
                     nanmean((ticket_compare_new(idx,1) - ticket_compare_new(idx,2)) ./ ticket_compare_new(idx,1) * 100), ...
                     2)), '%');
        h = legend(title1, title2);
        set(h, 'location', 'northoutside');
        set(gca, 'fontsize', font_size); 
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(fig_path, 'ticket_reduction_', mat2str(moments(usage_tail_idx))));

    end
end