% Tenant Classification
% close all; clear; clc
% 
% load ../New_Data_7days/box_vm_time_series_summary_mem_only_with_zeros
% box_vm_time_series_mem = box_vm_time_series_summary;
% load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros
% 
% moment_path = '../New_Data_box_moments_fit/moment_coeff/';
% load (strcat(moment_path, 'unique_tenant'));
% load (strcat(moment_path, 'tenant_box_series'));

fig_path = strcat('../New_Data_box_moments_fit/','Tenant_Box_Dependency_Fig/');
mkdir(fig_path);

day_thres = 5;
wind_length = 96;

moments = 5 : 5 : 95;
moments = [moments, 98];

server_num_thres = 20;
ticket_thres = 60;
% usage_tail_idx = numel(moments) - 3;

% Zero Check Parameters:
cc = {'g', 'b', 'r', 'k' ,'m', 'y', 'c'};
line_style = {'--', ':', '-.', '-', ':'};
marker_style = {'*', '<', '^', 'o', '>', 'x', 's', '+', 'd', 'p', 'h', '.'};
BOX_USAGE_MEAN_TAIL_TENANT = {};

% First Check Parameters:
BOX_TENANT_USAGE_TAIL = {}; MEAN_BOX_TENANT_USAGE_TAIL = {};

% Third Check Parameters:
BOX_MEAN = []; BOX_MOMENTS = [];
MOMENT_MEAN_COEFF = []; MOMENT_STD_COEFF = []; MOMENT_CONF_COEFF = [];
MOMENT_MEAN = []; MOMENT_STD = [];

for usage_tail_idx = numel(moments) - 3 : numel(moments)
    BOX_TENANT_USAGE_TAIL{end+1} = [];
    MEAN_BOX_TENANT_USAGE_TAIL{end+1} = [];
    BOX_USAGE_MEAN_TAIL_TENANT{end+1} = [];
    
    usage_tail = moments(usage_tail_idx);
    all_usage = {}; mean_box_usage_tail = [];
    
    % Second Check Parameters:
    TENANT_BOX_USAGE_TAIL = []; TENANT_USAGE_MOMENTS =[];
    TENANT_USAGE_PCTILES = [];
    for tenant_id = 1 : numel(unique_tenant)
        num_box = numel(tenant_box_series{tenant_id});
        if  num_box < server_num_thres
            break;
        end

        %%%%%%%%%%%%%%%%%%%%%% Zero Index to check %%%%%%%%%%%%%%%%%%%%%%%%
        % Motivation: Why we need to cluster Tenant?
        % reverse engineering: Do we need to classify Tenant?
        % Pick up the first *5* biggest tenants w/ about 1000 Boxes
        if tenant_id <= 5
            time_stamp = []; 
            for box_id = 1 : num_box
                box_usage = tenant_box_series{tenant_id}{box_id}{1}(:,end-2); 
                BOX_USAGE_MEAN_TAIL_TENANT{end}(end+1, :) = [nanmean(box_usage), ...
                                       prctile(box_usage, usage_tail), tenant_id];
                                   
                % find the longest the time stamp
                new_time_stamp = tenant_box_series{tenant_id}{box_id}{1}(:,1);
                time_stamp = union(new_time_stamp, time_stamp);                                   
            end
            mean_box_usage_tail(end+1) = nanmean(BOX_USAGE_MEAN_TAIL_TENANT{end}(:,2));
            % disp the begin and end
            disp(strcat('time length is ', {' '}, mat2str(numel(time_stamp))));
            
            % make up these holes
            tenant_demands_series = zeros(numel(time_stamp), 1);
            tenant_capacity = 0;
            for box_id = 1 : num_box
                box_demands = tenant_box_series{tenant_id}{box_id}{1}(:,end-3);
                box_usage = tenant_box_series{tenant_id}{box_id}{1}(:,end-2);
                new_time_stamp = tenant_box_series{tenant_id}{box_id}{1}(:,1);
                
                % [~, idx] = setdiff(time_stamp, new_time_stamp);
                [~, ~, idx] = intersect(time_stamp, new_time_stamp); 
                tenant_demands_series(idx) = tenant_demands_series(idx) + box_demands;
                
                % non zero index
                non_zero_idx = find(box_usage ~= 0);
                tenant_capacity = tenant_capacity + nanmean(box_demands(non_zero_idx) ...
                                                    ./ box_usage(non_zero_idx) * 100);
            end           
            all_usage{end+1} = tenant_demands_series / tenant_capacity * 100;         
        end

        %%%%%%%%%%%%%%%%%%%%%% Find the same time stamp %%%%%%%%%%%%%%%%%%%
        common_time_series = tenant_box_series{tenant_id}{1}{1}(:,1);
        for box_id = 2 : num_box
            temp_time = tenant_box_series{tenant_id}{box_id}{1}(:,1);
            common_time_series = intersect(common_time_series, temp_time);
        end

        % we only pick the long enough data
%         if numel(common_time_series) < wind_length * day_thres
%             % disp(strcat('tenant id is ', mat2str(tenant_id), 'num of box is ', mat2str(num_box)));
%             continue;
%         end

        %%%%%%%%%%%%%%%%%%%%%% First Index to check %%%%%%%%%%%%%%%%%%%%%%%
        % the box usage tail *correlates* with tenant usage tail?
        box_usage_tail = []; box_usage_mean = [];
        tenant_demand_series = zeros(numel(common_time_series), 1); 
        tenant_capacity = 0;
        for box_id = 1 : num_box
            temp_time = tenant_box_series{tenant_id}{box_id}{1}(:,1);
            [~, idx, ~] = intersect(temp_time, common_time_series);
            tenant_demand_series = tenant_demand_series + ...
                                   tenant_box_series{tenant_id}{box_id}{1}(idx, end-3);
            tenant_capacity = tenant_capacity + ...
                              nanmean(tenant_box_series{tenant_id}{box_id}{1}(idx, end-3) ...
                              ./ tenant_box_series{tenant_id}{box_id}{1}(idx, end-2) * 100);

            box_usage_tail(box_id) = prctile(tenant_box_series{tenant_id}{box_id}{1}(idx, end-2), usage_tail); 
            box_usage_mean(box_id) = nanmean(tenant_box_series{tenant_id}{box_id}{1}(idx, end-2));
        end
        tenant_usage_series = tenant_demand_series / tenant_capacity * 100; 
        tenant_usage_tail = prctile(tenant_usage_series, usage_tail);

        BOX_TENANT_USAGE_TAIL{end}(end+1:end+num_box, :) = [box_usage_tail', ...
                        ones(num_box, 1) * tenant_usage_tail, box_usage_mean'];
        MEAN_BOX_TENANT_USAGE_TAIL{end}(end+1, :) = [nanmean(box_usage_tail), ...
                        tenant_usage_tail, nanmean(box_usage_mean), nanmean(tenant_usage_series)];   

        %%%%%%%%%%%%%%%%%%%%% second things to check %%%%%%%%%%%%%%%%%%%%%%
        % cluster tenants based on the box usage tail 
        second_try = false;
        if second_try
            grat = 10;
            tenant_sum_box_tail_usage = zeros(100/grat + 1, 2);
            time_stamp = [];
            for box_id = 1 : num_box
                box_usage = tenant_box_series{tenant_id}{box_id}{1}(:, end-2);
                mean_box_usage_idx = round(nanmean(box_usage) / grat) + 1;

                tenant_sum_box_tail_usage(mean_box_usage_idx,1) = ...
                    tenant_sum_box_tail_usage(mean_box_usage_idx,1) + prctile(box_usage, usage_tail);
                tenant_sum_box_tail_usage(mean_box_usage_idx,2) = ...
                    tenant_sum_box_tail_usage(mean_box_usage_idx,2) + 1;

                % find the longest the time stamp
                new_time_stamp = tenant_box_series{tenant_id}{box_id}{1}(:,1);
                time_stamp = union(new_time_stamp, time_stamp);      
            end
            tenant_mean_box_tail_usage = tenant_sum_box_tail_usage(:,1) ./ tenant_sum_box_tail_usage(:,2);
            tenant_mean_box_tail_usage(isnan(tenant_mean_box_tail_usage)) = 0;
            TENANT_BOX_USAGE_TAIL(:, end+1) = tenant_mean_box_tail_usage;

            % calculate: mean, cov and different prctiles
            tenant_demands_series = zeros(numel(time_stamp), 1);
            tenant_capacity = 0;
            for box_id = 1 : num_box
                box_demands = tenant_box_series{tenant_id}{box_id}{1}(:,end-3);
                box_usage = tenant_box_series{tenant_id}{box_id}{1}(:,end-2);
                new_time_stamp = tenant_box_series{tenant_id}{box_id}{1}(:,1);

                % [~, idx] = setdiff(time_stamp, new_time_stamp);
                [~, ~, idx] = intersect(time_stamp, new_time_stamp); 
                tenant_demands_series(idx) = tenant_demands_series(idx) + box_demands;

                % non zero index
                non_zero_idx = find(box_usage ~= 0);
                tenant_capacity = tenant_capacity + nanmean(box_demands(non_zero_idx) ...
                                                        ./ box_usage(non_zero_idx) * 100);
            end        
            tenant_usage = tenant_demands_series / tenant_capacity * 100;
            TENANT_USAGE_MOMENTS(end+1,:) = [nanmean(tenant_usage), std(tenant_usage) / nanmean(tenant_usage)];
            TENANT_USAGE_PCTILES(end+1, :) = prctile(tenant_usage, moments);
        end
  
    end

    
    %%%%%%%%%%%%%%%%%%%%% Zero Index to check %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reverse engineering: Do we need to classify Tenant?
    fig = figure;
    set(fig, 'Position', [200 200 500 250]);
    font_size = 15;
    tenant = {};
    for tenant_id = 3 : 4
        [f, x] = ecdf(all_usage{tenant_id});
        plot(x, f, 'linestyle', line_style{tenant_id}, ...
             'color', cc{tenant_id}, 'linewidth', 2);
        hold on  
        str = strcat('Tenant', mat2str(tenant_id-2), ': Mean =', {' '}, ...
                     mat2str(round(nanmean(all_usage{tenant_id}), 1)), ...
                     ', C.O.V =', {' '}, ...
                     mat2str(round(std(all_usage{tenant_id}) / nanmean(all_usage{tenant_id}), 2)));
    end
    xlabel('Tenant CPU Usage');
    ylabel('CDF');
    % set(gca, 'xlim', [0 80], 'xtick', [0 : 10 : 80]);
    set(gca, 'ylim', [0 1], 'ytick', [0 : 0.1 : 1]);
    tenant1 = strcat('Tenant1: Mean =', {' '}, ...
                 mat2str(round(nanmean(all_usage{3}), 1)), ...
                 ', C.O.V =', {' '}, ...
                 mat2str(round(std(all_usage{3}) / nanmean(all_usage{3}), 2)));
    tenant2 = strcat('Tenant2: Mean =', {' '}, ...
                 mat2str(round(nanmean(all_usage{4}), 1)), ...
                 ', C.O.V =', {' '}, ...
                 mat2str(round(std(all_usage{4}) / nanmean(all_usage{4}), 2)));         

    h = legend(char(tenant1), char(tenant2));
    set(h, 'location', 'southeast', 'box', 'off');
    set(gca, 'fontsize', font_size);  
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'tenant_box_usage'));
    
    % PDF after discretization
    fig = figure;
    set(fig, 'Position', [200 200 600 500]);
    font_size = 15;
    subplot(2,2,1)
    for tenant_id = 3 : 4
        choose_id = find(BOX_USAGE_MEAN_TAIL_TENANT{end}(:, end) == tenant_id);
        scatter(BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 1), ...
                BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 2), ...
                marker_style{tenant_id}, cc{tenant_id});
        hold on   
    end
    xlabel('Mean Box Usage');
    ylabel(strcat('Box Usage Tail =', {' '}, mat2str(usage_tail), '%' ));
    set(gca, 'xlim', [0 80], 'xtick', [0 : 10 : 80]);
    set(gca, 'ylim', [0 80], 'ytick', [0 : 10 : 80]);
    h = legend('Tenant1 Size: 152', 'Tenant2 Size: 100');
    set(h, 'location', 'southeast', 'box', 'off');
    title('Original')
    set(gca, 'fontsize', font_size); 
    
    subplot(2,2,2)
    for tenant_id = 3 : 4
        choose_id = find(BOX_USAGE_MEAN_TAIL_TENANT{end}(:, end) == tenant_id);
        scatter(round(BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 1)/10)*10, ...
                BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 2), ...
                marker_style{tenant_id}, cc{tenant_id});
        hold on   
    end
    xlabel('Mean Box Usage');
    ylabel(strcat('Box Usage Tail =', {' '}, mat2str(usage_tail), '%' ));
    set(gca, 'xlim', [0 80], 'xtick', [0 : 10 : 80]);
    set(gca, 'ylim', [0 80], 'ytick', [0 : 10 : 80]);
    h = legend('Tenant1 Size: 152', 'Tenant2 Size: 100');
    set(h, 'location', 'southeast', 'box', 'off');
    title('Discretization')
    set(gca, 'fontsize', font_size);
    
    subplot(2,2,3)
    pdf = [];
    for tenant_id = 3 : 4
        choose_id = find(BOX_USAGE_MEAN_TAIL_TENANT{end}(:, end) == tenant_id);
        mean_usage = round(BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 1)/10)*10;
        tail_usage =BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 2);
        
        choose_box_idx = find(mean_usage == 10);
        bin = 10 : 5 : 40;
        [N, edges] = histcounts(tail_usage(choose_box_idx), bin);
        pdf(:, end+1) = N' / sum(N);
    end
    h = bar(pdf, 'grouped');
    set(h(1), 'facecolor', cc{3}); 
    set(h(2), 'facecolor', cc{4}); 
    set(gca, 'xticklabel', edges(1:end-1));
    xlabel('Box Usage Tail'); ylabel('PDF');
    h = legend('Tenant1 Size: 152', 'Tenant2 Size: 100');
    set(h, 'location', 'northeast', 'box', 'off');
    title('Mean Usage = 10')
    set(gca, 'fontsize', font_size);
    
    subplot(2,2,4)
    pdf = [];
    for tenant_id = 3 : 4
        choose_id = find(BOX_USAGE_MEAN_TAIL_TENANT{end}(:, end) == tenant_id);
        mean_usage = round(BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 1)/10)*10;
        tail_usage =BOX_USAGE_MEAN_TAIL_TENANT{end}(choose_id, 2);
        
        choose_box_idx = find(mean_usage == 20);
        bin = 10 : 5 : 40;
        [N, edges] = histcounts(tail_usage(choose_box_idx), bin);
        pdf(:, end+1) = N' / sum(N);
    end
    h = bar(pdf, 'grouped');
    set(h(1), 'facecolor', cc{3}); 
    set(h(2), 'facecolor', cc{4}); 
    set(gca, 'xticklabel', edges(1:end-1));
    xlabel('Box Usage Tail'); ylabel('PDF');
    h = legend('Tenant1 Size: 152', 'Tenant2 Size: 100');
    set(h, 'location', 'northeast', 'box', 'off');
    title('Mean Usage = 20')
    set(gca, 'fontsize', font_size);
    
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'pdf_tenant_box_tail_scatter_', mat2str(usage_tail)));
    
    %%%%%%%%%%%%%%%%%%%%% First Index to check %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the box usage tail *correlates* with tenant usage tail?
    fig = figure;
    set(fig, 'Position', [200 200 300 250]);
    font_size = 15;
%     subplot(1,2,1)
%     scatter(BOX_TENANT_USAGE_TAIL{end}(:,2), BOX_TENANT_USAGE_TAIL{end}(:,1));
%     xlabel(strcat('Tenant Usage Tail =', {' '}, mat2str(usage_tail), '%' ));
%     ylabel(strcat('Box Usage Tail =', {' '}, mat2str(usage_tail), '%' ));
%     set(gca, 'xlim', [0 100], 'xtick', [0 : 10 : 100]);
%     set(gca, 'ylim', [0 100], 'ytick', [0 : 10 : 100]);
%     title('All Boxes')
%     set(gca, 'fontsize', font_size);
    %view(30, -10)
    scatter(MEAN_BOX_TENANT_USAGE_TAIL{end}(:,2), MEAN_BOX_TENANT_USAGE_TAIL{end}(:,1), '*');
    hold on
    plot(prctile(all_usage{3}, usage_tail), mean_box_usage_tail(3), 'marker', '^', 'color', 'r');
    hold on
    plot(prctile(all_usage{4}, usage_tail), mean_box_usage_tail(4), 'marker', 'o', 'color', 'k');
    
    xlabel(strcat('Tenant Usage Tail =', {' '}, mat2str(usage_tail), '%' ));
    ylabel(strcat('Box Usage Tail =', {' '}, mat2str(usage_tail), '%' ));
    set(gca, 'xlim', [0 100], 'xtick', [0 : 20 : 100]);
    set(gca, 'ylim', [0 100], 'ytick', [0 : 20 : 100]);
    h = legend('All Tenants', 'Tenant1', 'Tenant2');
    set(h, 'location', 'northeast', 'box', 'off');
    % title('Mean Boxes')
    set(gca, 'fontsize', font_size);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'tenant_box_tail_scatter_', mat2str(usage_tail)));

    %%%%%%%%%%%%%%%%%%%%% Second Index to check %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cluster tenants based on the box tail
    if second_try
        total_metric_num = size(TENANT_BOX_USAGE_TAIL, 2);
        d = zeros(total_metric_num); d_half = zeros(total_metric_num);
        for row = 2 : total_metric_num
            idx1 = find(TENANT_BOX_USAGE_TAIL(:, row) ~= 0);
            for col = 1 : row -1
                idx2 = find(TENANT_BOX_USAGE_TAIL(:, col) ~= 0);
                idx = intersect(idx1, idx2);
                if numel(idx) == 0
                    % we don't have anything common
                    % I'm not sure if it is safe to do this
                    d(row, col) = 1;
                else
                    norm_row = norm(TENANT_BOX_USAGE_TAIL(idx, row));
                    norm_col = norm(TENANT_BOX_USAGE_TAIL(idx, col));
                    d(row, col) = 1 - dot(TENANT_BOX_USAGE_TAIL(idx, row), ...
                                      TENANT_BOX_USAGE_TAIL(idx, col)) / (norm_row * norm_col);
                end

                d(col, row) = d(row, col);
                d_half(row, col) = d(row, col);
                % in case of totally same case
                if d_half(row, col) == 0
                    d_half(row, col) = 0.0001;
                end
            end
        end 

        d_half_vec = reshape(d_half, 1, total_metric_num^2);
        d_vec = nonzeros(d_half_vec)';
        d_vec(find(d_vec == 0.0001)) = 0;
        z = linkage(d_vec, 'single'); % hierarchical clustering

        % Verify dissimilarity
        c = cophenet(z, d_vec);
        disp(strcat('The accuracy of clustering ([0,1]) is ', mat2str(c)));

         % Create clusters and pick up the optimal number of clusters
    %     s_diff_cluster = [];
    %     for num_cluster = 2 : ceil(sqrt(total_metric_num/2))
    %         cidx = cluster(z, 'maxclust', num_cluster);
    %         s_diff_cluster(num_cluster-1) = mean(silhouette([], cidx, d_vec));
    %     end   
    %     [max_s, idx] = max(abs(s_diff_cluster));
    %     disp(max_s);
    %     disp(num_of_cluster);

        num_of_cluster = 2;

        % do the optimal clustering
        cidx = cluster(z, 'maxclust', num_of_cluster);
        usage = 0 : 10 : 100;
        fig = figure; set(fig, 'Position', [200 200 300 250]);
        for cluster_id = 1 : num_of_cluster
            idx = find(cidx == cluster_id);
            for usage_id = 1 : numel(usage)
                usage_used_tail = TENANT_BOX_USAGE_TAIL(usage_id, idx);
                real_usage_tail = usage_used_tail(usage_used_tail ~= 0);
                if numel(real_usage_tail) == 0
                    continue;
                end
                mean_tail = mean(real_usage_tail); std_tail = std(real_usage_tail);
                plot(ones(3,1) * usage(usage_id) , [mean_tail-std_tail, mean_tail, mean_tail+std_tail], ...
                        'marker', marker_style{cluster_id}, 'color', cc{cluster_id});
    %             scatter(ones(1, numel(idx)) * usage(usage_id) , TENANT_BOX_USAGE_TAIL(usage_id, idx), ...
    %                     marker_style{cluster_id}, cc{cluster_id});    
                hold on
            end  
        end
        xlabel('Mean Box Usage');
        ylabel(strcat('Box Usage Tail =', {' '}, mat2str(usage_tail), '%' ));
        set(gca, 'xlim', [0 100], 'xtick', [0 : 10 : 100]);
        set(gca, 'ylim', [0 100], 'ytick', [0 : 10 : 100]);
        set(gca, 'fontsize', font_size);   
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(fig_path, 'tenant_box_tail_class_', mat2str(usage_tail)));
    end
    
end