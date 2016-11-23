% Tenant Classification
close all; clc

% clear
% load ../New_Data_7days/box_vm_time_series_summary_mem_only_with_zeros
% box_vm_time_series_mem = box_vm_time_series_summary;
% load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros
% % 
% moment_path = '../New_Data_box_moments_fit/moment_coeff/';
% load (strcat(moment_path, 'unique_tenant'));
% load (strcat(moment_path, 'tenant_box_series'));


day_thres = 5;
wind_length = 96;

moments = 5 : 5 : 95;
moments = [moments, 98];

server_num_thres = 10;
ticket_thres = 60;
usage_tail_idx = numel(moments) - 1;

% Check Parameters:
BOX_MEAN = []; BOX_MOMENTS = [];
MOMENT_MEAN = []; MOMENT_STD = [];
MOMENT_MEAN_COEFF = []; 

MOMENT_STD_COEFF_CLUSTER = [];
MOMENT_MEAN_COEFF_CLUSTER = [];

fig_path = strcat('../New_Data_box_moments_fit/','Tenant_Cluster_Fig/');
mkdir(fig_path);

font_size = 15;

thres = [1.1, 1.6];

%%%%%%%%%%%%%%%%%%%%%%%%% Zero Index to check %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% box spatial dependency
zero_try = false;
if zero_try
    mean_corr_coeff = []; representative_box_series = {};
    for tenant_id = 1 : numel(tenant_box_series)
        num_box = numel(tenant_box_series{tenant_id});
        if  num_box < 4
            break;
        end
            
        box_usage_tail = []; num_used_box = 0;
        box_usage_all = {};
        for box_id = 1 : num_box
            box_usage = tenant_box_series{tenant_id}{box_id}{1}(:,end-2); 
            if numel(box_usage) <= wind_length * 3
                continue;
            end
            
            box_usage_tail(end+1, :) = prctile(box_usage, moments)/nanmean(box_usage); 
            box_usage_all{end+1} = box_usage/nanmean(box_usage);  
            num_used_box = num_used_box + 1;
        end
      
        if num_used_box / num_box < 0.5
            continue;
        end
        
        all_coeff = []; all_coeff_matric = ones(num_used_box, num_used_box);
        for box_id = 1 : num_used_box
            for cand_id = box_id + 1 : num_used_box
                [h, pval]= kstest2(box_usage_all{box_id}, box_usage_all{cand_id});
                all_coeff(end+1, :) = [h, pval];
                all_coeff_matric(box_id, cand_id) = h;
                all_coeff_matric(cand_id, box_id) = h;
            end
        end
        
        [~, idx] = min(mean(all_coeff_matric));
        
        representative_box_series{end+1} = box_usage_all{idx};
        
        mean_corr_coeff(end+1, :) = mean(all_coeff);
    end
    
    tenant_num = numel(representative_box_series);
    corr_coeff_tenant = [];
    all_corr_coeff_tenant = ones(tenant_num, tenant_num);
    for tenant_id = 1 : tenant_num
        for cand_id = tenant_id+1 : tenant_num
            [h, pval] = kstest2(representative_box_series{tenant_id}, ...
                                representative_box_series{cand_id});   
            corr_coeff_tenant(end+1,:) = [h, pval];  
            
            all_corr_coeff_tenant(tenant_id, cand_id) = pval;
            all_corr_coeff_tenant(cand_id, tenant_id) = pval;
        end
    end
    mean_corr_coeff_tenant = mean(all_corr_coeff_tenant);
    
    fig = figure;
    set(fig, 'Position', [200 200 400 300]);
    [f1, x1] = ecdf(mean_corr_coeff(:, 2));
    [f2, x2] = ecdf(mean_corr_coeff_tenant);
    plot(x1, f1, 'r-', 'linewidth', 1.5);
    hold on
    plot(x2, f2, 'k--', 'linewidth', 1.5);
    hold on
    plot([0.05 0.05], [0 1], 'b-.', 'linewidth', 1.5)
    xlabel('p value of distribution fitting'); ylabel('CDF');
    set(gca, 'xlim',[0 0.1]); set(gca, 'xtick',[0 : 0.01 : 0.1]);
    set(gca, 'ylim',[0 1]); set(gca, 'ytick',[0 :0.1 : 1]);
    h = legend('Across Tenant', 'Across System','Significance Level');
    set(h, 'location', 'southeast');
    set(gca, 'fontsize', font_size);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'cdf_coeff_for_moments'));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%% First Index to check %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deterministic clustering tenants based on the usage tail
first_try = true;
if first_try
    % collect the mean and usage tail
    for box_id = 1 : numel(box_vm_time_series_summary)
        num_vm = numel(box_vm_time_series_summary{box_id});

        if num_vm == 0
            continue;
        end

        total_points = numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
        if total_points < wind_length * day_thres
            continue;
        end

        box_usage = box_vm_time_series_summary{1, box_id}{1,1}(:,4);
        
        % feature set: mean usage
        mean_box_usage = nanmean(box_usage);
        box_usage_moment = [];
        for moment_id = 1 : numel(moments)
            box_usage_moment(moment_id) = prctile(box_usage, moments(moment_id)); 
        end

        % gather all the information
        BOX_MEAN(end+1) = mean_box_usage;
        BOX_MOMENTS(end+1, :) = box_usage_moment;                
    end

    % calculate the mean for the usage tail per
    mean_d = 10;
    BOX_MEAN_DIS = round(BOX_MEAN / mean_d) * mean_d;
    UNIQUE_BOX_MEAN = unique(BOX_MEAN_DIS);
    for moment_id = numel(moments) - 3 : numel(moments)
        for mean_id = 1 : numel(UNIQUE_BOX_MEAN)
            if UNIQUE_BOX_MEAN(mean_id) >= ticket_thres
                break;
            end
            box_id = find(BOX_MEAN_DIS == UNIQUE_BOX_MEAN(mean_id));
            moment = BOX_MOMENTS(box_id, moment_id);
            % get the parameter for the normal distribution
            [MOMENT_MEAN(moment_id, mean_id), MOMENT_STD(moment_id, mean_id)] = normfit(moment);
        end
    end

    % calculate the coefficient for the usage
    for moment_id = numel(moments) - 3 : numel(moments)
        x = UNIQUE_BOX_MEAN(1:mean_id - 1)'; y = MOMENT_MEAN(moment_id, 1:mean_id-1)';
        fitresult = fit(x, y, 'poly1');
        MOMENT_MEAN_COEFF(moment_id, :) = polyfit(x, y, 1);              
    end 
    
    % cluster tenants based on the mean
    num_of_cluster = 3; 
    cluster_results = []; tenant_cluster_mean_tail = cell(1, num_of_cluster);
    for tenant_id = 1 : numel(tenant_box_series)
        num_box = numel(tenant_box_series{tenant_id});
        if  num_box < server_num_thres - 6
            break;
        end
       
        cluster_summary = []; mean_usage_tail_all = [];
        for box_id = 1 : num_box
            box_usage = tenant_box_series{tenant_id}{box_id}{1}(:,end-2); 
            if numel(box_usage) <= wind_length * (day_thres - 2)              
                continue;
            end
            
            mean_usage = nanmean(box_usage);
            box_usage_tail = prctile(box_usage, moments(usage_tail_idx));
            
            mean_usage_tail_all(end+1, :) = [mean_usage, box_usage_tail];
            
            usage_tail_thres = MOMENT_MEAN_COEFF(usage_tail_idx,1) * mean_usage ...
                               + MOMENT_MEAN_COEFF(usage_tail_idx,2);
            if usage_tail_thres >= box_usage_tail
                cluster_summary(end+1) = 1;
            else
                cluster_summary(end+1) = 2;
            end           
        end
      
        num_used_box = numel(cluster_summary);
        if num_used_box / num_box < 0.5
            continue;
        end
        
        if sum(cluster_summary) >= thres(2) * num_used_box
            cluster_results(end+1,1:2) = [tenant_id, 2];
            tenant_cluster_mean_tail{2}(end+1 : end+ num_used_box, :) = mean_usage_tail_all;
        elseif sum(cluster_summary) <= thres(1) * num_used_box
            cluster_results(end+1,1:2) = [tenant_id, 3];
            tenant_cluster_mean_tail{3}(end+1 : end+ num_used_box, :) = mean_usage_tail_all;
        else
            cluster_results(end+1,1:2) = [tenant_id, 1];
            tenant_cluster_mean_tail{1}(end+1 : end+ num_used_box, :) = mean_usage_tail_all;
        end 
        
        cluster_results(end, 3) = sum(cluster_summary) / num_used_box;
    end
    save(strcat(moment_path, 'cluster_results_',mat2str(moments(usage_tail_idx))), 'cluster_results');
    
    %%%%%%%%%%%%%%%%%% Plots: Check Cluster distribution %%%%%%%%%%%%%%%%%%
    fig = figure;
    set(fig, 'Position', [200 200 400 300]);
    [f, x] = ecdf(cluster_results(:,3));
    plot(x, f, 'r-', 'linewidth', 1.5);
    hold on
    plot([thres(1) thres(1)], [0 1], 'k--','linewidth', 1.5);
    hold on
    plot([thres(2) thres(2)], [0 1], 'k--','linewidth', 1.5);
    xlabel('Class Ranking'); ylabel('CDF');
    set(gca, 'xlim',[1 2]); set(gca, 'xtick',[1 : 0.1 : 2]);
    set(gca, 'ylim',[0 1]); set(gca, 'ytick',[0 :0.1 : 1]);
    set(gca, 'fontsize', font_size);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'cdf_class_ranking_', mat2str(moments(usage_tail_idx))));

    disp(strcat('Class 1 ', mat2str(sum(cluster_results(:,2) == 1) / numel(cluster_results(:,2)))));
    disp(strcat('Class 2 ', mat2str(sum(cluster_results(:,2) == 2) / numel(cluster_results(:,2)))));
    disp(strcat('Class 3 ', mat2str(sum(cluster_results(:,2) == 3) / numel(cluster_results(:,2)))));
    
    %%%%%%%%%% Plots: Check Cluster Mean/Tail Usage distribution %%%%%%%%%%
    fig = figure;
    set(fig, 'Position', [200 200 800 300]);
    subplot(1,2,1)
    [f1, x1] = ecdf(tenant_cluster_mean_tail{1}(:,1));
    [f2, x2] = ecdf(tenant_cluster_mean_tail{2}(:,1));
    [f3, x3] = ecdf(tenant_cluster_mean_tail{3}(:,1));
    plot(x1, f1, 'g-', 'linewidth', 1.5);
    hold on
    plot(x2, f2, 'r--','linewidth', 1.5);
    hold on
    plot(x3, f3, 'k-.','linewidth', 1.5);
    xlabel('Mean Box Usage'); ylabel('CDF');
    set(gca, 'xlim',[0 100]); set(gca, 'xtick',[0 : 10 : 100]);
    set(gca, 'ylim',[0 1]); set(gca, 'ytick',[0 :0.1 : 1]);
    h = legend('Class 1', 'Class 2', 'Class 3');
    set(h, 'location', 'southeast');
    set(gca, 'fontsize', font_size);
    
    subplot(1,2,2)
    [f1, x1] = ecdf(tenant_cluster_mean_tail{1}(:,2));
    [f2, x2] = ecdf(tenant_cluster_mean_tail{2}(:,2));
    [f3, x3] = ecdf(tenant_cluster_mean_tail{3}(:,2));
    plot(x1, f1, 'g-', 'linewidth', 1.5);
    hold on
    plot(x2, f2, 'r--','linewidth', 1.5);
    hold on
    plot(x3, f3, 'k-.','linewidth', 1.5);
    xlabel('Box Usage Tail'); ylabel('CDF');
    set(gca, 'xlim',[0 100]); set(gca, 'xtick',[0 : 10 : 100]);
    set(gca, 'ylim',[0 1]); set(gca, 'ytick',[0 :0.1 : 1]);
    h = legend('Class 1', 'Class 2', 'Class 3');
    set(h, 'location', 'southeast'); 
    set(gca, 'fontsize', font_size);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'cdf_mean_tail_', mat2str(moments(usage_tail_idx))));
    
    %%%%%%%%%%%%%%%% Plots: Check PDF of usage tail on mean %%%%%%%%%%%%%%%
    % check the distribution of monents to the mean: 
    test_dist = 'normal';
    dist_fit = - ones(numel(tenant_cluster_mean_tail), 11);
    mean_fit = zeros(numel(tenant_cluster_mean_tail), 11);
    std_fit = zeros(numel(tenant_cluster_mean_tail), 11);
    fitting_tail = []; fitting_std = [];
    for cluster_id = 1 : numel(tenant_cluster_mean_tail)
        box_usage_cluster = tenant_cluster_mean_tail{cluster_id};
        BOX_MEAN_DIS = round(box_usage_cluster(:,1)/ mean_d) * mean_d;
        UNIQUE_BOX_MEAN = 0 : mean_d : 100;
        for mean_id = 1 : numel(UNIQUE_BOX_MEAN)
            if UNIQUE_BOX_MEAN(mean_id) >= ticket_thres
                break;
            end
            
            box_id = find(BOX_MEAN_DIS == UNIQUE_BOX_MEAN(mean_id));
            if numel(box_id) == 0
                continue;
            end
            
            moment = box_usage_cluster(box_id, 2);
            
            [mean_fit(cluster_id, mean_id), std_fit(cluster_id, mean_id)] = normfit(moment);
            fit_data = normrnd(mean_fit(cluster_id, mean_id), std_fit(cluster_id, mean_id), [1000,1]);
            
            [h, p] = kstest2(moment, fit_data);
            dist_fit(cluster_id, mean_id) = p;
            
            if h 
                dist_fit(cluster_id, mean_id) = 1;
            else
                dist_fit(cluster_id, mean_id) = 0;
            end
            
            % to check the distribution
            fig = figure;
            set(fig, 'Position', [200 200 400 300]);
            histfit(moment, ceil(sqrt(numel(box_id))), test_dist);
            title(strcat('# of Samples = ', mat2str(numel(box_id)), ...
                         ', Mean = ', mat2str(UNIQUE_BOX_MEAN(mean_id))));
            xlabel(strcat(mat2str(moments(usage_tail_idx)), '%ile')); 
            ylabel(strcat('PDF'));
            set(gca,'fontsize', font_size);
            set(gcf, 'paperpositionmode', 'auto');
            print('-depsc2','-r300', strcat(fig_path, 'pdf_tail_of_mean_', ...
                        mat2str(moments(usage_tail_idx)),'_', mat2str(UNIQUE_BOX_MEAN(mean_id)), ...
                        '_cluster_', mat2str(cluster_id)));
            
            
        end
        
        % fit the mean and std for each cluster
        x = UNIQUE_BOX_MEAN(1 : mean_id - 1)';
        y1 = mean_fit(cluster_id, 1 : mean_id -1)';
        y2 = std_fit(cluster_id, 1 : mean_id -1)';
        fitresult = fit(x, y1, 'poly1');
        MOMENT_MEAN_COEFF_CLUSTER(cluster_id, :) = polyfit(x, y1, 1);
        fitresult2 = fit(x, y2, 'poly2');
        MOMENT_STD_COEFF_CLUSTER(cluster_id, :) = polyfit(x, y2, 2);
        
        fitting_tail(:, cluster_id) =  MOMENT_MEAN_COEFF_CLUSTER(cluster_id,1) * x + MOMENT_MEAN_COEFF_CLUSTER(cluster_id,2);
        fitting_std(:, cluster_id) = MOMENT_STD_COEFF_CLUSTER(cluster_id, 1) * (x .* x) ...
                                     + MOMENT_STD_COEFF_CLUSTER(cluster_id, 2) * x ...
                                     + MOMENT_STD_COEFF_CLUSTER(cluster_id, 3);
        
        fig = figure;
        set(fig, 'Position', [200 200 800 300]);
        subplot(1, 2, 1);
        plot(fitresult, x, y1);
        xlabel('Mean CPU Usage'); ylabel(strcat('mu of ', {' '}, mat2str(moments(usage_tail_idx)), 'th PCTile of CPU Usage'));
        set(gca, 'xlim',[0 ticket_thres]); set(gca, 'xtick',[0 : 10 : ticket_thres]);
        set(gca, 'ylim',[0 100]); set(gca, 'ytick',[0 : 10 : 100]);
        h = legend('data', 'fitted curve'); set(h,'location', 'southeast')
        set(gca, 'fontsize', font_size);
        subplot(1, 2, 2);
        plot(fitresult2, x, y2);
        xlabel('Mean CPU Usage'); ylabel(strcat('sigma of ', {' '}, mat2str(moments(usage_tail_idx)), 'th PCTile of CPU Usage'));
        set(gca, 'xlim',[0 ticket_thres]); set(gca, 'xtick',[0 : 10 : ticket_thres]);
        set(gca, 'ylim',[0 20]); set(gca, 'ytick',[0 : 5 : 20]);
        h = legend('data', 'fitted curve'); set(h,'location', 'southeast')
        set(gca, 'fontsize', font_size);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(fig_path, 'cluster_fit_moment_with_mean_pct_',mat2str(moments(usage_tail_idx))));
             
    end
    
    save(strcat(moment_path, 'MOMENT_STD_COEFF_CLUSTER_',mat2str(moments(usage_tail_idx))), 'MOMENT_STD_COEFF_CLUSTER');
    save(strcat(moment_path, 'MOMENT_MEAN_COEFF_CLUSTER_',mat2str(moments(usage_tail_idx))), 'MOMENT_MEAN_COEFF_CLUSTER');
    
    fig = figure;
    set(fig, 'Position', [200 200 400 300]);
    cc = {'g', 'r', 'k'}; line_style = {'-', '--', '-.'}; marker = {'*', 'o', '^'};
    for id = 1 : numel(tenant_cluster_mean_tail)
        errorbar(x, fitting_tail(:, id), fitting_std(:,id), 'color', cc{id}, 'linestyle', line_style{id}, ...
             'linewidth', 1.5, 'marker', marker{id});
        hold on
    end
    xlabel('Mean CPU Usage'); ylabel(strcat(mat2str(moments(usage_tail_idx)), 'th PCTile of CPU Usage'));
    set(gca, 'xlim',[0 ticket_thres - mean_d]); set(gca, 'xtick',[0 : 10 : ticket_thres - mean_d]);
    set(gca, 'ylim',[0 100]); set(gca, 'ytick',[0 : 10 : 100]);
    h = legend('Class 1', 'Class 2', 'Class 3'); set(h,'location', 'southeast')
    set(gca, 'fontsize', font_size);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'cluster_fit_moment_with_mean_pct_',mat2str(moments(usage_tail_idx)), '_all'));
    
    dist = reshape(dist_fit(1:end,1:mean_id - 1), [cluster_id * (mean_id - 1), 1]);
    disp(strcat('The mean p-value of distribution fitting is ', mat2str(mean(dist))));
    
   
    %%%%%%%%%%%%%%%%%%%% Plots: Check Tail distribution %%%%%%%%%%%%%%%%%%%
    fig = figure;
    set(fig, 'Position', [200 200 400 300]);
    scatter(tenant_cluster_mean_tail{1}(:,1),tenant_cluster_mean_tail{1}(:,2), '*', 'g');
    hold on
    scatter(tenant_cluster_mean_tail{2}(:,1),tenant_cluster_mean_tail{2}(:,2), 'o', 'r');
    hold on
    scatter(tenant_cluster_mean_tail{3}(:,1),tenant_cluster_mean_tail{3}(:,2), '^', 'k');
    hold on
    x = 0 : mean_d : ticket_thres; y = [];
    for i = 1 : numel(x)
        y(i) = MOMENT_MEAN_COEFF(usage_tail_idx, 1) * x(i) + MOMENT_MEAN_COEFF(usage_tail_idx, 2);
    end
    plot(x, y, 'b--', 'linewidth', 1.5);
    xlabel('Mean CPU Usage'); ylabel(strcat(mat2str(moments(usage_tail_idx)), 'th PCTile of CPU Usage'));
    set(gca, 'xlim',[0 ticket_thres - mean_d]); set(gca, 'xtick',[0 : 10 : ticket_thres - mean_d]);
    set(gca, 'ylim',[0 100]); set(gca, 'ytick',[0 : 10 : 100]);
    h = legend('Class 1', 'Class 2', 'Class 3'); set(h,'location', 'southeast')
    set(gca, 'fontsize', font_size);
    box on
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(fig_path, 'cluster_', mat2str(moments(usage_tail_idx))));
end