% this script aims at charactering BOX tickets per box

close all; clear; clc

load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros

path = '../New_Data_box_moments_fit/';
mkdir(path);

mkdir(strcat(path, 'Fig/'));

moments = 5 : 5 : 95;
moments = [moments, 98];

cc = hsv(numel(moments));

time_grat = 900;
day_point = 96;
ticket_cand = 60;

% related features
BOX_MOMENTS = [];
BOX_MEAN = [];

BOX_ID = [];

box_num = 1;
size_box_vm = size(box_vm_time_series_summary);

test_lag = 0;

fuzzy_bound = 0;

sig_vm_thres = 0.8;

USED_BOX = 0;

%%%%%%%%%%%%%%% Step 1:  
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
    
    box_usage = box_vm_time_series_summary{1, box_id}{1,1}(:,4);
    box_demands = box_vm_time_series_summary{1, box_id}{1,1}(:,3);
    
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

font_size = 18;
mean_d = 10;
BOX_MEAN_DIS = round(BOX_MEAN / mean_d) * mean_d;
UNIQUE_BOX_MEAN = unique(BOX_MEAN_DIS);
for moment_id = numel(moments) - 1 : numel(moments) - 1 %numel(moments) - 3 : numel(moments)
    fig = figure;
    set(fig, 'Position', [200 200 800 300]);
    subplot(1,2,1)
    scatter(BOX_MEAN, BOX_MOMENTS(:, moment_id));
    box on; 
    grid on;
    ylabel(strcat(mat2str(moments(moment_id)), '%ile')); 
    xlabel('Mean');
    set(gca, 'xlim',[0 ticket_cand]); set(gca, 'xtick',[0 : 10 :ticket_cand]);
    set(gca, 'ylim',[0 100]); set(gca, 'ytick',[0 : 10 : 100]);
    set(gca,'fontsize', font_size);
    subplot(1,2,2)
    for mean_id = 1 : numel(UNIQUE_BOX_MEAN) 
        box_id = find(BOX_MEAN_DIS == UNIQUE_BOX_MEAN(mean_id));
        scatter(BOX_MEAN_DIS(box_id), BOX_MOMENTS(box_id, moment_id));
        hold on
    end   
    box on; 
    grid on
    ylabel(strcat(mat2str(moments(moment_id)), '%ile')); 
    xlabel('Mean');
    set(gca, 'xlim',[0 ticket_cand]); set(gca, 'xtick',[0 : 10 : ticket_cand]);
    set(gca, 'ylim',[0 100]); set(gca, 'ytick',[0 : 10 : 100]);
    set(gca,'fontsize', font_size);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'Fig/moment_mean_',mat2str(moments(moment_id))));    
    
    fig = figure;
    set(fig, 'Position', [200 200 400 300]);
    for mean_id = 1 : numel(UNIQUE_BOX_MEAN) 
        box_id = find(BOX_MEAN_DIS == UNIQUE_BOX_MEAN(mean_id));
        [f_moment, x_moment] = ksdensity(BOX_MOMENTS(box_id, moment_id));
        scatter3(x_moment, ones(1, numel(x_moment)) * UNIQUE_BOX_MEAN(mean_id), f_moment);
        hold on
    end   
    view(120, 30);
    box on; 
    grid on;
    xlabel(strcat(mat2str(moments(moment_id)), '%ile')); 
    ylabel('Mean'); zlabel('PDF');
    set(gca, 'xlim',[0 100]); set(gca, 'xtick',[0 : 20 : 100]);
    set(gca, 'ylim',[0 ticket_cand]); set(gca, 'ytick',[0 : 10 : ticket_cand]);
    set(gca, 'zlim',[0 0.16]); set(gca, 'ztick',[0 : 0.04 : 0.16]);
    set(gca,'fontsize', font_size-3);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'Fig/moment_mean_3d_',mat2str(moments(moment_id))));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%% See if distribution follows %%%%%%%%%%%%%%%%%%%%
% discretization for the mean to get the distribution overview
check_moment_dist = true;
test_dist = 'normal';
mean_d = 5;
BOX_MEAN_DIS = round(BOX_MEAN / mean_d) * mean_d;
UNIQUE_BOX_MEAN = unique(BOX_MEAN_DIS);
if check_moment_dist
    font_size = 18;
    MOMENT_MEAN = []; MOMENT_STD = [];
    
    % check the distribution of monents to the mean: 
    for moment_id = numel(moments) - 1 : numel(moments) - 1 %numel(moments) - 3 : numel(moments)
        ALL_ERROR = [];
        for mean_id = 1 : numel(UNIQUE_BOX_MEAN)
            if UNIQUE_BOX_MEAN(mean_id) >= ticket_cand
                break;
            end
            box_id = find(BOX_MEAN_DIS == UNIQUE_BOX_MEAN(mean_id));
            moment = BOX_MOMENTS(box_id, moment_id);
            % to check the distribution
%             fig = figure;
%             set(fig, 'Position', [200 200 400 300]);
% 
%             h = histfit(moment, 20, test_dist);
%             title(strcat('# of Samples = ', mat2str(numel(box_id)), ...
%                          ', Mean = ', mat2str(UNIQUE_BOX_MEAN(mean_id))));
%             xlabel(strcat(mat2str(moments(moment_id)), '%ile')); 
%             ylabel(strcat('PDF'));
%             set(gca,'fontsize', font_size);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path, 'Fig/pdf_moment_fit_pct_',...
%                         mat2str(moments(moment_id)),'_', mat2str(UNIQUE_BOX_MEAN(mean_id))));
                    
            % get the parameter for the normal distribution
            [MOMENT_MEAN(moment_id, mean_id), MOMENT_STD(moment_id, mean_id)] = normfit(moment);
            
            % check the fitting errors
            upper_bound = MOMENT_MEAN(moment_id, mean_id) + 1 * MOMENT_STD(moment_id, mean_id);
            error = upper_bound - moment;
            ALL_ERROR = [ALL_ERROR; error ./ moment];
        end
        
        % plot the error
        fig = figure;
        set(fig, 'Position', [200 200 400 300]);  
        [f, x] = ecdf(ALL_ERROR * 100);
        plot(x, f, 'r', 'linewidth', 2);  
        grid on
        xlabel('PCT Error (%)'); ylabel('CDF');
        set(gca, 'xlim', [0 1], 'xtick', [0 : 0.2 : 1]);
        set(gca, 'xlim', [-50 500], 'xtick', [-50 : 50 : 500]);
        h = legend('1-sigma upper bound guess');
        set(h, 'location', 'southeast');
        xticklabel_rotate([], 45);
        title(strcat('Mean =', mat2str(round(nanmean(abs(ALL_ERROR * 100)),1)), ...
                     '%, STD = ', mat2str(round(nanstd(abs(ALL_ERROR * 100)),1)), '%'), ...
                     'fontweight', 'normal');
        set(gca, 'fontsize',15);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/fitting_error_upper'));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Fit the mean and std %%%%%%%%%%%%%%%%%%%%%%%%
    MOMENT_MEAN_COEFF = []; MOMENT_STD_COEFF = []; MOMENT_CONF_COEFF = [];
    for moment_id = numel(moments) - 1 : numel(moments) - 1 % numel(moments) - 3 : numel(moments)
        x = UNIQUE_BOX_MEAN(1:mean_id - 1)'; y = MOMENT_MEAN(moment_id, 1:mean_id-1)';
        y2 = MOMENT_STD(moment_id, 1:mean_id-1)';
        fitresult = fit(x, y, 'poly1');
        MOMENT_MEAN_COEFF(moment_id, :) = polyfit(x, y, 2);
        fitresult2 = fit(x, y2, 'poly2');
        MOMENT_STD_COEFF(moment_id, :) = polyfit(x, y2, 2);
        fig = figure;
        set(fig, 'Position', [200 200 800 300]);
        subplot(1, 2, 1);
        plot(fitresult, x, y);
        grid on
        xlabel('Mean CPU Usage'); ylabel(strcat('mu of ', mat2str(moments(moment_id)), 'th PCTile'));
        set(gca, 'xlim',[0 ticket_cand]); set(gca, 'xtick',[0 : 10 : ticket_cand]);
        set(gca, 'ylim',[0 100]); set(gca, 'ytick',[0 : 10 : 100]);
        h = legend('data', 'fitted curve'); set(h,'location', 'southeast')
        set(gca, 'fontsize', font_size);
        subplot(1, 2, 2);
        plot(fitresult2, x, y2);
        grid on
        xlabel('Mean CPU Usage'); ylabel(strcat('sigma of ', mat2str(moments(moment_id)), 'th PCTile'));
        set(gca, 'xlim',[0 ticket_cand]); set(gca, 'xtick',[0 : 10 : ticket_cand]);
        set(gca, 'ylim',[0 20]); set(gca, 'ytick',[0 : 5 : 20]);
        h = legend('data', 'fitted curve'); set(h,'location', 'southeast')
        set(gca, 'fontsize', font_size);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/fit_moment_with_mean_pct_',mat2str(moments(moment_id))));
    
        % fitting error just using the mean
        predict = polyval(MOMENT_MEAN_COEFF(moment_id, :), BOX_MEAN);
        all_error = (predict' - BOX_MOMENTS(:, moment_id)) ./ BOX_MOMENTS(:, moment_id);
        
        % plot the error
        fig = figure;
        set(fig, 'Position', [200 200 400 300]);  
        [f, x] = ecdf(all_error * 100);
        plot(x, f, 'r', 'linewidth', 2);   
        grid on
        xlabel('PCT Error (%)'); ylabel('CDF');
        set(gca, 'xlim', [0 1], 'xtick', [0 : 0.2 : 1]);
        set(gca, 'xlim', [-50 150], 'xtick', [-50 : 50 : 150]);
        % xticklabel_rotate([], 45);
        title(strcat('Mean =', mat2str(round(nanmean(abs(all_error * 100)),1)), ...
                     '%, STD = ', mat2str(round(nanstd(abs(all_error * 100)),1)), '%'), ...
                     'fontweight', 'normal');
        h = legend('Fitting Mean');
        set(h, 'location', 'southeast');
        set(gca, 'fontsize',15);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/fitting_error_mean'));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% Fit the 5(95%) confidence interval %%%%%%%%%%%%%%%%%%%
    MOMENT_CONF_COEFF_UPPER = MOMENT_MEAN_COEFF + 2 * MOMENT_STD_COEFF;
    MOMENT_CONF_COEFF_LOWER = MOMENT_MEAN_COEFF - 2 * MOMENT_STD_COEFF;
    CONF_UPPER = []; CONF_LOWER = [];
    % calculate the confidence interval 
    for moment_id = numel(moments) - 3 : numel(moments)
        for cand_mean_id = 1 : mean_id - 1
            used_mean = UNIQUE_BOX_MEAN(cand_mean_id);
            % calculate upper 
            p2 = MOMENT_CONF_COEFF_UPPER(moment_id, 1);
            p1 = MOMENT_CONF_COEFF_UPPER(moment_id, 2);
            p0 = MOMENT_CONF_COEFF_UPPER(moment_id, 3);       
            CONF_UPPER(moment_id, cand_mean_id) = max(min(100, used_mean^2 * p2 + used_mean * p1 + p0), 0);
            % calculate lower
            p2 = MOMENT_CONF_COEFF_LOWER(moment_id, 1);
            p1 = MOMENT_CONF_COEFF_LOWER(moment_id, 2);
            p0 = MOMENT_CONF_COEFF_LOWER(moment_id, 3); 
            CONF_LOWER(moment_id, cand_mean_id) = max(min(100, used_mean^2 * p2 + used_mean * p1 + p0), 0);
        end
    end
    
    moment_path = strcat(path, 'moment_coeff/'); mkdir(moment_path);
%     save(strcat(moment_path, 'MOMENT_CONF_COEFF_UPPER'), 'MOMENT_CONF_COEFF_UPPER');
%     save(strcat(moment_path, 'MOMENT_CONF_COEFF_LOWER'), 'MOMENT_CONF_COEFF_LOWER');
end



%%%%%%%%%%%%%%%%%%%%%%%%% check the moments with the mean %%%%%%%%%%%%%%%%%
check_moment_mean = false;
if check_moment_mean
    font_size = 18;
    for moment_id = 1 : numel(moments)
        if rem(moment_id, 2) == 1
            fig = figure;
            set(fig, 'Position', [200 200 800 300])
            subplot(1,2,1)
            scatter(BOX_MEAN, BOX_MOMENTS(:, moment_id));
            xlabel('Mean'); ylabel(strcat(mat2str(moments(moment_id)), 'th PCTile'))
            set(gca,'fontsize', font_size);        
        else
            subplot(1,2,2)
            scatter(BOX_MEAN, BOX_MOMENTS(:, moment_id));
            xlabel('Mean'); ylabel(strcat(mat2str(moments(moment_id)), 'th PCTile'))
            set(gca,'fontsize', font_size);
            set(gcf, 'paperpositionmode', 'auto');
            print('-depsc2','-r300', strcat(path, 'Fig/moment_fit_pct_',mat2str(moments(moment_id))));
        end    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Fitting the moments with mean %%%%%%%%%%%%%%%%%%%
fit_moment_mean = false;
ticket_cand = 60;
if fit_moment_mean
    font_size = 18;
    mean_d = 5;
    BOX_MEAN_DIS = ceil(BOX_MEAN / mean_d) * mean_d;
    check_box_less_than_thres_id = BOX_MEAN_DIS < ticket_cand;
    for moment_id = numel(moments) - 3 : numel(moments)      
        x = BOX_MEAN_DIS(check_box_less_than_thres_id)';
        y = BOX_MOMENTS(check_box_less_than_thres_id, moment_id);
        fitresult = fit(x, y, 'poly1');      
        predict = predint(fitresult, x, 0.95, 'functional', 'on');
        fig = figure;
        set(fig, 'Position', [200 200 400 300]);
        plot(fitresult, x, y);
%         hold on
%         plot(x, min(100, max(0, predict)), 'm--');
        xlabel('Mean'); ylabel(strcat(mat2str(moments(moment_id)), 'th PCTile'));
        set(gca, 'xlim',[0 ticket_cand]); set(gca, 'xtick',[0 : 10 : ticket_cand]);
        set(gca, 'ylim',[0 100]); set(gca, 'ytick',[0 : 10 : 100]);
        h = legend('data', 'fitted curve'); set(h,'location', 'southeast')
        set(gca, 'fontsize', font_size);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/fit_moment_with_mean_pct_',mat2str(moments(moment_id))));
       
        
    end
end

save(strcat(moment_path, 'BOX_MEAN'), 'BOX_MEAN');
save(strcat(moment_path, 'BOX_MOMENTS'), 'BOX_MOMENTS');