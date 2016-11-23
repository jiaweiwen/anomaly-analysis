close all; clear; clc

load ../New_Data/box_vm_time_series_summary_cpu_only
load ../New_Data/measure_capacity_box_vm_cpu
load ../New_Data/coeff_determine_original
coeff_determine_original = coeff_determine;

mkdir('../New_Data_CPU_WEIGHTED_VM_BOX_Fit_FIGURE_20150909');
path = '../New_Data_CPU_WEIGHTED_VM_BOX_Fit_FIGURE_20150909/';

mkdir('../New_Data_heatmap_CPU_VM_BOX');
path_heat = '../New_Data_heatmap_CPU_VM_BOX/';

mkdir('../New_Data_acf_CPU_BOX_itself');
path_acf = '../New_Data_acf_CPU_BOX_itself/';


grat = 900;

max_lag = 4;
lag_labels = {};
for lag = 0 : max_lag
    lag_labels{lag+1} = strcat(mat2str(lag*grat/3600), 'h');
end

coeff_determine = {}; no_box = 1;

corr_thres = 0.5; compare_coeff = [];

% We set this test_len to limit the influence from the overall trend change
for box_id = 1 : numel(box_vm_time_series_summary)
    size_box =  numel(box_vm_time_series_summary{box_id});
    if size_box == 0
        continue;
    end
    
    pm_id = box_vm_time_series_summary{box_id}{1,1}(1,2);
    
    time = box_vm_time_series_summary{box_id}{1,1}(:,1); 
    time_len = numel(time);
    sum_vm_cpu = zeros(time_len, 1);
    X = ones(time_len, 1);
    weights_original = [];
    for vm_id = 2 : size_box
        sum_vm_cpu = sum_vm_cpu + box_vm_time_series_summary{box_id}{1, vm_id}(1:time_len, 4);
        X = [X, box_vm_time_series_summary{box_id}{1, vm_id}(1:time_len, 4) / measure_capacity_box_vm{box_id}{1,vm_id}(1)];
        weights_original = [weights_original; measure_capacity_box_vm{box_id}{1, vm_id}(1)/measure_capacity_box_vm{box_id}{1, 1}(1)];
    end
    sum_box_itself_cpu = box_vm_time_series_summary{box_id}{1, 1}(1:time_len, 3) - sum_vm_cpu;
    
    % Get rid of zeros 
    [row_idx, col_idx] = find(box_vm_time_series_summary{box_id}{1,1}(1:time_len,3) > 0);
    
    % Assume that we use the *mean* to represent the total capacity in box
    sum_vm_cpu_pct = sum_vm_cpu(row_idx,:) / measure_capacity_box_vm{box_id}{1,1}(1);
    sum_box_itself_cpu_pct = sum_box_itself_cpu(row_idx,:) / measure_capacity_box_vm{box_id}{1,1}(1);
    sum_box_cpu_pct = box_vm_time_series_summary{box_id}{1,1}(row_idx,3)/ measure_capacity_box_vm{box_id}{1,1}(1);
    X = X(row_idx, :);
    time = time(row_idx, :);
    
    %% New Step: Check the correlation between VMs and BOX itself
    xcf = [];
    labels = {};
    for vm_id = 2 : size_box
        [xcf(vm_id-1,:) , lags, bounds] = crosscorr(X(:,vm_id), sum_box_itself_cpu_pct, max_lag); 
        labels{vm_id -1} = strcat('VM', mat2str(vm_id -1));
    end
    
    % we only need to look at the postive lags that means X moves to right
    
    % Check the heatmap
    fig = figure;
    set(fig,'Position',[200,200,400, 300]);
    heatmap(xcf(:, max_lag + 1 : end), lag_labels , labels, '%0.2f', 'TextColor', 'w', ...
            'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
            'TickAngle', 0, 'TickFontSize',10);
    caxis([-1 1]);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path_heat, 'PM_', mat2str(pm_id)));
    
    % Find the correlated traces
    used_lags_xcf = xcf(:, max_lag + 1 : end);
    great_than_corr_thres = abs(used_lags_xcf) >= 0;
 
    % Due to time lag, all the current observations should be delayed to
    % maximum lag
    sum_vm_cpu_pct = sum_vm_cpu_pct(max_lag+1 : end);
    sum_box_itself_cpu_pct = sum_box_itself_cpu_pct(max_lag+1 : end);
    sum_box_cpu_pct = sum_box_cpu_pct(max_lag+1 : end);
    time = time(max_lag + 1 : end);
    
    % Pick up related labels
    X_final = ones(numel(sum_vm_cpu_pct),1);
    for vm_id = 1 : size_box-1
        for lag_id = 1 : max_lag + 1
            % Only pick up the highlt-correlated vectors as the
            % features/variables to train the linear regression model
            if great_than_corr_thres(vm_id, lag_id) == 1
                X_final(:, end+1) = X(max_lag+2 - lag_id : end - lag_id + 1, vm_id + 1);
            end
        end
    end

    % Linear fitting the rest part of BOX CPU: 
    % 1) BOX itself 2) Relates with VMs 3) Residual
    
    % We have two different cases for the linear fitting
    % 1) X_final only has one column, no need to do fitting
    % 2) X_final has several VMs time series
    if numel(X_final(1,:)) == 1
        b = mean(sum_box_itself_cpu_pct);
        r = sum_box_itself_cpu_pct - b;
    else    
        [b, bint, r, rint, stats] = regress(sum_box_itself_cpu_pct, X_final);
    end
    
    % Box itself Part w/ residual 
    box_itself_cpu_pct_after_fit = b(1,1) + r;
    
    % Box relates with VMs
    box_with_vm_cpu_pct_after_fit = sum_box_itself_cpu_pct - box_itself_cpu_pct_after_fit;
    
    % Check the ACF of the box itself
    % Need to notice that, to check the autocorrelation, we have to fill in
    % the holes, due to it is really lag sensitive
    box_itself_cpu_pct_after_fit_without_gaps = [];
    time_begin = time(1,1); time_end = time(end); t = 1;
    time_lag = 0; total_lag = (time_end-time_begin)/grat;
    while time_lag <= total_lag
        if time_lag*grat+time_begin  == time(t)
            box_itself_cpu_pct_after_fit_without_gaps(time_lag+1,:) = ...
                                           box_itself_cpu_pct_after_fit(t);
        else
            while time(t) > time_lag*grat+time_begin
                box_itself_cpu_pct_after_fit_without_gaps(time_lag+1,:) = NaN;
                time_lag = time_lag + 1;
            end
            box_itself_cpu_pct_after_fit_without_gaps(time_lag+1,:) = ...
                                           box_itself_cpu_pct_after_fit(t);
        end
        t = t + 1;
        time_lag = time_lag + 1;
    end
    
    sum_box_cpu_pct_without_gaps = [];
    time_begin = time(1,1); time_end = time(end); t = 1;
    time_lag = 0; total_lag = (time_end-time_begin)/grat;
    while time_lag <= total_lag
        if time_lag*grat+time_begin  == time(t)
            sum_box_cpu_pct_without_gaps(time_lag+1,:) = ...
                                           sum_box_cpu_pct(t);
        else
            while time(t) > time_lag*grat+time_begin
                sum_box_cpu_pct_without_gaps(time_lag+1,:) = NaN;
                time_lag = time_lag + 1;
            end
            sum_box_cpu_pct_without_gaps(time_lag+1,:) = ...
                                           sum_box_cpu_pct(t);
        end
        t = t + 1;
        time_lag = time_lag + 1;
    end
    
    acf_lag = floor(total_lag * 0.6);
    acf = nanautocorr(box_itself_cpu_pct_after_fit_without_gaps, acf_lag);
    acf_original = nanautocorr(sum_box_cpu_pct_without_gaps, acf_lag);
    lags = 0 : grat : acf_lag*grat;
    
    fig = figure;
    set(fig,'Position',[200, 200, 1500, 300]);
    set(gca, 'fontsize', 15);
    plot(lags/(3600*24), acf, 'k');
    hold on
    plot(lags/(3600*24), acf_original, 'r:')
    h = legend('BOX Itself', 'BOX Overall');
    set(h, 'location', 'northeast', 'box','on', 'fontsize', 15);
    xlabel('Lag (day)', 'fontsize', 15); ylabel('ACF', 'fontsize', 15);
    set(gca, 'xtick', [0 : min(ceil(lags(end)/(3600*24)), 1) : ceil(lags(end)/(3600*24))],'fontsize', 15);
    set(gca, 'xlim', [0 ceil(lags(end)/(3600*24))]);
    set(gca, 'ytick',[-1 : 0.25 : 1],'fontsize', 15);
    set(gca, 'ylim', [-1 1]);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path_acf,strcat('PM_', mat2str(pm_id), '_box_acf')));
    
    time = time/(3600*24);
    % Plot the overtime figure
    fig = figure;
    set(fig,'Position',[200, 200, 1500, 700]);
    set(gca, 'fontsize', 15);
    subplot(2,1,1)
    plot(time, sum_box_cpu_pct, 'r-');
    hold on
    plot(time, sum_vm_cpu_pct, 'k:');
    hold on
    plot(time, sum_box_itself_cpu_pct, 'g--');
    max_util = ceil(max(sum_box_cpu_pct) * 100) / 100;
    h = legend('BOX', 'VM', 'BOX except VM');
    set(h, 'location', 'northeast', 'box','on', 'fontsize', 15);
    xlabel('Time (day)', 'fontsize', 15); ylabel('CPU USED PCT', 'fontsize', 15)
    set(gca, 'xtick', [floor(time(1)) : min(floor(time(end)) - ceil(time(1)), 3) : ceil(time(end))],'fontsize', 15);
    set(gca, 'xlim', [floor(time(1)) ceil(time(end))]);
    set(gca, 'ytick',[0 : max_util/5 : max_util],'fontsize', 15);
    set(gca, 'ylim', [0 max_util]);
    
    subplot(2,1,2)
    plot(time, box_with_vm_cpu_pct_after_fit, 'b-.');
    hold on
    plot(time, box_itself_cpu_pct_after_fit, 'm--');
    h = legend('BOX w/ VM', 'BOX Alone');
    set(h, 'location', 'northeast', 'box','on', 'fontsize', 15);
    xlabel('Time (day)', 'fontsize', 15); ylabel('CPU USED PCT', 'fontsize', 15)
    set(gca, 'xtick', [floor(time(1)) : min(floor(time(end)) - ceil(time(1)), 3) : ceil(time(end))],'fontsize', 15);
    set(gca, 'xlim', [floor(time(1)) ceil(time(end))]);
    max_util = ceil(max(max(box_with_vm_cpu_pct_after_fit),max(box_itself_cpu_pct_after_fit))*100)/100;
    min_util = min(0, ceil(min(min(box_with_vm_cpu_pct_after_fit),min(box_itself_cpu_pct_after_fit))*100)/100);
    util_len = max_util - min_util;
    set(gca, 'ytick',[min_util : util_len/5 : max_util],'fontsize', 15);
    set(gca, 'ylim', [min_util max_util]);
%     set(gca, 'yscale', 'log');
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_overtime_all')));
    
    % Plot error of fitting
    [fit_mu, fit_sigma] = normfit(r);
    fit_data = (fit_mu - 5*fit_sigma) : fit_sigma/100 : (fit_mu + 5*fit_sigma);    
    [f, x] = ecdf(r);
    bincounts = histc(r, fit_data);
    fit_pdf = normpdf(fit_data, fit_mu, fit_sigma);
    fig = figure;
    set(fig,'Position',[200, 200, 800, 300]);
    set(gca, 'fontsize', 15);
    
    subplot(1,2,1);
    plot(x, f, 'k-', 'linewidth', 1.5);
    h = legend('Residual of Prediction');
    set(h, 'location', 'southeast', 'box','on');
    xlabel('Error of CPU USED PCT'); ylabel('CDF')
    set(gca,'ytick',[ 0 : 0.2 : 1]); set(gca, 'ylim', [0 1]);
    set(gca,'xlim', [-max(abs(r)), max(abs(r))]);
    
    subplot(1,2,2)
    plot(fit_data, bincounts/sum(bincounts),'b:', 'linewidth', 1.5);
    hold on
    plot(fit_data, fit_pdf/sum(fit_pdf), 'r-.', 'linewidth', 1.5);
    h = legend('Residual of Prediction', 'Fitted Normal Distribution');
    set(h, 'location', 'northeast', 'box','on');
    xlabel('Error of CPU USED PCT'); ylabel('PDF')
    set(gca,'xlim', [-max(abs(r)), max(abs(r))]);
    title(strcat('Mean = ', mat2str(round(fit_mu, 5)), ', STD = ', mat2str(round(fit_sigma, 5))));
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_Error_CPU_PDF')));
    
    % Plot the PDF and CDF of the allocation
    test_box_itself_cpu = box_itself_cpu_pct_after_fit;
    [fit_mu, fit_sigma] = normfit(test_box_itself_cpu);
    fit_data = (fit_mu - 5*fit_sigma) : fit_sigma/100 : (fit_mu + 5*fit_sigma);    
    [f, x] = ecdf(test_box_itself_cpu);
    bincounts = histc(test_box_itself_cpu, fit_data);
    fit_pdf = normpdf(fit_data, fit_mu, fit_sigma);
    fig = figure;
    set(fig,'Position',[200, 200, 800, 300]);
    set(gca, 'fontsize', 15);
    
    subplot(1,2,1);
    plot(x, f, 'k-', 'linewidth', 1.5);
    h = legend('BOX Alone');
    set(h, 'location', 'southeast', 'box','on');
    xlabel('CPU USED PCT'); ylabel('CDF')
    set(gca,'ytick',[ 0 : 0.2 : 1]); set(gca, 'ylim', [0 1]);
    set(gca,'xlim', [min(test_box_itself_cpu), max(test_box_itself_cpu)]);
    
    subplot(1,2,2)
    plot(fit_data, bincounts/sum(bincounts),'b:', 'linewidth', 1.5);
    hold on
    plot(fit_data, fit_pdf/sum(fit_pdf), 'r-.', 'linewidth', 1.5);
    h = legend('Box except VM', 'Fitted Normal Distribution');
    set(h, 'location', 'northeast', 'box','on');
    xlabel('CPU USED PCT'); ylabel('PDF')
    set(gca,'xlim', [min(test_box_itself_cpu), max(test_box_itself_cpu)]);
    title(strcat('Mean = ', mat2str(round(fit_mu, 5)), ', STD = ', mat2str(round(fit_sigma, 5))));
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_Alone_CPU_PDF')));

    
    % Determine the goodness of fitting using R^2
    ssresid = sum(r.^2); 
    sstotal = (length(sum_box_cpu_pct)-1) * var(sum_box_cpu_pct);
    rsq = 1 - ssresid/sstotal;
    
    coeff_determine{no_box}{1} = [pm_id, rsq];
    
    compare_coeff(end+1,:) = [coeff_determine_original{no_box}{1}(2), rsq];
    no_box = no_box + 1;
    % Show the pm id
%     disp(strcat('PM is',mat2str(pm_id)));
%     disp(strcat('R^2 is', mat2str(rsq)));
    close all
end

path = '../New_Data_VM_BOX_Figure/';
[F_new, X_new] = ecdf(compare_coeff(:,2));
[F_original, X_original] = ecdf(compare_coeff(:,1));
fig = figure;
plot(X_new, F_new, 'k');
hold on
plot(X_original, F_original, 'r:');
h = legend('Fit with different lags', 'Fit with original');
set(h, 'location', 'northeast', 'box','on');
xlabel('Coefficient of determination'); ylabel('CDF')
title('R^2 Comparison');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_Fitting_Different_Method_R2_Comparison')));

