close all
clear
clc

%% Step 0: Load all the raw data
load ../New_Data/box_differentID
% Grat = 15 min; 
% Time, BOX_ID, MEM_MB, MEM_PCT, MEM_MHZ, MEM_PCT
box = machine_cell;

load ../New_Data/vm_differentID
% Grat = 15 min; 
% Time, VM_ID, BOX_ID, MEM_MB, MEM_PCT, MEM_MHZ, MEM_PCT
% VM contains the BOX information, so we also need to divide based on BOX
vm = machine_cell;

%% Define the threshold to determine giving tickets or not
% Use the 80 percentile of PCT to serve as the threshold
% all_box_mem_pct = [];
% for box_id = 1 : numel(box)
%     all_box_mem_pct = [all_box_mem_pct; box{box_id}(:, 6)];
% end
% 
% all_vm_mem_pct = [];
% for vm_id = 1 : numel(vm)
%     all_vm_mem_pct = [all_vm_mem_pct; vm{vm_id}(:, 7)];
% end

% all_box_mem_pct_great_than_zero_index = find(all_box_mem_pct >= 0);
% all_vm_mem_pct_great_than_zero_index = find(all_vm_mem_pct >= 0);
% 
% all_mem_pct = [all_box_mem_pct(all_box_mem_pct_great_than_zero_index); ...
%                all_vm_mem_pct(all_vm_mem_pct_great_than_zero_index)];
% 
% pct = 80;
% mem_ticket_thres = prctile(all_mem_pct, pct)
% 
% mkdir('../New_Data_VM_BOX_Figure');
% path = '../New_Data_VM_BOX_Figure/';
% 
% % First, checking the distribution of VM and BOX
% [F_box, X_box] = ecdf(all_box_mem_pct(all_box_mem_pct_great_than_zero_index));
% [F_vm, X_vm] = ecdf(all_vm_mem_pct(all_vm_mem_pct_great_than_zero_index));
% [F_all, X_all] = ecdf(all_mem_pct);
% fig = figure;
% set(fig, 'Position', [200 200 400 300]);
% plot(X_box, F_box, 'r-');
% hold on
% plot(X_vm, F_vm, 'k--');
% hold on
% plot(X_all, F_all, 'b-.');
% set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0:0.1:1]);
% set(gca, 'xlim', [0 100]); set(gca, 'xtick',[0:10:100]);
% xlabel('MEM USED PCT (%)'); ylabel('CDF');
% h = legend('BOX', 'VM', 'BOX+VM');
% set(h, 'box', 'off', 'location', 'southeast');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'VM_BOX_MEM_PCT_CDF'));

% Clear data
all_mem_pct = []; 
all_box_mem_pct = []; all_vm_mem_pct = [];
all_box_mem_pct_great_than_zero_index = [];
all_vm_mem_pct_great_than_zero_index = [];
F_box = []; X_box = []; F_vm = []; X_vm = []; F_all = []; X_all =[];

%% Step 2: Group the different VMs into the same PM (with small granularity)
% The structure is: PM, VM0, VM1, ...
% The structure for PM is: PM_id, Time, MEM_PCT, MEM_PCT
% The structure for VM is: VM_id, Time, PM_id, MEM_PCT, MEM_PCT, DISK_PCT,
% NET_TX, NET_RX

box_vm_time_series = {}; box_idx = [];
% First put the PM into each cell and record the corresponding BOX ID
for box_id = 1 : numel(box)
    box_vm_time_series{box_id}{1} = box{box_id}(:,1:2);
    box_vm_time_series{box_id}{1} = [box_vm_time_series{box_id}{1}, box{box_id}(:,3:4)];
    
    % remove the holes   
    not_hole_idx = box_vm_time_series{box_id}{1}(:,3) >= 0;
    box_vm_time_series{box_id}{1} = box_vm_time_series{box_id}{1}(not_hole_idx, :);
    
    box_idx(box_id) = box{box_id}(1,2);
end

% Clear data
box = {};

% Second for each VM, seperate by VM ID and PM ID
% Each VM could be allocated to several PM
vm_idx = []; vm_box_idx = {};
summary_of_vm_diff_box = {};
for vm_id = 1 : numel(vm)
    % remove holes
    not_hole_idx = vm{vm_id}(:, 4) >= 0;
    vm{vm_id} = vm{vm_id}(not_hole_idx, :);
    
    temp_vm = sortrows(vm{vm_id}, [3 1]);
    time_len = numel(temp_vm(:,1));
    box_id = temp_vm(1,3); box_no = 1;
    time_idx = 1; start_idx = 1;
    while time_idx <= time_len
        if box_id ~= temp_vm(time_idx, 3)
            summary_of_vm_diff_box{vm_id}{box_no} = temp_vm(start_idx:time_idx-1,1:3);
            summary_of_vm_diff_box{vm_id}{box_no}(:,4:5) = temp_vm(start_idx:time_idx-1,4:5);
            % Record the corresponding BOX ID
            vm_box_idx{vm_id}(box_no) = box_id;  
            box_no = box_no + 1; start_idx = time_idx;
            box_id = temp_vm(time_idx, 3);
        end
        time_idx = time_idx + 1;
    end
    % Record the end of left metrics
    vm_box_idx{vm_id}(box_no) = box_id; 
    summary_of_vm_diff_box{vm_id}{box_no} = temp_vm(start_idx:time_idx-1,1:3);
    summary_of_vm_diff_box{vm_id}{box_no}(:,4:5) = temp_vm(start_idx:time_idx-1,4:5);
    
    % Record the corresponding VM ID
    vm_idx(vm_id) = temp_vm(1,2);
end

% Clear data
vm = {};

%% Step 3: Combine VM into PM
for vm_id = 1 : numel(vm_idx)
    for box_id = 1 : numel(vm_box_idx{vm_id})
        box_id_final = find(box_idx == vm_box_idx{vm_id}(box_id));
        % Ideally, the box_id_final shoud only has one element
        if numel(box_id_final) ~= 1
            %disp('we could not locate the box id');
            %disp(numel(box_id_final));
            continue;
        else
            box_vm_time_series{box_id_final}{end+1} = summary_of_vm_diff_box{vm_id}{box_id};
        end
    end    
end

save('../New_Data/box_vm_time_series_mem_only','box_vm_time_series')

% Clear data
summary_of_vm_diff_box = {};

%% Step 4: Combine VM and PM into one *small* time series
% This part makes the VMs and PM share the same time stamps to calculate
% the correlation among the big time series. The first try is to get some
% examples that all the VM have same resindual time.
box_vm_time_series_summary = {}; box_vm_time_series_ticket_summary = {};
for box_id = 1 : numel(box_vm_time_series)
    if numel(box_vm_time_series{box_id}) == 1
        disp('This box has no VM residing on it');
        continue;
    end
    
    % Based on time stamp do intersection
    common_time = box_vm_time_series{box_id}{1}(:,1);
    for cand_id = 2 : numel(box_vm_time_series{box_id})
        common_time = intersect(common_time, box_vm_time_series{box_id}{cand_id}(:,1));
    end
    
    if numel(common_time) == 0
        disp('This box has no VM concurrently residing on it');
        continue;
    end
    
    common_time_idx = {}; 
    for cand_id = 1 : numel(box_vm_time_series{box_id})
        [common_time_series, common_idx, common_time_idx{cand_id}] = ...
        intersect(common_time, box_vm_time_series{box_id}{cand_id}(:,1));
        box_vm_time_series_summary{box_id}{cand_id} = ...
            box_vm_time_series{box_id}{cand_id}(common_time_idx{cand_id},:); 
    end
    
end

save('../New_Data/box_vm_time_series_summary_mem_only','box_vm_time_series_summary');
% save('../box_vm_time_series_ticket_summary_mem_only', 'box_vm_time_series_ticket_summary');

% Clear data
box_vm_time_series = {};

%% Step 5: Validation of Capacity in both BOX an VM are constant
capacity_box_vm = {}; measure_capacity_box_vm = {};
for box_id = 1 : numel(box_vm_time_series_summary);
    size_box =  numel(box_vm_time_series_summary{box_id});
    if size_box == 0
        continue;
    end                              
        
   for vm_id = 1 : size_box
       % First validate the BOX capacity
        box_vm_cand = box_vm_time_series_summary{box_id}{1, vm_id}(:, end - 1) ./ ...
                      box_vm_time_series_summary{box_id}{1, vm_id}(:, end) * 100;

        capacity_box_vm{box_id}{1, vm_id} = box_vm_cand;

        % Summary BOX capacity
        measure_capacity_box_vm{box_id}{1, vm_id} = [mean(capacity_box_vm{box_id}{1, vm_id}), ...
                                           median(capacity_box_vm{box_id}{1, vm_id}), ...
                                           prctile(capacity_box_vm{box_id}{1, vm_id}, 10), ...
                                           prctile(capacity_box_vm{box_id}{1, vm_id}, 90), ...
                                           std(capacity_box_vm{box_id}{1, vm_id})/...
                                           mean(capacity_box_vm{box_id}{1, vm_id})];  
   end
end

save('../New_Data/measure_capacity_box_vm_mem','measure_capacity_box_vm');
disp('hahha');

% Plot the figure of capacity for BOX and VMs
% BOX first
% fig = figure;
% set(fig,'Position',[200, 200, 400, 300]);
% set(gca, 'fontsize', 18);
% cc = hsv(4); marker_cand = {'*', 'o', '^', '>'};
% num_box = 1;
% for box_id = 1 : numel(box_vm_time_series_summary)
%     size_box =  numel(box_vm_time_series_summary{box_id});
%     if size_box == 0
%         continue;
%     end  
%    
%     for metric_no = 1 : 4
%         plot(num_box, measure_capacity_box_vm{box_id}{1, 1}(metric_no), 'color', ...
%              cc(metric_no, :), 'marker', marker_cand{metric_no}, 'markersize', 6);
%         hold on
%     end
%     num_box = num_box + 1;
%     h = legend('Mean', 'Median', '10 PRCTILE', '90 PRCTILE');
%     set(h, 'location', 'northeast', 'box','on');
%     %disp(strcat('BOX COV = ', mat2str(measure_capacity_box_vm{box_id}{1, 1}(5))));
% end
% set(gca, 'yscale', 'log');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'BOX_MEM_VALIDATION'));
% 
% fig = figure;
% set(fig,'Position',[200, 200, 400, 300]);
% set(gca, 'fontsize', 18);
% cc = hsv(4); marker_cand = {'*', 'o', '^', '>'};
% num_box = 1;
% for box_id = 1 : numel(box_vm_time_series_summary)
%     size_box =  numel(box_vm_time_series_summary{box_id});
%     if size_box == 0
%         continue;
%     end  
%    
%     for vm_id = 2 : size_box
%         for metric_no = 1 : 4
%             plot(num_box, measure_capacity_box_vm{box_id}{1, vm_id}(metric_no), 'color', ...
%                  cc(metric_no, :), 'marker', marker_cand{metric_no}, 'markersize', 6);
%             hold on
%         end
%         num_box = num_box + 1;
%         %disp(strcat('VM COV = ', mat2str(measure_capacity_box_vm{box_id}{1, vm_id}(5))));
%     end
%     h = legend('Mean', 'Median', '10 PRCTILE', '90 PRCTILE');
%     set(h, 'location', 'northeast', 'box','on');
% end
% set(gca, 'yscale', 'log');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'VM_MEM_VALIDATION'));

%% Validation of VM and BOX Relationship and Compute the R^2 (coefficient 
%  of determination)
mkdir('../New_Data_MEM_WEIGHTED_VM_BOX_Fit_FIGURE');
path = '../New_Data_MEM_WEIGHTED_VM_BOX_Fit_FIGURE/';

coeff_determine = {}; no_box = 1; zero_box = 0;
all_average_ape = []; all_coeff_determin = [];
box_all_ratio = [];
box_activity_ratio = []; vm_activity_ratio = [];
% We set this test_len to limit the influence from the overall trend change
for box_id = 1 : numel(box_vm_time_series_summary)
    size_box =  numel(box_vm_time_series_summary{box_id});
    if size_box == 0
        zero_box = zero_box + 1;
        continue;
    end
    
    pm_id = box_vm_time_series_summary{box_id}{1,1}(1,2);
    
    time = box_vm_time_series_summary{box_id}{1,1}(:,1)/(3600*24); 
    time_len = numel(time);
    % The following method is not safe, this is because that the memory
    % allocated to each VM is variate overtime or just errors?!
    sum_vm_mem = zeros(time_len, 1);
    X = ones(time_len, 1);
%     weights_original = [];
    for vm_id = 2 : size_box
        sum_vm_mem = sum_vm_mem + box_vm_time_series_summary{box_id}{1, vm_id}(1:time_len, 4);
        X = [X, box_vm_time_series_summary{box_id}{1, vm_id}(1:time_len, 5)/100];
%         weights_original = [weights_original; measure_capacity_box_vm{box_id}{1, vm_id}(1)/measure_capacity_box_vm{box_id}{1, 1}(1)];
    end
%     sum_box_itself_mem = box_vm_time_series_summary{box_id}{1, 1}(1:time_len, 3) - sum_vm_mem;
%     
    % Get rid of zeros 
    [row_idx, col_idx] = find(box_vm_time_series_summary{box_id}{1,1}(1:time_len,4) > 0);
    
    % Assume that we use the *mean* to represent the total capacity in box
%     sum_vm_mem_pct = sum_vm_mem(row_idx,:) / measure_capacity_box_vm{box_id}{1,1}(1);
%     sum_box_itself_mem_pct = sum_box_itself_mem(row_idx,:) / measure_capacity_box_vm{box_id}{1,1}(1);
    sum_box_mem_pct = box_vm_time_series_summary{box_id}{1,1}(row_idx,4)/100;
    X = X(row_idx, :);
    time = time(row_idx, :);

    % Linear fitting BOX MEM: 
    % 1) BOX itself 2) Relates with VMs (sum weighted + influence) 3) Residual
    [b, bint, r, rint, stats] = regress(sum_box_mem_pct, X);
    
    % Box itself Part w/ residual 
    box_itself_mem_pct_after_fit = b(1,1);
    
    % Box relates with VMs
    % box_with_vm_mem_pct_after_fit = sum_box_itself_mem_pct - box_itself_mem_pct_after_fit;
    
    % Box activity ratio
    box_activity_ratio_temp = box_itself_mem_pct_after_fit ./ sum_box_mem_pct;
    vm_activity_ratio_temp = 1 - box_activity_ratio_temp;
    
    if mean(box_activity_ratio_temp) < 0
        continue;
    end
    
    box_all_ratio_temp = sum_box_mem_pct(sum_box_mem_pct >= 0);
    box_all_ratio_temp = box_all_ratio_temp(box_all_ratio_temp <= 1);
    
    box_activity_ratio_temp = box_activity_ratio_temp(box_activity_ratio_temp>=0);
    vm_activity_ratio_temp = vm_activity_ratio_temp(box_activity_ratio_temp>=0);
    vm_activity_ratio_temp = vm_activity_ratio_temp(vm_activity_ratio_temp <= 1);
    
    % Plot the overtime figure
%     fig = figure;
%     set(fig,'Position',[200, 200, 1500, 700]);
%     set(gca, 'fontsize', 15);
%     subplot(2,1,1)
%     plot(time, sum_box_mem_pct, 'r-');
%     hold on
%     plot(time, sum_vm_mem_pct, 'k:');
%     hold on
%     plot(time, sum_box_itself_mem_pct, 'g--');
%     max_util = ceil(max(sum_box_mem_pct) * 100) / 100;
%     h = legend('BOX', 'VM', 'BOX except VM');
%     set(h, 'location', 'northeast', 'box','on', 'fontsize', 15);
%     xlabel('Time (day)', 'fontsize', 15); ylabel('MEM USED PCT', 'fontsize', 15)
%     set(gca, 'xtick', [floor(time(1)) : min(floor(time(end)) - ceil(time(1)), 3) : ceil(time(end))],'fontsize', 15);
%     set(gca, 'xlim', [floor(time(1)) ceil(time(end))]);
%     set(gca, 'ytick',[0 : max_util/5 : max_util],'fontsize', 15);
%     set(gca, 'ylim', [0 max_util]);
%     
%     subplot(2,1,2)
%     plot(time, box_with_vm_mem_pct_after_fit, 'b-.');
%     hold on
%     plot(time, box_itself_mem_pct_after_fit, 'm--');
%     h = legend('BOX w/ VM', 'BOX Alone');
%     set(h, 'location', 'northeast', 'box','on', 'fontsize', 15);
%     xlabel('Time (day)', 'fontsize', 15); ylabel('MEM USED PCT', 'fontsize', 15)
%     set(gca, 'xtick', [floor(time(1)) : min(floor(time(end)) - ceil(time(1)), 3) : ceil(time(end))],'fontsize', 15);
%     set(gca, 'xlim', [floor(time(1)) ceil(time(end))]);
%     max_util = ceil(max(max(box_with_vm_mem_pct_after_fit),max(box_itself_mem_pct_after_fit))*100)/100;
%     min_util = min(0, ceil(min(min(box_with_vm_mem_pct_after_fit),min(box_itself_mem_pct_after_fit))*100)/100);
%     util_len = max_util - min_util;
%     set(gca, 'ytick',[min_util : util_len/5 : max_util],'fontsize', 15);
%     set(gca, 'ylim', [min_util max_util]);
% %     set(gca, 'yscale', 'log');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_overtime_all')));
%     
%     % Plot error of fitting
%     [fit_mu, fit_sigma] = normfit(r);
%     fit_data = (fit_mu - 5*fit_sigma) : fit_sigma/100 : (fit_mu + 5*fit_sigma);    
%     [f, x] = ecdf(r);
%     bincounts = histc(r, fit_data);
%     fit_pdf = normpdf(fit_data, fit_mu, fit_sigma);
%     fig = figure;
%     set(fig,'Position',[200, 200, 800, 300]);
%     set(gca, 'fontsize', 15);
%     
%     subplot(1,2,1);
%     plot(x, f, 'k-', 'linewidth', 1.5);
%     h = legend('Residual of Prediction');
%     set(h, 'location', 'southeast', 'box','on');
%     xlabel('Error of MEM USED PCT'); ylabel('CDF')
%     set(gca,'ytick',[ 0 : 0.2 : 1]); set(gca, 'ylim', [0 1]);
%     set(gca,'xlim', [-max(abs(r)), max(abs(r))]);
%     
%     subplot(1,2,2)
%     plot(fit_data, bincounts/sum(bincounts),'b:', 'linewidth', 1.5);
%     hold on
%     plot(fit_data, fit_pdf/sum(fit_pdf), 'r-.', 'linewidth', 1.5);
%     h = legend('Residual of Prediction', 'Fitted Normal Distribution');
%     set(h, 'location', 'northeast', 'box','on');
%     xlabel('Error of MEM USED PCT'); ylabel('PDF')
%     set(gca,'xlim', [-max(abs(r)), max(abs(r))]);
%     title(strcat('Mean = ', mat2str(round(fit_mu, 5)), ', STD = ', mat2str(round(fit_sigma, 5))));
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_Error_MEM_PDF')));
%     
%      % Plot the PDF and CDF of the allocation
%     test_box_itself_mem = box_itself_mem_pct_after_fit;
%     [fit_mu, fit_sigma] = normfit(test_box_itself_mem);
%     fit_data = (fit_mu - 5*fit_sigma) : fit_sigma/100 : (fit_mu + 5*fit_sigma);    
%     [f, x] = ecdf(test_box_itself_mem);
%     bincounts = histc(test_box_itself_mem, fit_data);
%     fit_pdf = normpdf(fit_data, fit_mu, fit_sigma);
%     fig = figure;
%     set(fig,'Position',[200, 200, 800, 300]);
%     set(gca, 'fontsize', 15);
%     
%     subplot(1,2,1);
%     plot(x, f, 'k-', 'linewidth', 1.5);
%     h = legend('BOX Alone');
%     set(h, 'location', 'southeast', 'box','on');
%     xlabel('MEM USED PCT'); ylabel('CDF')
%     set(gca,'ytick',[ 0 : 0.2 : 1]); set(gca, 'ylim', [0 1]);
%     set(gca,'xlim', [min(test_box_itself_mem), max(test_box_itself_mem)]);
%     
%     subplot(1,2,2)
%     plot(fit_data, bincounts/sum(bincounts),'b:', 'linewidth', 1.5);
%     hold on
%     plot(fit_data, fit_pdf/sum(fit_pdf), 'r-.', 'linewidth', 1.5);
%     h = legend('Box except VM', 'Fitted Normal Distribution');
%     set(h, 'location', 'northeast', 'box','on');
%     xlabel('MEM USED PCT'); ylabel('PDF')
%     set(gca,'xlim', [min(test_box_itself_mem), max(test_box_itself_mem)]);
%     title(strcat('Mean = ', mat2str(round(fit_mu, 5)), ', STD = ', mat2str(round(fit_sigma, 5))));
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_Alone_MEM_PDF')));
 
    ape = abs(r) ./ sum_box_mem_pct;
    all_average_ape(end+1,1:2) = [mean(ape), prctile(ape, 90)];
    
    box_activity_ratio(end+1,1:2) = [mean(box_activity_ratio_temp),prctile(box_activity_ratio_temp, 90)];
      
    vm_activity_ratio(end+1,1:2) = [mean(vm_activity_ratio_temp),prctile(vm_activity_ratio_temp, 90)];
    
    box_all_ratio(end+1, 1:2) = [mean(box_all_ratio_temp),prctile(box_all_ratio_temp, 90)];

    % Determine the goodness of fitting using R^2
    ssresid = sum(r.^2); 
    sstotal = (length(sum_box_mem_pct)-1) * var(sum_box_mem_pct);
    rsq = 1 - ssresid/sstotal;
    
    if isnan(rsq)
        rsq = 0;
    end
    
    all_coeff_determin(end+1) = rsq;
    
    coeff_determine{no_box}{1} = [pm_id, rsq]; no_box = no_box + 1;
    % Show the pm id
%     disp(strcat('PM is',mat2str(pm_id)));
%     disp(strcat('R^2 is', mat2str(rsqf)));
    close all
end
disp(no_box);

path = '../New_Data_VM_BOX_Figure/';
fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca,'FontSize',18)
[f_r, x_r] = ecdf(all_coeff_determin);
plot(x_r, 1-f_r, 'k', 'linewidth', 2)
set(gca,'xlim', [0 1]); set(gca, 'xtick', [0 : 0.1 : 1], 'fontsize', 18);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1]);
ylabel('CCDF', 'fontsize' ,18); 
xlabel('Coefficient of Determination for Linear Fitting (R^2)','fontsize' ,18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_MEM_WEIGHTED_FITTING_R^2'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_ave, x_ave] = ecdf(all_average_ape(:,1) * 100);
[f_prt, x_prt] = ecdf(all_average_ape(:,2) * 100);
plot(x_ave, f_ave, 'k', 'linewidth', 2)
hold on
plot(x_prt, f_prt, 'r--', 'linewidth', 2)
h = legend('Mean APE', '90%ile APE');
set(h, 'box','on','location','southeast','fontsize',18);
set(gca,'xlim', [0 20]); set(gca, 'xtick', [0 : 5 : 20]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('Absolute Percentage Error (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_MEM_WEIGHTED_FITTING_ERROR'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_ave, x_ave] = ecdf(box_activity_ratio(:,1) * 100);
[f_prt, x_prt] = ecdf(box_activity_ratio(:,2) * 100);
plot(x_ave, f_ave, 'k', 'linewidth', 2)
hold on
plot(x_prt, f_prt, 'r--', 'linewidth', 2)
h = legend('Mean PCT of BOX Activity', '90%ile PCT of BOX Activity');
set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('PCT of BOX Activity on Total Usage (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_MEM_WEIGHTED_BOX_RATIO'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_ave, x_ave] = ecdf(vm_activity_ratio(:,1) * 100);
[f_prt, x_prt] = ecdf(vm_activity_ratio(:,2) * 100);
plot(x_ave, f_ave, 'k', 'linewidth', 2)
hold on
plot(x_prt, f_prt, 'r--', 'linewidth', 2)
h = legend('Mean PCT of VM Activity', '90%ile PCT of VM Activity');
set(h, 'box','on','location','southeast','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('PCT of VM Activity on Total Usage (%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_MEM_WEIGHTED_VM_RATIO'));

fig = figure;
set(fig,'Position',[200, 200, 600, 400]);
set(gca, 'fontsize', 18);
[f_ave, x_ave] = ecdf(box_all_ratio(:,1) * 100);
[f_prt, x_prt] = ecdf(box_all_ratio(:,2) * 100);
plot(x_ave, f_ave, 'k', 'linewidth', 2)
hold on
plot(x_prt, f_prt, 'r--', 'linewidth', 2)
h = legend('Mean BOX USED PCT ', '90%ile BOX USED PCT');
set(h, 'box','on','location','northwest','fontsize',18);
set(gca,'xlim', [0 100]); set(gca, 'xtick', [0 : 20 : 100]);
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0 : 0.1 : 1], 'fontsize', 18);
% set(gca, 'xscale', 'log');
ylabel('CDF', 'fontsize', 18); 
xlabel('BOX MEM USED PCT(%)', 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_MEM_WEIGHTED_BOX_ALL_RATIO'));
