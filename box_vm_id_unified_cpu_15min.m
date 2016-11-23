close all
clear
clc

%% Step 0: Load all the raw data
load ../New_Data/box_differentID
% Grat = 15 min; 
% Time, BOX_ID, MEM_MB, MEM_PCT, CPU_MHZ, CPU_PCT
box = machine_cell;

load ../New_Data/vm_differentID
% Grat = 15 min; 
% Time, VM_ID, BOX_ID, MEM_MB, MEM_PCT, CPU_MHZ, CPU_PCT
% VM contains the BOX information, so we also need to divide based on BOX
vm = machine_cell;

%% Define the threshold to determine giving tickets or not
% Use the 80 percentile of PCT to serve as the threshold
all_box_cpu_pct = [];
for box_id = 1 : numel(box)
    all_box_cpu_pct = [all_box_cpu_pct; box{box_id}(:, 6)];
end

all_vm_cpu_pct = [];
for vm_id = 1 : numel(vm)
    all_vm_cpu_pct = [all_vm_cpu_pct; vm{vm_id}(:, 7)];
end

all_box_cpu_pct_great_than_zero_index = find(all_box_cpu_pct >= 0);
all_vm_cpu_pct_great_than_zero_index = find(all_vm_cpu_pct >= 0);

all_cpu_pct = [all_box_cpu_pct(all_box_cpu_pct_great_than_zero_index); ...
               all_vm_cpu_pct(all_vm_cpu_pct_great_than_zero_index)];

pct = 80;
cpu_ticket_thres = prctile(all_cpu_pct, pct)

mkdir('../New_Data_VM_BOX_Figure');
path = '../New_Data_VM_BOX_Figure/';

% First, checking the distribution of VM and BOX
[F_box, X_box] = ecdf(all_box_cpu_pct(all_box_cpu_pct_great_than_zero_index));
[F_vm, X_vm] = ecdf(all_vm_cpu_pct(all_vm_cpu_pct_great_than_zero_index));
[F_all, X_all] = ecdf(all_cpu_pct);
fig = figure;
set(fig, 'Position', [200 200 400 300]);
plot(X_box, F_box, 'r-');
hold on
plot(X_vm, F_vm, 'k--');
set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0:0.1:1]);
set(gca, 'xlim', [0 100]); set(gca, 'xtick',[0:10:100]);
xlabel('CPU USED PCT (%)'); ylabel('CDF');
h = legend('BOX', 'VM');
set(h, 'box', 'off', 'location', 'southeast');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_BOX_CPU_PCT_CDF'));

%% Step 1: Fill time gaps for BOX and VM
grat_small = 900; 
box_fill_gaps = {};
vm_fill_gaps = {};
for box_id = 1 : numel(box)
   [box_fill_gaps{box_id}] = fill_gaps(box{box_id}, grat_small);
end

for vm_id = 1 : numel(vm)
   [vm_fill_gaps{vm_id}] = fill_gaps(vm{vm_id}, grat_small);
end

%% Step 2: Group the different VMs into the same PM (with small granularity)
% The structure is: PM, VM0, VM1, ...
% The structure for PM is: PM_id, Time, CPU_PCT, MEM_PCT
% The structure for VM is: VM_id, Time, PM_id, CPU_PCT, MEM_PCT, DISK_PCT,
% NET_TX, NET_RX

box_vm_time_series = {}; box_idx = [];
% First put the PM into each cell and record the corresponding BOX ID
for box_id = 1 : numel(box_fill_gaps)
    box_vm_time_series{box_id}{1} = box_fill_gaps{box_id}(:,1:2);
    box_vm_time_series{box_id}{1} = [box_vm_time_series{box_id}{1}, box_fill_gaps{box_id}(:,5:6)];
    box_idx(box_id) = box_fill_gaps{box_id}(1,2);
end

% Second for each VM, seperate by VM ID and PM ID
% Each VM could be allocated to several PM
vm_fill_gaps_diff_box = {};
vm_idx = []; vm_box_idx = {};
summary_of_vm_diff_box = {};
for vm_id = 1 : numel(vm_fill_gaps)
    no_box = 1; time_id = 1; box_id = vm_fill_gaps{vm_id}(1,3);
    start_time_id = 1; end_time_id = numel(vm_fill_gaps{vm_id}(:,1));
    while time_id <= end_time_id
        if (box_id ~= vm_fill_gaps{vm_id}(time_id,3) ...
            && vm_fill_gaps{vm_id}(time_id,3) ~= -1) || time_id == end_time_id
            vm_fill_gaps_diff_box{vm_id}{no_box} = vm_fill_gaps{vm_id}...
                                                   (start_time_id:time_id-1,:);
            start_time_id = time_id; 
            box_id = vm_fill_gaps{vm_id}(time_id,3); no_box = no_box + 1;
        end
        time_id = time_id + 1;
    end
    vm_fill_gaps_diff_box{vm_id}{no_box-1} = ...
    [vm_fill_gaps_diff_box{vm_id}{no_box-1}; vm_fill_gaps{vm_id}(time_id-1, :)];
    
    for box_id = 1 : no_box - 1
        % Pick up the needed performance metrics
        summary_of_vm_diff_box{vm_id}{box_id} = ...
            vm_fill_gaps_diff_box{vm_id}{box_id}(:,1:3);
        summary_of_vm_diff_box{vm_id}{box_id}(:,4:5) = ...
            vm_fill_gaps_diff_box{vm_id}{box_id}(:,6:7);

%         % Pre-fill the empty fields for disk and net
%         summary_of_vm_diff_box{vm_id}{box_id}(:, 6:8) = 0;
        
        % Record the corresponding BOX ID
        vm_box_idx{vm_id}(box_id) = vm_fill_gaps_diff_box{vm_id}{box_id}(1,3);      
    end
    % Record the corresponding VM ID
    vm_idx(vm_id) = vm_fill_gaps{vm_id}(1,2);
end

%% Step 3: Combine VM into PM
for vm_id = 1 : numel(vm_idx)
    for box_id = 1 : numel(vm_box_idx{vm_id})
        box_id_final = find(box_idx == vm_box_idx{vm_id}(box_id));
        % Ideally, the box_id_final shoud only has one element
        if numel(box_id_final) ~= 1
            disp('we could not locate the box id');
            disp(numel(box_id_final));
            continue;
        else
            box_vm_time_series{box_id_final}{end+1} = summary_of_vm_diff_box{vm_id}{box_id};
        end
    end    
end

save('../New_Data/box_vm_time_series_cpu_only','box_vm_time_series')

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
        
        % Change the time series into tickets
        % we only need to change CPU USED PCT into ticket
%         all_metric_no = numel(box_vm_time_series_summary{box_id}{cand_id}(1,:));
%         box_vm_time_series_ticket_summary{box_id}{cand_id}(:, 1:all_metric_no-2) = ...
%             box_vm_time_series_summary{box_id}{cand_id}(:,1:all_metric_no-2);
%         box_vm_time_series_ticket_summary{box_id}{cand_id}(:, all_metric_no-1) = ...
%             box_vm_time_series_summary{box_id}{cand_id}(:,all_metric_no-1) >= cpu_ticket_thres;  
%         box_vm_time_series_ticket_summary{box_id}{cand_id}(:, all_metric_no) = ...
%             box_vm_time_series_summary{box_id}{cand_id}(:,all_metric_no);  
    end
    
end

save('../New_Data/box_vm_time_series_summary_cpu_only','box_vm_time_series_summary');
% save('../box_vm_time_series_ticket_summary_cpu_only', 'box_vm_time_series_ticket_summary');

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

        capacity_box_vm{box_id}{1, vm_id} = box_vm_cand(...
                 box_vm_time_series_summary{box_id}{1, vm_id}(:, end) > 0);

        % Summary BOX capacity
        measure_capacity_box_vm{box_id}{1, vm_id} = [mean(capacity_box_vm{box_id}{1, vm_id}), ...
                                           median(capacity_box_vm{box_id}{1, vm_id}), ...
                                           prctile(capacity_box_vm{box_id}{1, vm_id}, 10), ...
                                           prctile(capacity_box_vm{box_id}{1, vm_id}, 90), ...
                                           std(capacity_box_vm{box_id}{1, vm_id})/...
                                           mean(capacity_box_vm{box_id}{1, vm_id})];  
   end
end

save('../New_Data/measure_capacity_box_vm_cpu','measure_capacity_box_vm');
disp('hahha');

% Plot the figure of capacity for BOX and VMs
% BOX first
fig = figure;
set(fig,'Position',[200, 200, 400, 300]);
set(gca, 'fontsize', 18);
cc = hsv(4); marker_cand = {'*', 'o', '^', '>'};
num_box = 1;
for box_id = 1 : numel(box_vm_time_series_summary)
    size_box =  numel(box_vm_time_series_summary{box_id});
    if size_box == 0
        continue;
    end  
   
    for metric_no = 1 : 4
        plot(num_box, measure_capacity_box_vm{box_id}{1, 1}(metric_no), 'color', ...
             cc(metric_no, :), 'marker', marker_cand{metric_no}, 'markersize', 6);
        hold on
    end
    num_box = num_box + 1;
    h = legend('Mean', 'Median', '10 PRCTILE', '90 PRCTILE');
    set(h, 'location', 'northeast', 'box','on');
    %disp(strcat('BOX COV = ', mat2str(measure_capacity_box_vm{box_id}{1, 1}(5))));
end
set(gca, 'yscale', 'log');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'BOX_CPU_VALIDATION'));

fig = figure;
set(fig,'Position',[200, 200, 400, 300]);
set(gca, 'fontsize', 18);
cc = hsv(4); marker_cand = {'*', 'o', '^', '>'};
num_box = 1;
for box_id = 1 : numel(box_vm_time_series_summary)
    size_box =  numel(box_vm_time_series_summary{box_id});
    if size_box == 0
        continue;
    end  
   
    for vm_id = 2 : size_box
        for metric_no = 1 : 4
            plot(num_box, measure_capacity_box_vm{box_id}{1, vm_id}(metric_no), 'color', ...
                 cc(metric_no, :), 'marker', marker_cand{metric_no}, 'markersize', 6);
            hold on
        end
        num_box = num_box + 1;
        %disp(strcat('VM COV = ', mat2str(measure_capacity_box_vm{box_id}{1, vm_id}(5))));
    end
    h = legend('Mean', 'Median', '10 PRCTILE', '90 PRCTILE');
    set(h, 'location', 'northeast', 'box','on');
end
set(gca, 'yscale', 'log');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_CPU_VALIDATION'));

%% Validation of VM and BOX Relationship
mkdir('../New_Data_CPU_WEIGHTED_VM_BOX_Fit_FIGURE');
path = '../New_Data_CPU_WEIGHTED_VM_BOX_Fit_FIGURE/';
% test_len = 7 * 24 * 3600 / grat_small;
% for box_id = 1 : numel(box_vm_time_series_summary)
%     size_box =  numel(box_vm_time_series_summary{box_id});
%     if size_box == 0
%         continue;
%     end 
%     pm_id = box_vm_time_series_summary{box_id}{1,1}(1,1);
%     time = box_vm_time_series_summary{box_id}{1,1}(:,2)/(3600*24); 
%     time_len = numel(time);
%     sum_vm_cpu = zeros(time_len, 1);
%     for vm_id = 2 : size_box
%         sum_vm_cpu = sum_vm_cpu + box_vm_time_series_summary{box_id}{1, vm_id}(:, 5);
%     end
%     sum_box_itself_cpu = box_vm_time_series_summary{box_id}{1, 1}(:, 4) - sum_vm_cpu;
%     
%     % Assume that we use the *mean* to represent the total capacity in box
%     sum_vm_cpu_pct = sum_vm_cpu / measure_capacity_box_vm{box_id}{1,1}(1);
%     sum_box_itself_cpu_pct = sum_box_itself_cpu / measure_capacity_box_vm{box_id}{1,1}(1);
%     sum_box_cpu_pct = box_vm_time_series_summary{box_id}{1,1}(:,4) / measure_capacity_box_vm{box_id}{1,1}(1);
%     
%     % Plot the overtime for summed CPU and comparison
%     fig = figure;
%     set(fig,'Position',[200, 200, 1500, 300]);
%     set(gca, 'fontsize', 15);
%     plot(time(1:end), sum_box_cpu_pct(1:end), 'r-');
%     hold on
%     plot(time(1:end), sum_vm_cpu_pct(1:end), 'k:');
%     hold on
%     plot(time(1:end), sum_box_itself_cpu_pct(1:end), 'b-.');
%     max_util = ceil(max(sum_box_cpu_pct) * 100) / 100;
%     h = legend('BOX', 'VM Alone', 'BOX Alone');
%     set(h, 'location', 'northeast', 'box','on', 'fontsize', 15);
%     xlabel('Time (hour)', 'fontsize', 15); ylabel('CPU USED PCT', 'fontsize', 15)
%     set(gca, 'xtick', [time(1) : 5 : time(end)],'fontsize', 15);
%     set(gca, 'xlim', [time(1) time(end)]);
%     set(gca, 'ytick',[0 : max_util/5 : max_util],'fontsize', 15);
%     set(gca, 'ylim', [0 max_util]);
% %     set(gca, 'yscale', 'log');
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_overtime_all')));
%     
%     % Plot the PDF and CDF of the allocation
%     [non_zero_row, non_zero_col] = find(sum_box_cpu_pct);
%     test_box_itself_cpu = sum_box_itself_cpu_pct(non_zero_row);
%     [f, x] = ecdf(test_box_itself_cpu);
%     [n_1, x_1] = pdfplot(test_box_itself_cpu, 20, 0, max(test_box_itself_cpu));
%     fig = figure;
%     set(fig,'Position',[200, 200, 800, 300]);
%     set(gca, 'fontsize', 15);
%     
%     subplot(1,2,1);
%     plot(x, f, 'k-', 'linewidth', 1.5);
%     h = legend('BOX Alone');
%     set(h, 'location', 'southeast', 'box','on');
%     xlabel('CPU USED PCT'); ylabel('CDF')
%     set(gca,'ytick',[ 0 : 0.2 : 1]); set(gca, 'ylim', [0 1]);
%     set(gca,'xlim', [0, max(test_box_itself_cpu)]);
%     
%     subplot(1,2,2)
%     plot(x_1, n_1/sum(n_1),'k-', 'linewidth', 1.5);
%     h = legend('BOX Alone');
%     set(h, 'location', 'northeast', 'box','on');
%     xlabel('CPU USED PCT'); ylabel('PDF')
%     set(gca,'xlim', [0, max(test_box_itself_cpu)]);
%     set(gcf, 'paperpositionmode', 'auto');
%     print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_Alone_CPU_PDF')));
% end

%% Compute the R^2 (coefficient of determination)
test_len = 35 * 24 * 3600 / grat_small;
coeff_determine = {}; no_box = 1;
% We set this test_len to limit the influence from the overall trend change
for box_id = 1 : numel(box_vm_time_series_summary)
    size_box =  numel(box_vm_time_series_summary{box_id});
    if size_box == 0
        continue;
    end
    
    pm_id = box_vm_time_series_summary{box_id}{1,1}(1,2);
    
    time = box_vm_time_series_summary{box_id}{1,1}(1:test_len,1)/(3600*24); 
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

    % Linear fitting the rest part of BOX CPU: 
    % 1) BOX itself 2) Relates with VMs 3) Residual
    [b, bint, r, rint, stats] = regress(sum_box_itself_cpu_pct, X);
    
    % Box itself Part w/ residual 
    box_itself_cpu_pct_after_fit = b(1,1) + r;
    
    % Box relates with VMs
    box_with_vm_cpu_pct_after_fit = sum_box_itself_cpu_pct - box_itself_cpu_pct_after_fit;
    
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
    set(gca, 'xtick', [time(1) : ceil(time(end)/7) : time(end)],'fontsize', 15);
    set(gca, 'xlim', [time(1) time(end)]);
    set(gca, 'ytick',[0 : max_util/5 : max_util],'fontsize', 15);
    set(gca, 'ylim', [0 max_util]);
    
    subplot(2,1,2)
    plot(time, box_with_vm_cpu_pct_after_fit, 'b-.');
    hold on
    plot(time, box_itself_cpu_pct_after_fit, 'm--');
    h = legend('BOX w/ VM', 'BOX Alone');
    set(h, 'location', 'northeast', 'box','on', 'fontsize', 15);
    xlabel('Time (day)', 'fontsize', 15); ylabel('CPU USED PCT', 'fontsize', 15)
    set(gca, 'xtick', [time(1) : ceil(time(end)/7) : time(end)],'fontsize', 15);
    set(gca, 'xlim', [time(1) time(end)]);
    max_util = ceil(max(max(box_with_vm_cpu_pct_after_fit),max(box_itself_cpu_pct_after_fit))*100)/100;
    set(gca, 'ytick',[0 : max_util/5 : max_util],'fontsize', 15);
    set(gca, 'ylim', [0 max_util]);
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
    set(gca,'xlim', [0, max(test_box_itself_cpu)]);
    title(strcat('Mean = ', mat2str(round(fit_mu, 5)), ', STD = ', mat2str(round(fit_sigma, 5))));
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path,strcat('PM_', mat2str(pm_id), '_Alone_CPU_PDF')));

    
    % Determine the goodness of fitting using R^2
    ssresid = sum(r.^2); 
    sstotal = (length(sum_box_cpu_pct)-1) * var(sum_box_cpu_pct);
    rsq = 1 - ssresid/sstotal;
    
    coeff_determine{no_box}{1} = [pm_id, rsq]; no_box = no_box + 1;
    % Show the pm id
%     disp(strcat('PM is',mat2str(pm_id)));
%     disp(strcat('R^2 is', mat2str(rsq)));
    close all
end

path = '../New_Data_VM_BOX_Figure/';
fig = figure;
set(fig,'Position',[200, 200, 800, 400]);
set(gca, 'fontsize', 18); tick_name = {};
for box_id = 1 : no_box -1
    plot(box_id, coeff_determine{box_id}{1,1}(1,2), 'ro', 'markersize', 6)
    hold on
    tick_name{box_id} = mat2str(coeff_determine{box_id}{1,1}(1,1));
end
set(gca,'xtick', 1 : no_box -1, 'xticklabel', tick_name);
set(gca, 'ylim', [0.4 1]); set(gca, 'ytick', [0.4 : 0.1 : 1]);
xlabel('PM ID'); ylabel('Coefficient of Determination');
h = legend('Weighted Sum + Linear Regression');
set(h,'location', 'southeast','box','on');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path,'VM_CPU_WEIGHTED_FITTING_COMPARISON'));


% %% Step 6: Validate the sum weighted with the linear regression fitting
% % Attention: before do the linear fitting, we need to remove all the zeros
% compare_coeff = {};
% test_len = 35 * 24 * 3600 / grat_small; no_box = 1;
% % We set this test_len to limit the influence from the overall trend change
% for box_id = 1 : numel(box_vm_time_series_summary)
%     size_box =  numel(box_vm_time_series_summary{box_id});
%     if size_box == 0
%         continue;
%     end
%     
%     pm_id = box_vm_time_series_summary{box_id}{1,1}(1,1);
%     
%     % Find all the non-zeros index
%     [row_idx, col_idx] = find(box_vm_time_series_summary{box_id}{1,1}(1:test_len,3) > 0);
%     
%     Y = box_vm_time_series_summary{box_id}{1, 1}(row_idx,3);
%     X = ones([numel(row_idx), 1]); weights_original = [];
%     for vm_id = 2 : size_box
%         X = [X, box_vm_time_series_summary{box_id}{1, vm_id}(row_idx, 4)];
%         weights_original = [weights_original; measure_capacity_box_vm{box_id}{1, vm_id}(1)/measure_capacity_box_vm{box_id}{1, 1}(1)];
%     end
%     
%     % Do the multi-linear fitting
%     [b, bint, r, rint, stats] = regress(Y, X);
%     
%     % Compare with the original weight
%     compare_coeff = {compare_coeff{:}, [b(2:end), weights_original]};  
%     
%     coeff_determine{no_box}{2} = [pm_id, stats(1)]; no_box = no_box + 1;
%     
%     % Show the pm id
% %     disp(strcat('PM is',mat2str(pm_id)));
% %     disp(strcat('R^2 is', mat2str(stats(1))));
%     
% end
% 
% % Plot the coefficient of determination
% path = '../BOX_VM_CPU_CAPACITY_VALIDATION_FIGURE/';
% fig = figure;
% set(fig,'Position',[200, 200, 800, 400]);
% set(gca, 'fontsize', 18); tick_name = {};
% for box_id = 1 : no_box -1
%     plot(box_id, coeff_determine{box_id}{1,1}(1,2), 'ro', 'markersize', 6)
%     hold on
%     plot(box_id, coeff_determine{box_id}{1,2}(1,2), 'k*', 'markersize', 6)
%     hold on
%     tick_name{box_id} = mat2str(coeff_determine{box_id}{1,1}(1,1));
% end
% set(gca,'xtick', 1 : no_box -1, 'xticklabel', tick_name);
% set(gca, 'ylim', [0.4 1]); set(gca, 'ytick', [0.4 : 0.1 : 1]);
% xlabel('PM ID'); ylabel('Coefficient of Determination');
% h = legend('Weighted Sum',' Linear Regression');
% set(h,'location', 'southeast','box','on');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path,'VM_CPU_FITTING_COMPARISON'));
% 
% 
% 
% 
% 
% 
% 
