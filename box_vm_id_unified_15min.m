close all
clear
clc

%% Step 0: Load all the raw data
load ../box_differentID
% Grat = 15 min; CPU_USED_PCT 	CPU_USED_MHZ	CPU_COUNT	MEM_USED_PCT	MEM_USED_MB	
box = machine_cell;

load ../disk_differentID
% Grat = 15 min; DISK_USED_PCT	DISK_USED_BYTES
disk = machine_cell;

load ../vm_differentID
% Grat = 15 min; CPU_USED_PCT	CPU_USED_MHZ	CPU_COUNT	MEM_USED_PCT	MEM_USED_MB	
% VM contains the BOX information, so we also need to divide based on BOX
vm = machine_cell;

%% Define the threshold to determine giving tickets or not
% Use the 80 percentile of PCT to serve as the threshold
all_box_cpu_pct = []; all_box_mem_pct = [];
for box_id = 1 : numel(box)
    all_box_cpu_pct = [all_box_cpu_pct; box{box_id}(:,3)];
    all_box_mem_pct = [all_box_mem_pct; box{box_id}(:,6)];
end

all_vm_cpu_pct = []; all_vm_mem_pct = [];
for vm_id = 1 : numel(vm)
    all_vm_cpu_pct = [all_vm_cpu_pct; vm{vm_id}(:,4)];
    all_vm_mem_pct = [all_vm_mem_pct; vm{vm_id}(:,7)];
end

all_cpu_pct = [all_box_cpu_pct; all_vm_cpu_pct];
all_mem_pct = [all_box_mem_pct; all_vm_mem_pct];

pct = 80;
cpu_ticket_thres = prctile(all_cpu_pct, pct);
mem_ticket_thres = prctile(all_mem_pct, pct);
% disp(prctile(all_box_cpu_pct, pct));
% disp(prctile(all_vm_cpu_pct, pct));
% disp(prctile(all_cpu_pct, pct));
% disp(prctile(all_box_mem_pct, pct));
% disp(prctile(all_vm_mem_pct, pct));
% disp(prctile(all_mem_pct, pct));

figure;
[f_1, x_1] = ecdf(all_box_cpu_pct);
[f_2, x_2] = ecdf(all_vm_cpu_pct);
[f_3, x_3] = ecdf(all_cpu_pct);
plot(x_1, f_1, 'r--', 'linewidth', 2);
hold on
plot(x_2, f_2, 'b-', 'linewidth', 2);
hold on 
plot(x_3, f_3, 'k-.', 'linewidth', 2);

set(gca, 'fontsize', 15);
title('CPU CDF');
h = legend('BOX', 'VM', 'BOX+VM');
set(h, 'location', 'southeast', 'box','off');
set(gcf, 'paperpositionmode', 'auto');
set(gca, 'xlim', [0 100]);
print('-depsc2','-r300', 'CPU_PCT_CDF');

figure;
[f_1, x_1] = ecdf(all_box_mem_pct);
[f_2, x_2] = ecdf(all_vm_mem_pct);
[f_3, x_3] = ecdf(all_mem_pct);
plot(x_1, f_1, 'r--', 'linewidth', 2);
hold on
plot(x_2, f_2, 'b-', 'linewidth', 2);
hold on
plot(x_3, f_3, 'k-.', 'linewidth', 2);
title('MEM CDF')
set(gca,'xscale','log');
set(gca, 'fontsize', 15);
set(gca, 'xlim', [0 100]);
title('MEM CDF');
h = legend('BOX', 'VM', 'BOX + VM');
set(h, 'location', 'southeast', 'box','off');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', 'MEM_PCT_CDF');

%% Step 1: Change the granularity for all the metrics
grat_small = 900; 
box_fill_gaps = {};
disk_fill_gaps = {};
vm_fill_gaps = {};
for box_id = 1 : numel(box)
   [box_fill_gaps{box_id}] = fill_gaps(box{box_id}, grat_small);
end

for vm_id = 1 : numel(disk)
   [disk_fill_gaps{vm_id}] = fill_gaps(disk{vm_id}, grat_small);
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
    box_vm_time_series{box_id}{1} = box_fill_gaps{box_id}(:,1:3);
    box_vm_time_series{box_id}{1}(:, 4) = box_fill_gaps{box_id}(:, 6);
    box_idx(box_id) = box_fill_gaps{box_id}(1,1);
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
            && vm_fill_gaps{vm_id}(time_id,3) ~= 0) || time_id == end_time_id
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
            vm_fill_gaps_diff_box{vm_id}{box_id}(:,1:4);
        summary_of_vm_diff_box{vm_id}{box_id}(:,5) = ...
            vm_fill_gaps_diff_box{vm_id}{box_id}(:,7);
        
%         % Pre-fill the empty fields for disk and net
%         summary_of_vm_diff_box{vm_id}{box_id}(:, 6:8) = 0;
        
        % Record the corresponding BOX ID
        vm_box_idx{vm_id}(box_id) = vm_fill_gaps_diff_box{vm_id}{box_id}(1,3);      
    end
    % Record the corresponding VM ID
    vm_idx(vm_id) = vm_fill_gaps{vm_id}(1,1);
end

% Merge the DISK into the summary VM
% for vm_id = 1 : numel(disk_fill_gaps)
%     vm_id_in_vm_idx = find(vm_idx == disk_fill_gaps{vm_id}(1,1));
%     if numel(vm_id_in_vm_idx) == 0
%         continue
%     end
%     for box_id = 1 : numel(summary_of_vm_diff_box{vm_id_in_vm_idx})
%         [common_time, vm_series, disk_series] = ...
%         intersect(summary_of_vm_diff_box{vm_id_in_vm_idx}{box_id}(:,2), ...
%                       disk_fill_gaps{vm_id}(:,2));
%         if numel(common_time) ~= 0
%             summary_of_vm_diff_box{vm_id_in_vm_idx}{box_id}(vm_series, 6) = ...
%                       disk_fill_gaps{vm_id}(disk_series, 3);
%         end
%     end
% end

% Merge the NET into the summary VM
% for vm_id = 1 : numel(net_fill_gaps)
%     vm_id_in_vm_idx = find(vm_idx == net_fill_gaps{vm_id}(1,1));
%     if numel(vm_id_in_vm_idx) == 0
%         continue
%     end
%     for box_id = 1 : numel(summary_of_vm_diff_box{vm_id_in_vm_idx})
%         [common_time, vm_series, disk_series] = ...
%         intersect(summary_of_vm_diff_box{vm_id_in_vm_idx}{box_id}(:,2), ...
%                       net_fill_gaps{vm_id}(:,2));
%         if numel(common_time) ~= 0
%             summary_of_vm_diff_box{vm_id_in_vm_idx}{box_id}(vm_series, 7:8) = ...
%                       net_fill_gaps{vm_id}(disk_series, 3:4);
%         end
%     end
% end

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

save('../box_vm_time_series','box_vm_time_series')

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
    
    common_time = box_vm_time_series{box_id}{1}(:,2);
    for cand_id = 2 : numel(box_vm_time_series{box_id})
        common_time = intersect(common_time, box_vm_time_series{box_id}{cand_id}(:,2));
    end
    
    if numel(common_time) == 0
        disp('This box has no VM concurrently residing on it');
        continue;
    end
    
    common_time_idx = {}; 
    for cand_id = 1 : numel(box_vm_time_series{box_id})
        [common_time_series, common_idx, common_time_idx{cand_id}] = ...
        intersect(common_time, box_vm_time_series{box_id}{cand_id}(:,2));
        box_vm_time_series_summary{box_id}{cand_id} = ...
            box_vm_time_series{box_id}{cand_id}(common_time_idx{cand_id},:);
        
        % Change the time series into tickets
        all_metric_no = numel(box_vm_time_series_summary{box_id}{cand_id}(1,:));
        box_vm_time_series_ticket_summary{box_id}{cand_id}(:,1:all_metric_no-2) = ...
            box_vm_time_series_summary{box_id}{cand_id}(:,1:all_metric_no-2);
        box_vm_time_series_ticket_summary{box_id}{cand_id}(:, all_metric_no-1) = ...
            box_vm_time_series_summary{box_id}{cand_id}(:,all_metric_no-1) >= cpu_ticket_thres;
        box_vm_time_series_ticket_summary{box_id}{cand_id}(:, all_metric_no) = ...
            box_vm_time_series_summary{box_id}{cand_id}(:,all_metric_no) >= mem_ticket_thres;         
    end
    
end

save('../box_vm_time_series_summary','box_vm_time_series_summary');
save('../box_vm_time_series_ticket_summary', 'box_vm_time_series_ticket_summary');

%% Step 5: Check the correlation among different time series
corr_different_series = {}; corr_different_series_ticket = {};
max_lag = 2; metric_vm_no = 2;
for box_id = 1 : numel(box_vm_time_series_summary)
    machine_num = numel(box_vm_time_series_summary{box_id});
    if machine_num == 0
        continue;
    end
    
    start_time = 1; time_id = 1; slot_num = 1;
    len = numel(box_vm_time_series_summary{box_id}{1,1}(:,1));
    while time_id <= len
        if box_vm_time_series_summary{box_id}{1,1}(time_id, 3) == 0
            end_time = time_id - 1;
            
            % If too few samples
            if end_time - start_time <= 4*3600/grat_small %2* max_lag
                % Update the start time
                while time_id <= len
                    if box_vm_time_series_summary{box_id}{1,1}(time_id, 3) ~= 0
                        break;
                    end
                    time_id = time_id + 1;
                end
                start_time = time_id;
                continue
            end
            
            corr_different_series{box_id}{1,slot_num} = {};
            corr_different_series_ticket{box_id}{1, slot_num} = {};
            corr_different_series{box_id}{2,slot_num} = [start_time, end_time];
            corr_different_series_ticket{box_id}{2, slot_num} = [start_time, end_time];          
            
            for machine_id = 1 : machine_num
                % Find the time series index
                if machine_id == 1
                    begin_id_itself = 3;
                    corr_id_row = 1;
                else
                    begin_id_itself = 4;
                    corr_id_row = 3 + (machine_id-2) * metric_vm_no;
                end
                end_id_itself = numel(box_vm_time_series_summary{box_id}{machine_id}(1,:));
                
                for other_machine_id = 1 : machine_num
                    % Find the compared time series index
                    if other_machine_id == 1
                        begin_id_other = 3;
                        corr_id_col = 1;
                    else
                        begin_id_other = 4;
                        corr_id_col = 3 + (other_machine_id-2) * metric_vm_no;
                    end
                    end_id_other = numel(box_vm_time_series_summary{box_id}{other_machine_id}(1,:));
                    
                    % Two for loops to calculate the XCF
                    for time_series_itself_id = begin_id_itself : end_id_itself
                        time_series_itself = box_vm_time_series_summary...
                            {box_id}{machine_id}(start_time:end_time, time_series_itself_id);
                        time_series_itself_ticket = box_vm_time_series_ticket_summary ...
                            {box_id}{machine_id}(start_time:end_time, time_series_itself_id);
                        for time_series_other_id = begin_id_other : end_id_other
                            time_series_other = box_vm_time_series_summary...
                            {box_id}{other_machine_id}(start_time:end_time, time_series_other_id);
                            time_series_other_ticket = box_vm_time_series_ticket_summary ...
                            {box_id}{other_machine_id}(start_time:end_time, time_series_other_id);
                            
                            [xcf, lags, bounds] = crosscorr(time_series_itself, time_series_other, max_lag);
                            [xcf_ticket, lags, bounds] = crosscorr(time_series_itself_ticket, time_series_other_ticket, max_lag);                         
                            
                            row_id = time_series_itself_id - begin_id_itself + corr_id_row;
                            col_id = time_series_other_id - begin_id_other + corr_id_col;
                            for lag_idx = 1 : 2 * max_lag + 1                              
                                corr_different_series{box_id}{1, slot_num}...
                                    {lag_idx}(row_id, col_id)...
                                    = xcf(lag_idx);
                                corr_different_series_ticket{box_id}{1, slot_num}...
                                    {lag_idx}(row_id, col_id)...
                                    = xcf_ticket(lag_idx);
                            end                           
                        end
                    end
                    
                end
                
            end
            
            % Update the start time
            while time_id <= len
                if box_vm_time_series_summary{box_id}{1,1}(time_id, 3) ~= 0
                    break;
                end
                time_id = time_id + 1;
            end
            start_time = time_id;
            time_id = time_id - 1;
            
            % Update the slot
            slot_num = slot_num + 1;
        end
        time_id = time_id + 1;
    end
    
end

save('../corr_different_series','corr_different_series');
save('../corr_different_series_ticket','corr_different_series_ticket');

%% Step 6: Plot the heat map to check cross-correlation
mkdir('../XCF_Small_Actual'); mkdir('../XCF_Small_Ticket');
path1 = '../XCF_Small_Actual/'; path2 = '../XCF_Small_Ticket/';
% BOX - VM:  CPU v.s. CPU, MEM v.s. MEM
% BOX - BOX: CPU v.s. MEM
% VM - VM:   CPU v.s. MEM
all_corr_summary = {[], [], [], []}; 
all_corr_ticket_summary = {[], [], [], []};
for box_id = 1 : numel(corr_different_series)
    if numel(corr_different_series{box_id}) == 0
        continue
    end
    box_labels = {'BOX-CPU','BOX-MEM'};
    num_vm = numel(box_vm_time_series_summary{box_id}) -1;
    cpu_index = 1 : 2 : 2 * (num_vm + 1);
    mem_index = 2 : 2 : 2 * (num_vm + 1);
    labels = box_labels;
    labels_cpu = {box_labels{1}}; labels_mem = {box_labels{2}};
    for vm_id = 1 : num_vm
        vm_label_no = {strcat('VM',mat2str(vm_id),'-CPU'), ...
                       strcat('VM',mat2str(vm_id),'-MEM')};
        labels = {labels{:}, vm_label_no{:}};
        labels_cpu = {labels_cpu{:}, vm_label_no{1}};
        labels_mem = {labels_mem{:}, vm_label_no{2}};
    end
    for time_slot = 1 : numel(corr_different_series{box_id}(1,:))
   
        for lag_id = 1 : max_lag + 1
            % Change Nans to zeros
            nan_index = isnan(corr_different_series{box_id}{1, time_slot}{lag_id});
            nan_index_ticket = isnan(corr_different_series_ticket{box_id}{1, time_slot}{lag_id});
            len_metric = sqrt(numel(corr_different_series_ticket{box_id}{1, time_slot}{lag_id}));
            corr_different_series{box_id}{1, time_slot}{lag_id}(nan_index) = 0.001 * rand([sum(sum(nan_index)),1]);           
            corr_different_series_ticket{box_id}{1, time_slot}{lag_id}(nan_index_ticket) = 0.001 * rand([sum(sum(nan_index_ticket)),1]);
            
            if lag_id == max_lag + 1
               % Fill in the all correlation summary
               all_corr_summary{1} = [all_corr_summary{1}, corr_different_series{box_id}{1, time_slot}{lag_id}(1, cpu_index(2:end))];
               all_corr_summary{2} = [all_corr_summary{2}, corr_different_series{box_id}{1, time_slot}{lag_id}(2, mem_index(2:end))];
               all_corr_summary{3} = [all_corr_summary{3}, corr_different_series{box_id}{1, time_slot}{lag_id}(1, 2)];
               for cpu_id = 2 : num_vm + 1
                   all_corr_summary{4} = [all_corr_summary{4}, corr_different_series{box_id}{1, time_slot}{lag_id}(cpu_index(cpu_id), cpu_index(cpu_id)+1)]; 
               end
            end
            
            % In case of all the XCF are same
            % The following formula needs to be added if we plot heatmap
            corr_different_series_ticket{box_id}{1, time_slot}{lag_id} = ...
                corr_different_series_ticket{box_id}{1, time_slot}{lag_id} + 0.001 * rand([len_metric, len_metric]);
            
%             % Plot all the XCF for CPU + MEM
%             fig = figure;
%             set(fig,'Position',[200,200,1000, 800]);
%             heatmap(corr_different_series{box_id}{1, time_slot}{lag_id}, labels, labels, '%0.2f', 'TextColor', 'w', ...
%                     'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
%                     'TickAngle', 45, 'TickFontSize',15, 'FontSize', 20);
%             caxis([-1 1]);
%             time_len = corr_different_series{box_id}{2, time_slot}(2) - ...
%                        corr_different_series{box_id}{2, time_slot}(1);
%             title_name = strcat('BOX ID=', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)), ...
%                 ', Time Length is=', mat2str(int32(time_len*grat_small/3600)), 'h, Lag =', mat2str((max_lag-lag_id+1)*grat_small/60), ... 
%                 ' min');
%             title(title_name);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path1, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(corr_different_series{box_id}{2, time_slot}(1)), ...
%                    '_Length_', mat2str(int32(time_len*grat_small/3600)), '_Lag_', mat2str(max_lag - lag_id + 1), '_CPU+MEM'));
%             
%             % Plot XCF only for CPU
%             fig = figure;
%             set(fig,'Position',[200,200,1000, 800]);
%             heatmap(corr_different_series{box_id}{1, time_slot}{lag_id}(cpu_index, cpu_index), ...
%                     labels_cpu, labels_cpu, '%0.2f', 'TextColor', 'w', ...
%                     'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
%                     'TickAngle', 45, 'TickFontSize',15, 'FontSize', 20);
%             caxis([-1 1]);
%             title_name = strcat('BOX ID=', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)), ...
%                 ', CPU USD PCT, Time Length is=', mat2str(int32(time_len*grat_small/3600)), ...
%                 'h, Lag =', mat2str((max_lag-lag_id+1)*grat_small/60), ' min');
%             title(title_name);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path1, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(corr_different_series{box_id}{2, time_slot}(1)), ...
%                    '_Length_', mat2str(int32(time_len*grat_small/3600)), '_Lag_', mat2str(max_lag - lag_id + 1), '_CPU'));
%                
%             % Plot XCF only for MEM
%             fig = figure;
%             set(fig,'Position',[200,200,1000, 800]);
%             heatmap(corr_different_series{box_id}{1, time_slot}{lag_id}(mem_index, mem_index), ...
%                     labels_mem, labels_mem, '%0.2f', 'TextColor', 'w', ...
%                     'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
%                     'TickAngle', 45, 'TickFontSize',15, 'FontSize', 20);
%             caxis([-1 1]);
%             title_name = strcat('BOX ID=', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)), ...
%                 ', MEM USED PCT, Time Length is=', mat2str(int32(time_len*grat_small/3600)), ...
%                 'h, Lag =', mat2str((max_lag-lag_id+1)*grat_small/60), ' min');
%             title(title_name);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path1, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(corr_different_series{box_id}{2, time_slot}(1)), ...
%                    '_Length_', mat2str(int32(time_len*grat_small/3600)), '_Lag_', mat2str(max_lag - lag_id + 1), '_MEM'));
%              
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%             % Plot CPU + MEM
%             fig = figure;
%             set(fig,'Position',[200,200,1000, 800]);
%             heatmap(corr_different_series_ticket{box_id}{1, time_slot}{lag_id}, labels, labels, '%0.2f', 'TextColor', 'w', ...
%                     'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
%                     'TickAngle', 45, 'TickFontSize',15, 'FontSize', 20);
%             caxis([-1 1]);
%             title(title_name);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path2, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(corr_different_series{box_id}{2, time_slot}(1)), ...
%                    '_Length_', mat2str(int32(time_len*grat_small/3600)), '_Lag_', mat2str(max_lag - lag_id + 1), '_CPU+MEM'));
%                
%             % Plot XCF only for CPU
%             fig = figure;
%             set(fig,'Position',[200,200,1000, 800]);
%             heatmap(corr_different_series_ticket{box_id}{1, time_slot}{lag_id}(cpu_index, cpu_index), ...
%                     labels_cpu, labels_cpu, '%0.2f', 'TextColor', 'w', ...
%                     'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
%                     'TickAngle', 45, 'TickFontSize',15, 'FontSize', 20);
%             caxis([-1 1]);
%             title_name = strcat('BOX ID=', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)), ...
%                 ', CPU USD PCT, Time Length is=', mat2str(int32(time_len*grat_small/3600)), ...
%                 'h, Lag =', mat2str((max_lag-lag_id+1)*grat_small/60), ' min');
%             title(title_name);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path2, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(corr_different_series{box_id}{2, time_slot}(1)), ...
%                    '_Length_', mat2str(int32(time_len*grat_small/3600)), '_Lag_', mat2str(max_lag - lag_id + 1), '_CPU'));
%                
%             % Plot XCF only for MEM
%             fig = figure;
%             set(fig,'Position',[200,200,1000, 800]);
%             heatmap(corr_different_series_ticket{box_id}{1, time_slot}{lag_id}(mem_index, mem_index), ...
%                     labels_mem, labels_mem, '%0.2f', 'TextColor', 'w', ...
%                     'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
%                     'TickAngle', 45, 'TickFontSize',15, 'FontSize', 20);
%             caxis([-1 1]);
%             title_name = strcat('BOX ID=', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)), ...
%                 ', MEM USED PCT, Time Length is=', mat2str(int32(time_len*grat_small/3600)), ...
%                 'h, Lag =', mat2str((max_lag-lag_id+1)*grat_small/60), ' min');
%             title(title_name);
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path2, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(corr_different_series{box_id}{2, time_slot}(1)), ...
%                    '_Length_', mat2str(int32(time_len*grat_small/3600)), '_Lag_', mat2str(max_lag - lag_id + 1), '_MEM'));
%             
%             close all
        end
    end
end

% Plot the CDF of the XCF summary
fig = figure;
set(fig,'Position',[200,200,800, 600]);
[f_1, x_1] = ecdf(all_corr_summary{1});
[f_2, x_2] = ecdf(all_corr_summary{2});
[f_3, x_3] = ecdf(all_corr_summary{3});
[f_4, x_4] = ecdf(all_corr_summary{4});
plot(x_1, f_1, 'r--', 'linewidth', 2);
hold on
plot(x_2, f_2, 'b-', 'linewidth', 2);
hold on 
plot(x_3, f_3, 'k-.', 'linewidth', 2);
hold on
plot(x_4, f_4, 'm:', 'linewidth', 2);

set(gca, 'fontsize', 18);
h = legend('BOX-VM: CPU PCT V.S CPU PCT', 'BOX-VM: MEM PCT V.S. MEM PCT', ...
           'BOX-BOX: CPU PCT V.S. MEM PCT', 'VM-VM: CPU PCT V.S. MEM PCT');
set(h, 'location', 'northwest', 'box','on');
xlabel('XCF (Lag = 0 min)'); ylabel('CDF')
set(gcf, 'paperpositionmode', 'auto');
set(gca, 'xlim', [-1 1]);
print('-depsc2','-r300', '../CPU_MEM_ALL_Figure/XCF_CDF');

% Plot the CDF of the (absolute) XCF summary
fig = figure;
set(fig,'Position',[200,200,800, 600]);
[f_1, x_1] = ecdf(abs(all_corr_summary{1}));
[f_2, x_2] = ecdf(abs(all_corr_summary{2}));
[f_3, x_3] = ecdf(abs(all_corr_summary{3}));
[f_4, x_4] = ecdf(abs(all_corr_summary{4}));
plot(x_1, f_1, 'r--', 'linewidth', 2);
hold on
plot(x_2, f_2, 'b-', 'linewidth', 2);
hold on 
plot(x_3, f_3, 'k-.', 'linewidth', 2);
hold on
plot(x_4, f_4, 'm:', 'linewidth', 2);

set(gca, 'fontsize', 18);
h = legend('BOX-VM: CPU PCT V.S CPU PCT', 'BOX-VM: MEM PCT V.S. MEM PCT', ...
           'BOX-BOX: CPU PCT V.S. MEM PCT', 'VM-VM: CPU PCT V.S. MEM PCT');
set(h, 'location', 'southeast', 'box','on');
xlabel('XCF (Lag = 0 min)'); ylabel('CDF')
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', '../CPU_MEM_ALL_Figure/XCF_CDF_Positive');

% Plot the PDF of the XCF summary
num_bin = 20;
[n_1, x_1] = pdfplot(all_corr_summary{1}, num_bin, -1, 1);
[n_2, x_2] = pdfplot(all_corr_summary{2}, num_bin, -1, 1);
[n_3, x_3] = pdfplot(all_corr_summary{3}, num_bin, -1, 1);
[n_4, x_4] = pdfplot(all_corr_summary{4}, num_bin, -1, 1);

fig = figure;
set(fig,'Position',[200,200,800, 600]);
plot(x_1, n_1/sum(n_1), 'r--', 'linewidth', 2);
hold on
plot(x_2, n_2/sum(n_2), 'b-', 'linewidth', 2);
hold on 
plot(x_3, n_3/sum(n_3), 'k-.', 'linewidth', 2);
hold on
plot(x_4, n_4/sum(n_4), 'm:', 'linewidth', 2);
set(gca, 'fontsize', 18);
h = legend('BOX-VM: CPU PCT V.S CPU PCT', 'BOX-VM: MEM PCT V.S. MEM PCT', ...
           'BOX-BOX: CPU PCT V.S. MEM PCT', 'VM-VM: CPU PCT V.S. MEM PCT');
set(h, 'location', 'northwest', 'box','on');
xlabel('XCF (Lag = 0 min)'); ylabel('PDF')
set(gcf, 'paperpositionmode', 'auto');
set(gca, 'xlim', [-1 1]);
print('-depsc2','-r300', '../CPU_MEM_ALL_Figure/XCF_PDF');

%% Step 7: Check the time series of tickets or actual
% Plot the time series of CPU and Memory for BOX and VMs 
% mkdir('../Time_Series_Actual'); mkdir('../Time_Series_Ticket');
% path1 = '../Time_Series_Actual/'; path2 = '../Time_Series_Ticket/';
% for box_id = 1 : numel(box_vm_time_series_summary)
%     machine_num = numel(box_vm_time_series_summary{box_id});
%     if machine_num == 0
%         continue
%     end
%     
%     start_time = 1; time_id = 1; slot_num = 1;
%     len = numel(box_vm_time_series_summary{box_id}{1,1}(:,1));
%     while time_id <= len
%         if box_vm_time_series_summary{box_id}{1,1}(time_id, 3) == 0
%             end_time = time_id - 1;
%             
%             % If too few samples
%             if end_time - start_time <= 4 * 3600 / grat_small
%                 % Update the start time
%                 while time_id <= len
%                     if box_vm_time_series_summary{box_id}{1,1}(time_id, 3) ~= 0
%                         break;
%                     end
%                     time_id = time_id + 1;
%                 end
%                 start_time = time_id;
%                 continue
%             end
%             
%             % Plot Actual Overtime
%             fig = figure;
%             set(fig,'Position',[200, 200, 2000, 2000]);
%             set(gca, 'fontsize', 15);
%             time_len = box_vm_time_series_summary{box_id}{1,1}(start_time:end_time, 2)...
%                        / 3600;
% 
%             subplot(machine_num, 2, 1)
%             plot(time_len, box_vm_time_series_summary{box_id}{1,1}(start_time:end_time,3),'r-');
%             max_y = ceil(max(box_vm_time_series_summary{box_id}{1,1}(start_time:end_time,3))/10)*10;
%             h = legend('BOX');
%             set(h,'location','northeast','box','off')
%             title('CPU USE PCT');
%             set(gca, 'ylim', [0 max_y]);
%             set(gca, 'ytick', [0 : max_y/4 : max_y]);
%             xlabel('Time (hour)'); ylabel('CPU USED PCT')
% 
%             subplot(machine_num, 2, 2)
%             plot(time_len, box_vm_time_series_summary{box_id}{1,1}(start_time:end_time,4),'b-');
%             max_y = ceil(max(box_vm_time_series_summary{box_id}{1,1}(start_time:end_time,4))/10)*10;
%             h = legend('BOX');
%             set(h,'location','northeast','box','off')
%             title('MEM USE PCT');
%             set(gca, 'ylim', [0 max_y]);
%             set(gca, 'ytick', [0 : max_y/4 : max_y]);
%             xlabel('Time (hour)'); ylabel('MEM USED PCT')
% 
%             for machine_id = 2 : machine_num
%                 subplot(machine_num, 2, (machine_id -1)*2+1)
%                 plot(time_len, box_vm_time_series_summary{box_id}{1,machine_id}(start_time:end_time,4),'r-');
%                 max_y = ceil(max(box_vm_time_series_summary{box_id}{1,machine_id}(start_time:end_time,4))/10)*10;
%                 h = legend(strcat('VM', mat2str(machine_id -1)));
%                 set(h,'location','northeast','box','off')
%                 set(gca, 'ylim', [0 max_y]);
%                 set(gca, 'ytick', [0 : max_y/4 : max_y]);
%                 xlabel('Time (hour)'); ylabel('CPU USED PCT')
% 
%                 subplot(machine_num, 2, (machine_id -1)*2+2)
%                 plot(time_len, box_vm_time_series_summary{box_id}{1,machine_id}(start_time:end_time,5),'b-');
%                 max_y = ceil(max(box_vm_time_series_summary{box_id}{1,machine_id}(start_time:end_time,5))/10)*10;
%                 h = legend(strcat('VM', mat2str(machine_id -1)));
%                 set(h,'location','northeast','box','off')
%                 set(gca, 'ylim', [0 max_y]);
%                 set(gca, 'ytick', [0 : max_y/4 : max_y]);
%                 xlabel('Time (hour)'); ylabel('MEM USED PCT')               
%             end
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path1, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(start_time), ...
%                    '_Length_', mat2str(int32((end_time-start_time)*grat_small/3600)), '_Actual'));
% 
%             % Plot Ticket Overtime
%             fig = figure;
%             set(fig,'Position',[200, 200, 2000, 2000]);
%             set(gca, 'fontsize', 15);
%             subplot(machine_num, 2, 1)
%             plot(time_len, box_vm_time_series_ticket_summary{box_id}{1,1}(start_time:end_time,3),'r-');
%             h = legend('BOX');
%             set(h,'location','northeast','box','off')
%             title('CPU USE PCT (Ticket(1) or Not(0))');
%             set(gca, 'ylim', [0 2]);
%             set(gca, 'ytick', [0 : 1 : 2]);
%             xlabel('Time (hour)'); ylabel('CPU USED PCT')
% 
%             subplot(machine_num, 2, 2)
%             plot(time_len, box_vm_time_series_ticket_summary{box_id}{1,1}(start_time:end_time,4),'b-');
%             h = legend('BOX');
%             set(h,'location','northeast','box','off')
%             title('MEM USE PCT (Ticket(1) or Not(0))');
%             set(gca, 'ylim', [0 2]);
%             set(gca, 'ytick', [0 : 1: 2]);
%             xlabel('Time (hour)'); ylabel('MEM USED PCT')
% 
%             for machine_id = 2 : machine_num
%                 subplot(machine_num, 2, (machine_id -1)*2+1)
%                 plot(time_len, box_vm_time_series_ticket_summary{box_id}{1,machine_id}(start_time:end_time,4),'r-');
%                 h = legend(strcat('VM', mat2str(machine_id -1)));
%                 set(h,'location','northeast','box','off')
%                 set(gca, 'ylim', [0 2]);
%                 set(gca, 'ytick', [0 : 1 : 2]);
%                 xlabel('Time (hour)'); ylabel('CPU USED PCT')
% 
%                 subplot(machine_num, 2, (machine_id -1)*2+2)
%                 plot(time_len, box_vm_time_series_ticket_summary{box_id}{1,machine_id}(start_time:end_time,5),'b-');
%                 h = legend(strcat('VM', mat2str(machine_id -1)));
%                 set(h,'location','northeast','box','off')
%                 set(gca, 'ylim', [0 2]);
%                 set(gca, 'ytick', [0 : 1 : 2]);
%                 xlabel('Time (hour)'); ylabel('MEM USED PCT')               
%             end
%             set(gcf, 'paperpositionmode', 'auto');
%             print('-depsc2','-r300', strcat(path2, 'PM_', mat2str(box_vm_time_series_summary{box_id}{1}(1,1)),...
%                    '_Start_', mat2str(start_time), ...
%                    '_Length_', mat2str(int32((end_time-start_time)*grat_small/3600)), '_Ticket'));
%             close all
%             
%             % Update the start time
%             while time_id <= len
%                 if box_vm_time_series_summary{box_id}{1,1}(time_id, 3) ~= 0
%                     break;
%                 end
%                 time_id = time_id + 1;
%             end
%             start_time = time_id;
%             time_id = time_id - 1;
%             
%             % Update the slot
%             slot_num = slot_num + 1;
%         end
%         time_id = time_id + 1;
%     end
% end      

