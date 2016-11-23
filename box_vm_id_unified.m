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
load ../net_differentID
% Grat = 60 min; NET_TX_MBPS	NET_RX_MBPS
net = machine_cell;
load ../vm_differentID
% Grat = 15 min; CPU_USED_PCT	CPU_USED_MHZ	CPU_COUNT	MEM_USED_PCT	MEM_USED_MB	
% VM contains the BOX information, so we also need to divide based on BOX
vm = machine_cell;

%% Step 1: Change the granularity for all the metrics
grat_small = 900; grat_big = 3600;
box_fill_gaps = {}; box_fill_gaps_big = {};
disk_fill_gaps = {}; disk_fill_gaps_big = {};
net_fill_gaps_big = {}; 
vm_fill_gaps = {}; vm_fill_gaps_big = {};
for box_id = 1 : numel(box)
   [box_fill_gaps{box_id}] = fill_gaps(box{box_id}, grat_small);
   [box_fill_gaps_big{box_id}] = change_grat(box_fill_gaps{box_id}, grat_small, grat_big, false);
end

for vm_id = 1 : numel(disk)
   [disk_fill_gaps{vm_id}] = fill_gaps(disk{vm_id}, grat_small);
   [disk_fill_gaps_big{vm_id}] = change_grat(disk_fill_gaps{vm_id}, grat_small, grat_big, false);
end

for vm_id = 1 : numel(net)
   [net_fill_gaps_big{vm_id}] = fill_gaps(net{vm_id}, grat_big);
end

for vm_id = 1 : numel(vm)
   [vm_fill_gaps{vm_id}] = fill_gaps(vm{vm_id}, grat_small);
   [vm_fill_gaps_big{vm_id}] = change_grat(vm_fill_gaps{vm_id}, grat_small, grat_big, true);
end

%% Step 2: Group the different VMs into the same PM (with big granularity)
% The structure is: PM, VM0, VM1, ...
% The structure for PM is: PM_id, Time, CPU_PCT, MEM_PCT
% The structure for VM is: VM_id, Time, PM_id, CPU_PCT, MEM_PCT, DISK_PCT,
% NET_TX, NET_RX

box_vm_time_series_big = {}; box_idx = [];
% First put the PM into each cell and record the corresponding BOX ID
for box_id = 1 : numel(box_fill_gaps_big)
    box_vm_time_series_big{box_id}{1} = box_fill_gaps_big{box_id}(:,1:3);
    box_vm_time_series_big{box_id}{1}(:, 4) = box_fill_gaps_big{box_id}(:, 6);
    box_idx(box_id) = box_fill_gaps_big{box_id}(1,1);
end

% Second for each VM, seperate by VM ID and PM ID
% Each VM could be allocated to several PM
vm_fill_gaps_big_diff_box = {};
vm_idx = []; vm_box_idx = {};
summary_of_vm_big_diff_box = {};
for vm_id = 1 : numel(vm_fill_gaps_big)
    no_box = 1; time_id = 1; box_id = vm_fill_gaps_big{vm_id}(1,3);
    start_time_id = 1; end_time_id = numel(vm_fill_gaps_big{vm_id}(:,1));
    while time_id <= end_time_id
        if (box_id ~= vm_fill_gaps_big{vm_id}(time_id,3) ...
            && vm_fill_gaps_big{vm_id}(time_id,3) ~= 0) || time_id == end_time_id
            vm_fill_gaps_big_diff_box{vm_id}{no_box} = vm_fill_gaps_big{vm_id}...
                                                   (start_time_id:time_id-1,:);
            start_time_id = time_id; 
            box_id = vm_fill_gaps_big{vm_id}(time_id,3); no_box = no_box + 1;
        end
        time_id = time_id + 1;
    end
    vm_fill_gaps_big_diff_box{vm_id}{no_box-1} = ...
    [vm_fill_gaps_big_diff_box{vm_id}{no_box-1}; vm_fill_gaps_big{vm_id}(time_id-1, :)];
    
    for box_id = 1 : no_box - 1
        % Pick up the needed performance metrics
        summary_of_vm_big_diff_box{vm_id}{box_id} = ...
            vm_fill_gaps_big_diff_box{vm_id}{box_id}(:,1:4);
        summary_of_vm_big_diff_box{vm_id}{box_id}(:,5) = ...
            vm_fill_gaps_big_diff_box{vm_id}{box_id}(:,7);
        
        % Pre-fill the empty fields for disk and net
        summary_of_vm_big_diff_box{vm_id}{box_id}(:, 6:8) = 0;
        
        % Record the corresponding BOX ID
        vm_box_idx{vm_id}(box_id) = vm_fill_gaps_big_diff_box{vm_id}{box_id}(1,3);      
    end
    % Record the corresponding VM ID
    vm_idx(vm_id) = vm_fill_gaps{vm_id}(1,1);
end

% Merge the DISK into the summary VM
for vm_id = 1 : numel(disk_fill_gaps_big)
    vm_id_in_vm_idx = find(vm_idx == disk_fill_gaps_big{vm_id}(1,1));
    if numel(vm_id_in_vm_idx) == 0
        continue
    end
    for box_id = 1 : numel(summary_of_vm_big_diff_box{vm_id_in_vm_idx})
        [common_time, vm_series, disk_series] = ...
        intersect(summary_of_vm_big_diff_box{vm_id_in_vm_idx}{box_id}(:,2), ...
                      disk_fill_gaps_big{vm_id}(:,2));
        if numel(common_time) ~= 0
            summary_of_vm_big_diff_box{vm_id_in_vm_idx}{box_id}(vm_series, 6) = ...
                      disk_fill_gaps_big{vm_id}(disk_series, 3);
        end
    end
end

% Merge the NET into the summary VM
for vm_id = 1 : numel(net_fill_gaps_big)
    vm_id_in_vm_idx = find(vm_idx == net_fill_gaps_big{vm_id}(1,1));
    if numel(vm_id_in_vm_idx) == 0
        continue
    end
    for box_id = 1 : numel(summary_of_vm_big_diff_box{vm_id_in_vm_idx})
        [common_time, vm_series, disk_series] = ...
        intersect(summary_of_vm_big_diff_box{vm_id_in_vm_idx}{box_id}(:,2), ...
                      net_fill_gaps_big{vm_id}(:,2));
        if numel(common_time) ~= 0
            summary_of_vm_big_diff_box{vm_id_in_vm_idx}{box_id}(vm_series, 7:8) = ...
                      net_fill_gaps_big{vm_id}(disk_series, 3:4);
        end
    end
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
            box_vm_time_series_big{box_id_final}{end+1} = summary_of_vm_big_diff_box...
                                                          {vm_id}{box_id};
        end
    end    
end

save('../box_vm_time_series_big','box_vm_time_series_big')

%% Step 4: Combine VM and PM into one *big* time series
% This part makes the VMs and PM share the same time stamps to calculate
% the correlation among the big time series. The first try is to get some
% examples that all the VM have same resindual time.
box_vm_time_series_big_summary = {};
for box_id = 1 : numel(box_vm_time_series_big)
    if numel(box_vm_time_series_big{box_id}) == 1
        disp('This box has no VM residing on it');
        continue;
    end
    
    common_time = box_vm_time_series_big{box_id}{1}(:,2);
    for cand_id = 2 : numel(box_vm_time_series_big{box_id})
        common_time = intersect(common_time, box_vm_time_series_big{box_id}{cand_id}(:,2));
    end
    
    if numel(common_time) == 0
        disp('This box has no VM concurrently residing on it');
        continue;
    end
    
    common_time_idx = {}; 
    for cand_id = 1 : numel(box_vm_time_series_big{box_id})
        [common_time_series, common_idx, common_time_idx{cand_id}] = ...
        intersect(common_time, box_vm_time_series_big{box_id}{cand_id}(:,2));
        box_vm_time_series_big_summary{box_id}{cand_id} = ...
            box_vm_time_series_big{box_id}{cand_id}(common_time_idx{cand_id},:);
    end
    
end

save('../box_vm_time_series_big_summary','box_vm_time_series_big_summary');

%% Step 5: Check the correlation among different time series
corr_different_series = {}; max_lag = 2; metric_vm_no = 5;
for box_id = 1 : numel(box_vm_time_series_big_summary)
    machine_num = numel(box_vm_time_series_big_summary{box_id});
    if machine_num == 0
        continue;
    end
    
    start_time = 1; time_id = 1; slot_num = 1;
    len = numel(box_vm_time_series_big_summary{box_id}{1,1}(:,1));
    while time_id <= len
        if box_vm_time_series_big_summary{box_id}{1,1}(time_id, 3) == 0
            end_time = time_id - 1;
            
            % If too few samples
            if end_time - start_time <= 2* max_lag
                % Update the start time
                while time_id <= len
                    if box_vm_time_series_big_summary{box_id}{1,1}(time_id, 3) ~= 0
                        break;
                    end
                    time_id = time_id + 1;
                end
                start_time = time_id;
                continue
            end
            
            corr_different_series{box_id}{1,slot_num} = {};
            corr_different_series{box_id}{2,slot_num} = [start_time, end_time];
            for machine_id = 1 : machine_num
                % Find the time series index
                if machine_id == 1
                    begin_id_itself = 3;
                    corr_id_row = 1;
                else
                    begin_id_itself = 4;
                    corr_id_row = 3 + (machine_id-2) * metric_vm_no;
                end
                end_id_itself = numel(box_vm_time_series_big_summary{box_id}{machine_id}(1,:));
                
                for other_machine_id = 1 : machine_num
                    % Find the compared time series index
                    if other_machine_id == 1
                        begin_id_other = 3;
                        corr_id_col = 1;
                    else
                        begin_id_other = 4;
                        corr_id_col = 3 + (other_machine_id-2) * metric_vm_no;
                    end
                    end_id_other = numel(box_vm_time_series_big_summary{box_id}{other_machine_id}(1,:));
                    
                    % Two for loops to calculate the XCF
                    for time_series_itself_id = begin_id_itself : end_id_itself
                        time_series_itself = box_vm_time_series_big_summary...
                            {box_id}{machine_id}(start_time:end_time, time_series_itself_id);
                        for time_series_other_id = begin_id_other : end_id_other
                            time_series_other = box_vm_time_series_big_summary...
                            {box_id}{other_machine_id}(start_time:end_time, time_series_other_id);
                            [xcf, lags, bounds] = crosscorr(time_series_itself, time_series_other, max_lag);
                            row_id = time_series_itself_id - begin_id_itself + corr_id_row;
                            col_id = time_series_other_id - begin_id_other + corr_id_col;
                            for lag_idx = 1 : 2 * max_lag + 1                              
                                corr_different_series{box_id}{1, slot_num}...
                                    {lag_idx}(row_id, col_id)...
                                    = xcf(lag_idx);
                            end
                        end
                    end
                    
                end
                
            end
            
            % Update the start time
            while time_id <= len
                if box_vm_time_series_big_summary{box_id}{1,1}(time_id, 3) ~= 0
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

%% Step 6: Plot the heat map to check cross-correlation
mkdir('../XCF');
path = '../XCF/';
for box_id = 1 : numel(corr_different_series)
    if numel(corr_different_series{box_id}) == 0
        continue
    end
    box_labels = {'BOX-CPU','BOX-MEM'};
    vm_labels = {'VM-CPU','VM-MEM','VM-DISK','VM-NET-TX','VM-NET-RX'};
    num_vm = numel(box_vm_time_series_big_summary{box_id}) -1;
    labels = box_labels;
    for vm_id = 1 : num_vm
        labels = {labels{:}, vm_labels{:}};
    end
    for time_slot = 1 : numel(corr_different_series{box_id}(1,:))
        for lag_id = 1 : max_lag + 1
            fig = figure;
            set(fig,'Position',[200,200,1000, 800]);
            heatmap(corr_different_series{box_id}{1, time_slot}{lag_id}, labels, labels, '%0.2f', 'TextColor', 'w', ...
                    'Colormap', 'copper', 'Colorbar', true, 'ShowAllTicks',true,...
                    'TickAngle', 45, 'TickFontSize',10);
            caxis([-1 1]);
            time_len = corr_different_series{box_id}{2, time_slot}(2) - ...
                       corr_different_series{box_id}{2, time_slot}(1);
            title_name = strcat('BOX ID=', mat2str(box_vm_time_series_big_summary{box_id}{1}(1,1)), ...
                ', Time Length is=', mat2str(time_len), 'h, Lag =', mat2str(max_lag - lag_id + 1), ... 
                ' h');
            title(title_name);
            set(gcf, 'paperpositionmode', 'auto');
            print('-dpng2','-r300', strcat(path, 'PM_', mat2str(box_vm_time_series_big_summary{box_id}{1}(1,1)),...
                   '_Start_', mat2str(corr_different_series{box_id}{2, time_slot}(1)), ...
                   '_Length_', mat2str(time_len), '_Lag_', mat2str(max_lag - lag_id + 1)));
    
        end
    end
    close all
end


