function [original_ticket, reduced_ticket, original_vm_position, final_vm_position,...
          optimal_ticket, optimal_weight, optimal_copy] = ...
                            migration(box_vm_time_series_summary_mem,...
                                      box_vm_time_series_summary, test_box,...
                                      ticket_threshold)
% use the migration/placement, to get the least number of tickets

    path = '../New_Data_Migration_Motivation_Fig/';
    mkdir(path);
    
    box_id = test_box(:,1); box_cap = test_box(:, 4);
    box_number_vm = round(test_box(:,end-1) ./ test_box(:, end));
    vm_time_series = []; box_activity = []; 
    vm_time_series_mem = []; box_mem = [];
    for box_cand = 1 : numel(box_id) 
        box_idx = box_id(box_cand); 
        
        box_mem(end+1) = nanmean(box_vm_time_series_summary_mem{box_idx}{1}(:,4));
        
        box_demands = box_vm_time_series_summary{box_idx}{1}(:,3);
        
        vm_demands = [];
        size_box = numel(box_vm_time_series_summary{box_idx});
        for vm_id = 2 : size_box
            vm_demands(:, end+1) = box_vm_time_series_summary{box_idx}{vm_id}(:,4);
            vm_time_series_mem(:, end + 1) = box_vm_time_series_summary_mem{box_idx}{vm_id}(:,5) ...
                                            .* box_vm_time_series_summary_mem{box_idx}{vm_id}(:,4) / 100;
        end
        vm_time_series(:, end + 1 : end + size_box - 1) = vm_demands;
        
        if size_box == 2
            box_activity(:, end+1) = box_demands - vm_demands;
        else
            box_activity(:, end+1) = box_demands - sum(vm_demands')';
        end
        
        cc = hsv(size_box);
        fig = figure;
        set(fig, 'Position', [200 200 300 200]);
        plot(box_demands, 'color', cc(1,:));
        hold on
        for vm_id = 2 : size_box
            plot(vm_demands(:, vm_id-1), 'color', cc(vm_id,:));
            hold on
        end
        thres = box_cap(box_cand)*ticket_threshold/100;
        plot(get(gca,'xlim'), [thres thres], 'k')
        xlabel('Time (15 min)'); ylabel('CPU Demand (GHZ)');
        set(gca, 'ylim', [0 box_cap(box_cand)]);
        set(gca, 'ytick', [0 : box_cap(box_cand) / 5 : box_cap(box_cand)]);
        set(gca, 'fontsize', 10);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'orginal_overtime_plots_', mat2str(box_cand)));
    
    end
    
    % first check if we have enough box capacity
    if numel(vm_time_series(1,:)) == 1
        total_vm_demands = vm_time_series;
    else
        total_vm_demands = sum(vm_time_series')';
    end
    
    if numel(box_activity(1,:)) == 1
        total_box_demands = box_activity;
    else
        total_box_demands = sum(box_activity')';
    end
    
    total_demands = total_vm_demands + total_box_demands;
    
    fig = figure;
    set(fig, 'Position', [200 200 300 200]);
    plot(sum(box_cap) * ticket_threshold /100 - total_demands)
    xlabel('Time (15 min)'); ylabel('Enough Capacity (+) or Not (-)');
    set(gca, 'fontsize', 10);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'total_available_75_compare_used'));
    
    % just assume that the box activity won't change after migration
    box_capacity_all = [];
    for box_idx = 1 : numel(box_id)
        box_capacity_all(:, end+1) = box_cap(box_idx) - box_activity(:, box_idx);
    end
    
    % do the optimal migration
    % first is the original case
    number_of_ticket = sum(test_box(:,2));
    original_vm_position = [];
    for box_idx = 1 : numel(box_id)
        original_vm_position(end+1 : end + box_number_vm(box_idx)) = ones(1, box_number_vm(box_idx)) * box_idx;
    end
    
    % second migration
    [num_ticket, final_vm_position] = min_ticket(numel(box_id), numel(vm_time_series(1,:)), box_capacity_all,...
                                                number_of_ticket, [], [], ...
                                                numel(vm_time_series(1,:)), vm_time_series,...
                                                box_capacity_all, ticket_threshold);
                                                                                       
    if numel(final_vm_position) == 0
        disp('The original case is the best case');
    else
        disp({'The original tickets is ', mat2str(number_of_ticket)});
        disp({'The migration tickets is ', mat2str(num_ticket)});
        disp(original_vm_position);
        disp(final_vm_position);
    end
    
    original_ticket = number_of_ticket;
    reduced_ticket = num_ticket;
    
    unique_box = unique(final_vm_position);
    str_name = {'VM1', 'VM2', 'VM3', 'VM4', 'VM5'};
    a11 = 0.5; a12 = 1 - a11; idx_box1 = 3;
    a21 = 0.5; a22 = 1 - a21; idx_box2 = 2;
    for box_idx = 1 : numel(unique_box)
        idx = find(final_vm_position == unique_box(box_idx));
        size_box = numel(idx); 
        cc = hsv(size_box + 1);
        fig = figure;
        set(fig, 'Position', [200 200 300 200]);
        
        for vm_id = 1 : size_box
            if vm_id == 1                             
                total_vm_demands = vm_time_series(:, idx(vm_id));
            else
                total_vm_demands = total_vm_demands + vm_time_series(:, idx(vm_id));
            end
            
            if vm_id == 1 && box_idx == 1
                plot(a11 * vm_time_series(:, idx(vm_id)), 'color', cc(vm_id+1,:), 'linestyle', '--');
            elseif vm_id == 3 && box_idx == 2
                plot(a22 * vm_time_series(:, idx(vm_id)), 'color', cc(vm_id+1,:), 'linestyle', '--');
            else
                plot(vm_time_series(:, idx(vm_id)), 'color', cc(vm_id+1,:));
            end
            hold on
        end
        
        if box_idx == 1
            total_vm_demands = total_vm_demands - a12 * vm_time_series(:, idx(idx_box1));
            clone_vm_cand = find(final_vm_position == unique_box(box_idx+1));
            clone_vm_idx = clone_vm_cand(idx_box2);
            total_vm_demands = total_vm_demands + a21 * vm_time_series(:, clone_vm_idx);
            
            plot(a21 * vm_time_series(:, clone_vm_idx), 'color',[0.9, 0.5, 0], 'linestyle', '-.');
            hold on
        end
        
        if box_idx == 2    
            total_vm_demands = total_vm_demands - a21 * vm_time_series(:, idx(idx_box2));
            clone_vm_cand = find(final_vm_position == unique_box(box_idx-1));
            clone_vm_idx = clone_vm_cand(idx_box1);
            total_vm_demands = total_vm_demands + a12 * vm_time_series(:, clone_vm_idx);
            
            plot(a12 * vm_time_series(:, clone_vm_idx), 'color',[0.9, 0.5, 0], 'linestyle', '-.');
            hold on
        end
        
        total_demands = total_vm_demands + box_activity(:,unique_box(box_idx));
        thres = box_cap(unique_box(box_idx))*ticket_threshold/100;
        plot(total_demands, 'color', cc(1,:));
        hold on
        plot(get(gca,'xlim'), [thres thres], 'k')
        xlabel('Time (15 min)'); ylabel('CPU Demand (GHZ)');
        % h = legend(str_name);
        set(gca, 'ylim', [0 box_cap(unique_box(box_idx))]);
        set(gca, 'ytick', [0 : box_cap(unique_box(box_idx)) / 5 : box_cap(unique_box(box_idx))]);
        set(gca, 'fontsize', 10);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'new_migration_overtime_plots_', mat2str(unique_box(box_idx))));
    end
    
    % third: best the solution from cloning
    time_thres = 2;
    [optimal_ticket, optimal_weight, optimal_copy] = clone(vm_time_series, ...
                                  vm_time_series_mem, box_cap, box_activity, ...
                                  box_mem, time_thres, ticket_threshold);
                              
    
    
end

