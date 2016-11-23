function [candidate, solution_or_not, left_capacity] = greedy_find_approximate(sample, ...
                                                counts, total_cpu_capacity)
% This is the greedy algorithm that find the approximate solution for the
% problem.

    % Calculate the ticket coefficient, which evaluate the efficiency of
    % decrease the demands for each vm in each step. Bigger the
    % *coeffcient* is, better to choose this as the candidate
    coeff_ticket = {};
    size_box = numel(sample);
    for machine_id = 1 : size_box
        decrease_demands = sample{machine_id}(1:end-1) - sample{machine_id}(2:end);
        increase_tickets = counts{machine_id}(2:end) - counts{machine_id}(1:end-1);
        if numel(increase_tickets) == 0
            coeff_ticket{machine_id} = 0;
        else
            coeff_ticket{machine_id} = decrease_demands' ./ increase_tickets; 
        end
        % disp(coeff_ticket{machine_id});
        coeff_ticket{machine_id} = [0, coeff_ticket{machine_id}];
    end

    % Start from the beginning
    
    candidate = ones(1,size_box);
    sum_capacity_cand = 0; stop_sign_total_number = 0;
    for machine_id = 1 : size_box
        sum_capacity_cand = sum_capacity_cand + sample{machine_id}(1);
        stop_sign_total_number = stop_sign_total_number + numel(counts{machine_id});
    end
   
    while sum_capacity_cand > total_cpu_capacity
        cand_coeff = [];
        for machine_id = 1 : size_box
            if candidate(machine_id) + 1 <= numel(coeff_ticket{machine_id})
                cand_coeff = [cand_coeff, coeff_ticket{machine_id}(candidate(machine_id) + 1)];
            else
                cand_coeff = [cand_coeff, 0];
            end
        end
        [max_coeff, max_idx] = max(cand_coeff);
        candidate(max_idx) = candidate(max_idx) + 1;
    
        sum_capacity_cand = sum_capacity_cand - ...
            (sample{max_idx}(candidate(max_idx)-1) - sample{max_idx}(candidate(max_idx)));
        
        if sum_capacity_cand <= total_cpu_capacity
            break;
        end
        
        % We may not find an optimal case (actually cannot happen), just in
        % case of extreme case that total capacity is zero
        if sum(candidate) == stop_sign_total_number
            break;
        end
    end
    
    if sum_capacity_cand > total_cpu_capacity
        solution_or_not = false;
%         disp('CUO LA');
%         sum_small_demand = 0;
%         for vm_id = 1 : numel(sample)
%             sum_small_demand = sum_small_demand + sample{vm_id}(end);
%         end
%         disp(sum_small_demand - total_cpu_capacity);
    else
        solution_or_not = true;
    end
    
    left_capacity = max(0, total_cpu_capacity - sum_capacity_cand);
end
