function [all_combination, sample, counts] = find_optimal(machine_id, capacity_left, ...
                                        size_box, sample, counts, ...
                                        all_combination, index)
% This function is to find all the possible combination of capacity
% planning
    if machine_id == size_box 
        for poss_id = 1 : numel(sample{machine_id})
            if sample{machine_id}(poss_id) <= capacity_left
                index = [index, poss_id];
                num_ticket = 0;
                for box_vm_id = 1 : machine_id
                    num_ticket = num_ticket + counts{box_vm_id}(index(box_vm_id));
                end
                if numel(all_combination) == 0
                    all_combination = {all_combination{:}, index};
                else
                    current_min_ticket = 0;
                    for box_vm_id = 1 : machine_id
                        current_min_ticket = current_min_ticket + counts{box_vm_id}(all_combination{end}(box_vm_id));
                    end
                    % The smaller one shouldn't be considered
                    if num_ticket == current_min_ticket
                        all_combination = {all_combination{:}, index};
                    elseif num_ticket < current_min_ticket
                        all_combination =  {index};
                    end
                end
                % Pop the last element
                index = index(1:end-1);
                break;
            end
        end
    else  
        for poss_id = 1 : numel(sample{machine_id})
            % Due to the update of sample size, we need to do double
            % checking of the *poss_id*
            if poss_id > numel(sample{machine_id})
                break;
            end
            
            if sample{machine_id}(poss_id) <= capacity_left
                index = [index, poss_id];
                capacity_left = capacity_left - sample{machine_id}(poss_id);
                
%                 disp(strcat('left capacity ', mat2str(capacity_left)));
%                 disp(strcat('machine id ', mat2str(machine_id)));
%                 disp(strcat('poss index is ', mat2str(poss_id)));

                % Update the current minimum size
                if machine_id == size_box - 1 && numel(all_combination) ~= 0
                    current_min_ticket = 0;
                    for box_vm_id = 1 : size_box
                        current_min_ticket = current_min_ticket + counts{box_vm_id}(all_combination{end}(box_vm_id));
                    end
                elseif machine_id == size_box - 1 && numel(all_combination) == 0
                    current_min_ticket = 0;
                    for box_vm_id = 1 : size_box
                        current_min_ticket = current_min_ticket + counts{box_vm_id}(end);
                    end
                end
               
%                 if machine_id == 1
%                     disp('This time work or not?');
%                     disp(counts{1}(end));
%                 end
                
                [all_combination, sample, counts] = find_optimal(machine_id + 1, capacity_left, ...
                                               size_box, sample, counts, ...
                                               all_combination, index);
                % Update the checking length
                % If we have new minimum tickets, we reduce the sample size
                % for each machidisp(counts{1}(end));ne
                if machine_id == size_box - 1 && numel(all_combination) ~= 0
                    new_min_ticket = 0;
                    for box_vm_id = 1 : size_box
                        new_min_ticket = new_min_ticket + counts{box_vm_id}(all_combination{end}(box_vm_id));
                    end
                    % Determine if we need to update the sample size
                    if new_min_ticket < current_min_ticket
                        for box_vm_id = 1 : size_box
                            cand_count_size = find(counts{box_vm_id} >= new_min_ticket, 1, 'first'); 
%                             disp(cand_count_size);
                            if numel(cand_count_size) ~= 0
                                counts{box_vm_id} = counts{box_vm_id}(1:cand_count_size);
                                sample{box_vm_id} = sample{box_vm_id}(1:cand_count_size);
                            end                            
                        end 
%                         disp(strcat(mat2str(machine_id), mat2str(poss_id)));
%                         disp(counts{1}(end));
                    end
                end
                index = index(1:end-1);
                capacity_left = capacity_left + sample{machine_id}(poss_id);
            end
        end
    end
end

