function [num_ticket, final_vm_position] = min_ticket(num_box, num_vm, box_left_cap,...
                                                num_ticket, vm_position, final_vm_position,...
                                                original_num_vm, vm_time_series,...
                                                box_cap, ticket_threshold)
%recursive function to calculate the minimum number of box tickets
    for box_id = 1 : num_box
        box_left_cap(:, box_id) = box_left_cap(:, box_id) - vm_time_series(:, original_num_vm - num_vm + 1);
        
        % check if the box cannot stand all the current VMs
        box_less_than_zero = sum(box_left_cap(box_id) < 0);
        if box_less_than_zero == 0
      
            vm_position(end+1) = box_id;
            
            % if we have the last vm 
            if num_vm == 1
                ticket_in_this_case = ...
                    sum(sum(box_left_cap ./ box_cap < 1 - ticket_threshold / 100));
                
                % disp(ticket_in_this_case);
                
                if ticket_in_this_case < num_ticket
                    num_ticket = ticket_in_this_case;
                    final_vm_position = vm_position;   
                end
                
            % if we still have some more VMs left
            else
                num_vm = num_vm - 1;
                [num_ticket, final_vm_position] = min_ticket(num_box, ...
                                                    num_vm, box_left_cap,...
                                                    num_ticket, vm_position, final_vm_position, ...
                                                    original_num_vm, vm_time_series,...
                                                    box_cap, ticket_threshold);
                num_vm = num_vm + 1;
            end 
            
            % recover to no placement for this VM
            vm_position = vm_position(1:end-1);
        end
        
        % recover to the original state
        box_left_cap(:, box_id) = box_left_cap(:, box_id) + vm_time_series(:, original_num_vm - num_vm + 1);
    end
    
    if num_vm == original_num_vm
        disp({'the best ticket number is ', mat2str(num_ticket)});
        return
    end
    
end