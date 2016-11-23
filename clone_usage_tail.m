function [reduced_ticket, optimal_weight, optimal_copy, ...
          optimal_ticket, optimal_target_ticket, ...
          num_of_clone, num_of_migrate, final_time_series] = ...
                   clone_usage_tail(vm_time_series, vm_time_series_mem, ...
                                    box_cap, box_activity, ...
                                    box_mem, box_mem_activity, ...
                                    test_time_thres, ...
                                    ticket_threshold, target_usage_tail, ...
                                    original_vm_placement)

    % use the cloning method, to get the least number of tickets

    T = min(test_time_thres, numel(vm_time_series(:,1)));                 
                                                    % number of time points
    N = numel(vm_time_series(1,:));                 % number of VMs
    M = numel(box_cap);                             % number of Boxes
    
    % variable format
    % X = [weights, clone_variable, ticket_variable, target_usage_variable]
    % weights = 
    %           [w_11, w_21, w_31, ..., w_N1,
    %            w_12, w_22, w_32, ..., w_N2,
    %            ...
    %            w_1M, w_2M, w_3M, ..., w_NM]
    % clone_variable = 
    %           [a_11, a_21, a_31, ..., a_N1,
    %            a_12, a_22, a_32, ..., a_N2,
    %            ...
    %            a_1M, a_2M, a_3M, ..., a_NM]
    % ticket_variable = 
    %           [I_11, I_12, I_13, ..., I_1T,
    %            I_21, I_22, I_23, ..., I_2T,
    %            ...
    %            I_M1, I_M2, I_M3, ..., I_MT]
    % target_usage_variable = 
    %           [Y_1, Y_2, Y_3, ..., Y_M]
    
    f = [zeros(2*N*M + M*T, 1); ones(M,1)];
    intcon = N*M+1 : 2*N*M + M*T + M;
    
    % constraints
    A = []; b = [];
  
    % contraints for the *Capacity*
    for box_idx = 1 : M
        % 1) CPU capacity constraints (*M * T* different inequality)
        A = [A; zeros(T, (box_idx-1) * N), vm_time_series(1:T, :), ...
                zeros(T, (M - box_idx) * N + N*M + M*T + M)];
        b = [b; ones(T, 1) * box_cap(box_idx) - box_activity(1:T, box_idx)];
        
        % 2) RAM capacity constraints (*M * T* different inequality)
        A = [A; zeros(T, (box_idx-1) * N), vm_time_series_mem(1:T, :), ...
                zeros(T, (M - box_idx) * N + N*M + M*T + M)];
        b = [b; ones(T, 1) * box_mem(box_idx) - box_mem_activity(1:T, box_idx)];
    end
       
    % constraints for *weights* and *clone_variable*
    Min_A = 0.1;
    for weight_id = 1 : N*M
        % each weight should be no greater than clone_variable
        A = [A; zeros(1, weight_id - 1), 1, zeros(1, N*M - weight_id), ...
                zeros(1, weight_id - 1), -1, zeros(1, N*M - weight_id), ...
                zeros(1, T*M + M)];
        b = [b; 0];
        
        % force clone_variable to be zero when w_ij is less than Min_A
        A = [A; zeros(1, weight_id - 1), -1, zeros(1, N*M - weight_id), ...
                zeros(1, weight_id - 1), 1, zeros(1, N*M - weight_id), ...
                zeros(1, T*M + M)];
        b = [b; 1 - Min_A];
    end
    
    % constraint for *weights* and *ticket_variable* and *target_usage_tail*
    Max_A = 10^8; Max_B = 10^8;
    Max_C = 10^8; Max_D = 10^8;
    all_ticket_cand = diag(ones(T,1));
    for box_idx = 1 : M
        % force to have tickets, when exceeding ticket_threshod
        A = [A; zeros(T, (box_idx-1) * N), vm_time_series(1:T, :), ...
                zeros(T, (M - box_idx) * N), zeros(T, N*M), ...
                zeros(T, (box_idx-1) * T), -all_ticket_cand * Max_A, ...
                zeros(T, (M - box_idx) * T + M)];
        b = [b; box_cap(box_idx) * ticket_threshold / 100 - box_activity(1:T, box_idx)];
        
        % force to have no ticket, when not exceeding ticket_threshold
        A = [A; zeros(T, (box_idx-1) * N), -vm_time_series(1:T, :), ...
                zeros(T, (M - box_idx) * N), zeros(T, N*M), ...
                zeros(T, (box_idx-1) * T), all_ticket_cand * Max_B, ...
                zeros(T, (M - box_idx) * T + M)];
        b = [b; Max_B - box_cap(box_idx) * ticket_threshold / 100 + box_activity(1:T, box_idx)];
    
        % force to have usage tail ticket, when exceeding target_usage_tail
        A = [A; zeros(1, 2*N*M),...
                zeros(1, (box_idx-1) * T), ones(1, T), zeros(1, (M - box_idx) * T), ...
                zeros(1, box_idx-1), -Max_C, zeros(1, M - box_idx)];
        b = [b; round((1-target_usage_tail/100) * T)];
        
        % force to have no usage tail ticket, when not exceeding the target
        % usage tail
        A = [A; zeros(1, 2*N*M),...
                zeros(1, (box_idx-1) * T), -ones(1, T), zeros(1, (M - box_idx) * T), ...
                zeros(1, box_idx-1), Max_D, zeros(1, M - box_idx)];
        b = [b; Max_D - round((1-target_usage_tail/100) * T)]; 
    end
    
    % Aeq and beq
    Aeq = []; beq = [];
    for vm_idx = 1 : N
        weight_vm = [];
        for box_idx = 1 : M
            weight_vm = [weight_vm, zeros(1, vm_idx - 1), 1, zeros(1, N - vm_idx)];
        end
        Aeq = [Aeq; weight_vm, zeros(1, N*M + T*M + M)];
        beq = [beq; 1];
    end
    
    % upper bound and lower bound
    lb = zeros(2*N*M + M*T + M, 1);
    ub = ones(2*N*M + M*T + M, 1);

    options = optimoptions('intlinprog','Display','off');
    [x, fval, exitflag] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
     
    reduced_ticket = 0; optimal_weight = []; optimal_copy = [];
    optimal_ticket = []; optimal_target_ticket = [];
    if numel(x) == 0
        disp('No feasible solution is found: No way!!! It is bug!!!');
    else
        reduced_ticket = fval;
        x(intcon) = round(x(intcon));
        optimal_weight = vec2mat(x(1: N*M), N);
        optimal_copy = vec2mat(x(1+N*M : 2*N*M), N);
        optimal_ticket = vec2mat(x(2*N*M+1 : 2*N*M + M*T), T);
        optimal_target_ticket = x(2*N*M + M* T + 1 : end);
        
        % Get the number of clones and migration
        num_of_clone = 0; num_of_migrate = 0;
        for j = 1 : M
            for i = 1 : N
                if optimal_weight(j, i) > 0 && original_vm_placement(j, i) == 0
                    if abs(optimal_weight(j, i) -  1) <= 0.01
                        num_of_migrate = num_of_migrate + 1;
                    else
                        num_of_clone = num_of_clone + 1;
                    end
                end
            end
        end
        
        % Get the final placement time series
        final_time_series = {};
        for j = 1 : M
            final_time_series{j}{1} = [box_activity(:, j), box_mem_activity(:, j)];
            for i = 1 : N
                if optimal_weight(j, i) > 0 
                    final_time_series{j}{end+1} = optimal_weight(j, i) * ...
                            [vm_time_series(:,i), vm_time_series_mem(:,i)];
                    final_time_series{j}{1} = final_time_series{j}{1} + ...
                                              final_time_series{j}{end};
                end
            end
            final_time_series{j}{1}(:,3) = final_time_series{j}{1}(:,1) / box_cap(j) * 100;
            final_time_series{j}{1}(:,4) = final_time_series{j}{1}(:,2) / box_mem(j) * 100;
        end

    end
end
