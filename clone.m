function [reduced_ticket, optimal_weight, optimal_copy] = ...
                            clone(vm_time_series, vm_time_series_mem, ...
                                  box_cap, box_activity, box_mem, ...
                                  time_thres, ticket_threshold)
                              
    % use the cloning method, to get the least number of tickets

    T = min(time_thres, numel(vm_time_series(:,1)));% number of time points
    N = numel(vm_time_series(1,:));                 % number of VMs
    M = numel(box_cap);                             % number of Boxes
    
    % variable format
    % X = [weights, clone_variable, ticket_variable]
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
    
    f = [zeros(2 * N*M, 1); ones(M*T,1)];
    intcon = N*M+1 : 2*N*M + M*T;
    
    % constraints
    A = []; b = [];
  
    % contraints for the *Capacity*
    for box_idx = 1 : M
        % 1) CPU capacity constraints (*M * T* different inequality)
        A = [A; zeros(T, (box_idx-1) * N), vm_time_series(1:T, :), ...
                zeros(T, (M - box_idx) * N + N*M + M*T)];
        b = [b; ones(T, 1) * box_cap(box_idx) - box_activity(1:T, box_idx)];
        
        % 2) RAM capacity constraints (*M * T* different inequality)
        A = [A; zeros(T, N*M), zeros(T, (box_idx-1) * N), ...
                vm_time_series_mem(1:T, :), zeros(T, (M - box_idx) * N + M*T)];
        b = [b; ones(T, 1) * box_mem(box_idx)];
    end
    
%     % constraints for the weights
%     for weight_id = 1 : N*M
%         A = [A; zeros(1, weight_id - 1),  -1, zeros(N*M - weight_id), ...
%                 zeros(1, N*M + T*M)];
%         b = [b; 0];    
%     end
    
    % constraints for *weights* and *clone_variable*
    Min_A = 0.01;
    for weight_id = 1 : N*M
        % each weight should be no greater than clone_variable
        A = [A; zeros(1, weight_id - 1), 1, zeros(1, N*M - weight_id), ...
                zeros(1, weight_id - 1), -1, zeros(1, N*M - weight_id), ...
                zeros(1, T*M)];
        b = [b; 0];
        
        % force clone_variable to be zero when w_ij is less than Min_A
        A = [A; zeros(1, weight_id - 1), -1, zeros(1, N*M - weight_id), ...
                zeros(1, weight_id - 1), 1, zeros(1, N*M - weight_id), ...
                zeros(1, T*M)];
        b = [b; 1 - Min_A];
    end
    
    % constraint for *weights* and *ticket_variable*
    Max_A = 10^8; Max_B = 10^8;
    for box_idx = 1 : M
        % force to have tickets, when exceeding ticket_threshod
        A = [A; zeros(T, (box_idx-1) * N), vm_time_series(1:T, :), ...
                zeros(T, (M - box_idx) * N), zeros(T, N*M), ...
                zeros(T, (box_idx-1) * T), -ones(T, T) * Max_A, ...
                zeros(T, (M - box_idx) * T)];
        b = [b; ones(T, 1) * box_cap(box_idx) * ticket_threshold - box_activity(1:T, box_idx)];
        
        % force to have no ticket, when not exceeding ticket_threshold
        A = [A; zeros(T, (box_idx-1) * N), -vm_time_series(1:T, :), ...
                zeros(T, (M - box_idx) * N), zeros(T, N*M), ...
                zeros(T, (box_idx-1) * T), ones(T, T) * Max_B, ...
                zeros(T, (M - box_idx) * T)];
        b = [b; ones(T, 1) * (Max_B - box_cap(box_idx) * ticket_threshold) + box_activity(1:T, box_idx)];
    end
    
    % Aeq and beq
    Aeq = []; beq = [];
    for vm_idx = 1 : N
        weight_vm = [];
        for box_idx = 1 : M
            weight_vm = [weight_vm, zeros(1, vm_idx - 1), 1, zeros(1, N - vm_idx)];
        end
        Aeq = [Aeq; weight_vm, zeros(1, N*M + T*M)];
        beq = [beq; 1];
    end
    
    % upper bound and lower bound
    lb = zeros(2*N*M + M*T, 1);
    ub = ones(2*N*M + M*T, 1);

    options = optimoptions('intlinprog','Display','off');
    [x, fval, exitflag] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
     
    reduced_ticket = 0; optimal_weight = []; optimal_copy = [];
    if numel(x) == 0
        disp('No feasible solution is found: No way!!! It is bug!!!');
    else
        reduced_ticket = fval;
        optimal_weight = x(1: N*M);
        optimal_copy = x(1+N*M : 2*N*M);
    end
end
