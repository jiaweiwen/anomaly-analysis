% calculate the conditional probability between box tickets and vm tickets

close all; clear; clc

load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat

mkdir('../New_Data_box_vm_tickets_characterization_new');
path = '../New_Data_box_vm_tickets_characterization_new/';

ticket_thres = [60, 70, 80];
test_lag_cand = [0];

grat = 0.1;

time_grat = 900;

for lag_id = 1 : numel(test_lag_cand)
    test_lag = test_lag_cand(lag_id);
    
    BOX_TICKET_VM_TICKET_CPU = {[],[],[]}; BOX_TICKET_VM_TICKET_MEM = {[],[],[]};
    
    BOX_TICKET_INTER_ARRIVAL_CPU = {[], [], []};
    BOX_TICKET_INTER_ARRIVAL_MEM = {[], [], []};
    
    box_num = 1;
    size_box_vm = size(box_vm_time_series_summary);
    
    for box_id = 1 : size_box_vm(2)

        size_box = numel(box_vm_time_series_summary{1, box_id});

        % If we don't have time series
        if size_box < 2
            continue;
        end
        
%         if box_num >= 1000
%             break;
%         end

        if numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1)) <= 10
            continue;
        end

        % First extract the number of tickets for CPU and RAM for different
        % thresholds
        box_cpu = box_vm_time_series_summary{1, box_id}{1,1}(:,4);
        box_mem = box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 4);
        time_stamp = box_vm_time_series_summary{1, box_id}{1,1}(:, 1);
        box_ticket_overtime_cpu = [];
        box_ticket_overtime_mem = [];
        for ticket_id = 1 : numel(ticket_thres)
            box_ticket_overtime_cpu(:, ticket_id) = box_cpu > ticket_thres(ticket_id);
            box_ticket_overtime_mem(:, ticket_id) = box_mem > ticket_thres(ticket_id);
        end  

        vm_ticket_overtime_cpu = {}; 
        vm_ticket_overtime_mem = {}; 
        for vm_id = 2 : size_box
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%consider VM%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vm_cpu = box_vm_time_series_summary{1, box_id}{1,vm_id}(:,5);
            vm_mem = box_vm_time_series_summary_mem{1, box_id}{1,vm_id}(:, 5);
            for ticket_id = 1 : numel(ticket_thres)
                % cpu
                temp_cpu = vm_cpu > ticket_thres(ticket_id);
                vm_cpu_lag = [];
                for lag = -test_lag : 0
                    vm_cpu_lag(:, end+1) = [temp_cpu(1 - lag : end); zeros(-lag,1)];
                end
                for lag = 1 : test_lag
                    vm_cpu_lag(:, end+1) = [zeros(lag, 1); temp_cpu(1 : end - lag)];
                end         
                if test_lag ~= 0
                    vm_ticket_overtime_cpu{ticket_id}(:, vm_id - 1) = max(vm_cpu_lag');
                else
                    vm_ticket_overtime_cpu{ticket_id}(:, vm_id - 1) = vm_cpu_lag;
                end
                
                % mem
                temp_mem = vm_mem > ticket_thres(ticket_id);
                vm_mem_lag = [];
                for lag = -test_lag : 0
                    vm_mem_lag(:, end+1) = [temp_mem(1 - lag : end); zeros(-lag,1)];
                end
                for lag = 1 : test_lag
                    vm_mem_lag(:, end+1) = [zeros(lag, 1); temp_mem(1 : end - lag)];
                end             
                if test_lag ~= 0
                    vm_ticket_overtime_mem{ticket_id}(:, vm_id - 1) = max(vm_mem_lag');
                else
                    vm_ticket_overtime_mem{ticket_id}(:, vm_id - 1) = vm_mem_lag;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        % check if box has ticket, how many vms have tickets
        for ticket_id = 1 : numel(ticket_thres)
            % CPU
            if numel(vm_ticket_overtime_cpu{ticket_id}(1,:)) == 1
                vm_ticket_overtime_cpu_lag_final = vm_ticket_overtime_cpu{ticket_id}';
            else
                vm_ticket_overtime_cpu_lag_final = sum(vm_ticket_overtime_cpu{ticket_id}');
            end
            BOX_TICKET_VM_TICKET_CPU{ticket_id} = [BOX_TICKET_VM_TICKET_CPU{ticket_id}; ...
                                                   box_ticket_overtime_cpu(:, ticket_id), ...
                                                   vm_ticket_overtime_cpu_lag_final', ...
                                                   ones(numel(vm_ticket_overtime_cpu_lag_final), 1) * (size_box-1)];
                                                
            non_zero_idx = find(box_ticket_overtime_cpu(:, ticket_id) == 1);         
            if numel(non_zero_idx) > 1
                chosen_time = time_stamp(non_zero_idx);
                BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id} = [BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id};  
                                                           chosen_time(2:end) - chosen_time(1:end-1) - time_grat,...
                                                           ones(numel(non_zero_idx)-1,1) * box_id];
            end

            % MEM
            if numel(vm_ticket_overtime_mem{ticket_id}(1,:)) == 1
                vm_ticket_overtime_mem_lag_final = vm_ticket_overtime_mem{ticket_id}';
            else
                vm_ticket_overtime_mem_lag_final = sum(vm_ticket_overtime_mem{ticket_id}');
            end           
            BOX_TICKET_VM_TICKET_MEM{ticket_id} = [BOX_TICKET_VM_TICKET_MEM{ticket_id}; ...
                                                   box_ticket_overtime_mem(:, ticket_id), ...
                                                   vm_ticket_overtime_mem_lag_final', ...
                                                   ones(numel(vm_ticket_overtime_mem_lag_final), 1) * (size_box-1)];
                                               
            non_zero_idx = find(box_ticket_overtime_mem(:, ticket_id) == 1);         
            if numel(non_zero_idx) > 1
                chosen_time = time_stamp(non_zero_idx);
                BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id} = [BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id};  
                                                           chosen_time(2:end) - chosen_time(1:end-1) - time_grat,...
                                                           ones(numel(non_zero_idx)-1,1) * box_id];
            end
        end
        
        box_num = box_num + 1;
    end
    
    % Calculate the overall probability to have BOX Tickets
    % Calculate the mean time to have a Box tickets
    prob_box_ticket = []; mean_time_have_box_ticket = [];
    for ticket_id = 1 : numel(ticket_thres)
        box_cpu = BOX_TICKET_VM_TICKET_CPU{ticket_id}(:,1);
        total_time = numel(box_cpu);
        ticket_time = sum(box_cpu == 1);
        prob_box_ticket(ticket_id,1) = ticket_time / total_time;
        mean_time_have_box_ticket(ticket_id,1) = (total_time - ticket_time) * time_grat / (60 * ticket_time);
        
        box_mem = BOX_TICKET_VM_TICKET_MEM{ticket_id}(:,1);
        total_time = numel(box_mem);
        ticket_time = sum(box_mem == 1);
        prob_box_ticket(ticket_id,2) = ticket_time / total_time;
        mean_time_have_box_ticket(ticket_id,2) = (total_time - ticket_time) * time_grat / (60 * ticket_time);    
    end
    
    % Plot the figures of conditional probability
    for ticket_id = 1 : numel(ticket_thres)
        %%%%%%%%%%%%%%%%%%%%%%%%%% cpu first %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check the distribution of the inter-arrival time
        upper_arrival = prctile(BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id}(:,1), 95);
        interest_idx = BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id}(:,1) < upper_arrival;
        max_inter = max(BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id}(interest_idx,1));
        min_inter = min(BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id}(interest_idx,1));
        check_inter = min_inter : time_grat : max_inter;
        [N, edges] = histcounts(BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id}(interest_idx,1), check_inter);
        % pd = fitdist(x,'Exponential')
        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        h = bar(N/sum(N));
        set(gca, 'xticklabel', edges(1:end-1)/60);
        set(gca, 'xlim', [0 numel(check_inter)]);
        xlabel('Inter-arrival of Box Tickets (min)'); ylabel('PDF');
        % set(gca, 'xscale', 'log');
        title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'pdf_inter_arrival_box_tickets_',mat2str(ticket_thres(ticket_id)),'_cpu_lag_', mat2str(test_lag)));
        
        % calculate *Prob(Any VM has Tickets | BOX has Tickets)
        box_has_ticket_idx = BOX_TICKET_VM_TICKET_CPU{ticket_id}(:, 1) == 1;
        vm_correspond = BOX_TICKET_VM_TICKET_CPU{ticket_id}(box_has_ticket_idx, 2);
        prob_vm_any_box = sum(vm_correspond >= 1) / sum(box_has_ticket_idx);
        disp(strcat('CPU Ticket thres = ', mat2str(ticket_thres(ticket_id))));
        disp(strcat('CPU Prob(Any VM has Tickets | BOX has Tickets) = ', mat2str(prob_vm_any_box)));
        
        % calculate *Prob(BOX has Tickets | Any VM has Tickets)
        vm_has_ticket_idx = BOX_TICKET_VM_TICKET_CPU{ticket_id}(:, 2) >= 1;
        box_correspond = BOX_TICKET_VM_TICKET_CPU{ticket_id}(vm_has_ticket_idx, 1);
        prob_box_vm_any = sum(box_correspond > 0) / sum(vm_has_ticket_idx);
        disp(strcat('CPU Ticket thres = ', mat2str(ticket_thres(ticket_id))));
        disp(strcat('CPU Prob(BOX has Tickets | Any VM has Tickets) = ', mat2str(prob_box_vm_any)));
        
        % vm_ratio_with_ticket = ceil(BOX_TICKET_VM_TICKET_CPU{ticket_id}(:,2)/grat) * grat;
        vm_ratio_with_ticket = BOX_TICKET_VM_TICKET_CPU{ticket_id}(:,2);
        
        % unique_vm_ratio = [unique(vm_ratio_with_ticket); max(vm_ratio_with_ticket) + grat];
        unique_vm_ratio = [0, 1, 2, 3, 4, 5, 6, 8, 16, 32, 64, 128, 256];
        [N, edges] = histcounts(vm_ratio_with_ticket, unique_vm_ratio);
        
        fig = figure;
        set(fig, 'Position', [200, 200, 600, 400]);
        h = bar(N);
        set(gca, 'xticklabel', unique_vm_ratio(1:end-1));
        set(gca, 'xlim', [0 numel(unique_vm_ratio)]);
        xlabel('Number of VM w/ Tickets'); ylabel('Histogram');
        set(gca, 'ylim', [1 1000000]); set(gca, 'yscale' ,'log');
        title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'hist_vm_with_tickets_',mat2str(ticket_thres(ticket_id)),'_cpu_lag_', mat2str(test_lag)));
        
        % calculate the probability
        prob_box_with_vm = []; prob_vm_with_box = [];
        for vm_ratio_id = 1 : numel(unique_vm_ratio)-1
            vm_idx = find(vm_ratio_with_ticket == unique_vm_ratio(vm_ratio_id));
            all_box_number = numel(vm_idx);
            box_with_ticket_number = sum(BOX_TICKET_VM_TICKET_CPU{ticket_id}(vm_idx,1));
            prob_box_with_vm(vm_ratio_id) = box_with_ticket_number / all_box_number;
             
            box_idx = find(BOX_TICKET_VM_TICKET_CPU{ticket_id}(:,1) == 1);
            all_box_number = numel(box_idx);
            vm_with_ticket_number = sum(vm_ratio_with_ticket(box_idx) == unique_vm_ratio(vm_ratio_id));
            prob_vm_with_box(vm_ratio_id) = vm_with_ticket_number / all_box_number;
        end
        disp(strcat('CPU Ticket thres = ', mat2str(ticket_thres(ticket_id))));   
        disp(strcat('CPU Prob(BOX has Tickets | different number VM has Tickets) = ', mat2str(prob_box_with_vm)));
        
        disp(strcat('CPU Ticket thres = ', mat2str(ticket_thres(ticket_id))));   
        disp(strcat('CPU Prob(different number VM has Tickets | BOX has Tickets) = ', mat2str(prob_vm_with_box)));
        
        prob_cpu = [prob_box_with_vm; prob_vm_with_box];
        fig = figure;
        set(fig, 'Position', [200, 200, 600, 400]);
        h = bar(prob_cpu', 'grouped');   
        set(h(1), 'facecolor', 'r');
        set(h(2), 'facecolor', 'g');
        xlabel('Ratio of VM w/ Tickets'); ylabel('Prob. w/ Tickets');
        set(gca, 'ytick', [0: 0.1 : 1]); set(gca, 'ylim', [0 1]); 
        set(gca, 'xticklabel', unique_vm_ratio(1:end-1));
        set(gca, 'xlim', [0 numel(unique_vm_ratio)]);
        title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        h = legend('Prob(Box | VM)', ...
                   'Prob(VM | Box)');
        set(h, 'box', 'on', 'location', 'northwest');
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'prob_vm_box_with_tickets_',mat2str(ticket_thres(ticket_id)),'_cpu_lag_', mat2str(test_lag)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% mem second %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check the distribution of the inter-arrival time
        upper_arrival = prctile(BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id}(:,1), 95);
        interest_idx = BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id}(:,1) <= upper_arrival;
        max_inter = max(BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id}(interest_idx,1));
        min_inter = min(BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id}(interest_idx,1));
        if min_inter == max_inter
            max_inter = min_inter + time_grat;
        end
        check_inter = min_inter : time_grat : max_inter;
        [N, edges] = histcounts(BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id}(interest_idx,1), check_inter);
        % pd = fitdist(x,'Exponential')
        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        h = bar(N'/sum(N));
        set(gca, 'xticklabel', edges(1:end)/60);
        set(gca, 'xlim', [0 numel(check_inter)]);
        xlabel('Inter-arrival of Box Tickets (min)'); ylabel('PDF');
        % set(gca, 'xscale', 'log');
        title(strcat('RAM, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'pdf_inter_arrival_box_tickets_',mat2str(ticket_thres(ticket_id)),'_mem_lag_', mat2str(test_lag)));
        
        % calculate *Prob(Any VM has Tickets | BOX has Tickets)
        box_has_ticket_idx = BOX_TICKET_VM_TICKET_MEM{ticket_id}(:, 1) == 1;
        vm_correspond = BOX_TICKET_VM_TICKET_MEM{ticket_id}(box_has_ticket_idx, 2);
        prob_vm_any_box = sum(vm_correspond >= 1) / sum(box_has_ticket_idx);
        disp(strcat('MEM Ticket thres = ', mat2str(ticket_thres(ticket_id))));
        disp(strcat('MEM Prob(Any VM has Tickets | BOX has Tickets) = ', mat2str(prob_vm_any_box)));
        
        % calculate *Prob(BOX has Tickets | Any VM has Tickets)
        vm_has_ticket_idx = BOX_TICKET_VM_TICKET_MEM{ticket_id}(:, 2) >= 1;
        box_correspond = BOX_TICKET_VM_TICKET_MEM{ticket_id}(vm_has_ticket_idx, 1);
        prob_box_vm_any = sum(box_correspond > 0) / sum(vm_has_ticket_idx);
        disp(strcat('MEM Ticket thres = ', mat2str(ticket_thres(ticket_id))));
        disp(strcat('MEM Prob(BOX has Tickets | Any VM has Tickets) = ', mat2str(prob_box_vm_any)));
        
        % vm_ratio_with_ticket = ceil(BOX_TICKET_VM_TICKET_MEM{ticket_id}(:,2)/grat) * grat;
        vm_ratio_with_ticket = BOX_TICKET_VM_TICKET_MEM{ticket_id}(:,2);
        
        % unique_vm_ratio = [unique(vm_ratio_with_ticket); max(vm_ratio_with_ticket) + grat];
        unique_vm_ratio = [0, 1, 2, 3, 4, 5, 6, 8, 16, 32, 64, 128, 256];
        [N, edges] = histcounts(vm_ratio_with_ticket, unique_vm_ratio);
        
        fig = figure;
        set(fig, 'Position', [200, 200, 600, 400]);
        h = bar(N);
        set(gca, 'xticklabel', unique_vm_ratio(1:end-1));
        set(gca, 'xlim', [0 numel(unique_vm_ratio)]);
        xlabel('Number of VM w/ Tickets'); ylabel('Histogram');
        set(gca, 'ylim', [1 1000000]); set(gca, 'yscale' ,'log');
        title(strcat('RAM, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'hist_vm_with_tickets_',mat2str(ticket_thres(ticket_id)),'_mem_lag_', mat2str(test_lag)));
        
        % calculate the probability
        prob_box_with_vm = []; prob_vm_with_box = [];
        for vm_ratio_id = 1 : numel(unique_vm_ratio)-1
            vm_idx = find(vm_ratio_with_ticket == unique_vm_ratio(vm_ratio_id));
            all_box_number = numel(vm_idx);
            box_with_ticket_number = sum(BOX_TICKET_VM_TICKET_MEM{ticket_id}(vm_idx,1));
            prob_box_with_vm(vm_ratio_id) = box_with_ticket_number / all_box_number;
             
            box_idx = find(BOX_TICKET_VM_TICKET_MEM{ticket_id}(:,1) == 1);
            all_box_number = numel(box_idx);
            vm_with_ticket_number = sum(vm_ratio_with_ticket(box_idx) == unique_vm_ratio(vm_ratio_id));
            prob_vm_with_box(vm_ratio_id) = vm_with_ticket_number / all_box_number;
        end
        
        disp(strcat('MEM Ticket thres = ', mat2str(ticket_thres(ticket_id))));   
        disp(strcat('MEM Prob(BOX has Tickets | different number VM has Tickets) = ', mat2str(prob_box_with_vm)));
        
        disp(strcat('MEM Ticket thres = ', mat2str(ticket_thres(ticket_id))));   
        disp(strcat('MEM Prob(different number VM has Tickets | BOX has Tickets) = ', mat2str(prob_vm_with_box)));       
        
        prob_mem = [prob_box_with_vm; prob_vm_with_box];
        fig = figure;
        set(fig, 'Position', [200, 200, 600, 400]);
        h = bar(prob_mem', 'grouped');
        set(h(1), 'facecolor', 'r');
        set(h(2), 'facecolor', 'g');
        xlabel('Ratio of VM w/ Tickets'); ylabel('Prob. w/ Tickets');
        set(gca, 'ytick', [0: 0.1 : 1]); set(gca, 'ylim', [0 1]);
        set(gca, 'xticklabel', unique_vm_ratio(1:end-1));
        set(gca, 'xlim', [0 numel(unique_vm_ratio)]);
        title(strcat('RAM, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        h = legend('Prob(Box | VM)', ...
                   'Prob(VM | Box)');
        set(h, 'box', 'on', 'location', 'north');
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'prob_vm_box_with_tickets_',mat2str(ticket_thres(ticket_id)),'_mem_lag_', mat2str(test_lag)));
        
        close all
    end

    save(strcat(path, 'BOX_TICKET_VM_TICKET_CPU', strcat(test_lag)), 'BOX_TICKET_VM_TICKET_CPU');
    save(strcat(path, 'BOX_TICKET_VM_TICKET_MEM', strcat(test_lag)), 'BOX_TICKET_VM_TICKET_MEM');
    save(strcat(path, 'BOX_TICKET_INTER_ARRIVAL_CPU', strcat(test_lag)), 'BOX_TICKET_INTER_ARRIVAL_CPU');
    save(strcat(path, 'BOX_TICKET_INTER_ARRIVAL_MEM', strcat(test_lag)), 'BOX_TICKET_INTER_ARRIVAL_MEM');
    
end