cclose all; clear; clc

path = '../New_Data_box_vm_tickets_characterization_new/';

ticket_thres = [60, 70, 80];
test_lag_cand = [0];

time_grat = 900/3600;

for lag_id = 1 : numel(test_lag_cand)
    
    test_lag = test_lag_cand(lag_id);
    
    load (strcat(path, 'BOX_TICKET_INTER_ARRIVAL_CPU', strcat(test_lag)))
    load (strcat(path, 'BOX_TICKET_INTER_ARRIVAL_MEM', strcat(test_lag)))
    
    for ticket_id = 1 : numel(ticket_thres)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% cpu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box_original_ticket = BOX_TICKET_INTER_ARRIVAL_CPU{ticket_id};
        box_bursty_ticket = [];
        
        col_num = 1; ticket_num = 1; box_id = box_original_ticket(1,2);
        
        total_col = numel(box_original_ticket(:,1));
        while col_num <= total_col - 1
            if box_original_ticket(col_num + 1, 2) ~= box_id
                if box_original_ticket(col_num, 1) == 0
                    box_bursty_ticket(end+1, 1:2) = [ticket_num + 1, box_id];
                else
                    box_bursty_ticket(end+1:end+2, 1:2) = [ticket_num, box_id; ticket_num, box_id];
                end                
                col_num = col_num + 1;
                ticket_num = 1;
                box_id = box_original_ticket(col_num,2);
            else 
                if box_original_ticket(col_num,1) == 0 && box_original_ticket(col_num + 1,1) == 0
                    ticket_num = ticket_num + 1;
                    col_num = col_num + 1;
                    continue;
                end
                
                if box_original_ticket(col_num, 1) == 0
                    box_bursty_ticket(end+1, 1:2) = [ticket_num + 1, box_id];
                else
                    box_bursty_ticket(end+1:end+2, 1:2) = [ticket_num, box_id; ticket_num, box_id];
                end
                col_num = col_num + 1;
                ticket_num = 1;                
            end
        end
        
        % Consider the last ticket
        if box_original_ticket(col_num, 1) == 0
            box_bursty_ticket(end+1, 1:2) = [ticket_num + 1, box_id];
        else
            box_bursty_ticket(end+1:end+2, 1:2) = [ticket_num, box_id; ticket_num, box_id];
        end 
        
        % Plot the distribution of duration of bursty tickets
        max_inter = max(box_bursty_ticket(:,1));
        min_inter = min(box_bursty_ticket(:,1));
        if min_inter == max_inter
            max_inter = min_inter + 1;
        end
        check_inter = min_inter : max_inter+1;
        [N, edges] = histcounts(box_bursty_ticket(:,1), check_inter);
        
        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        plot(edges(1:end-1)*time_grat, N/sum(N), 'ro-', 'linewidth', 1.5)
        set(gca, 'xlim', [0 edges(end-1)*time_grat]);
        set(gca, 'xtick', [0 : 3 : edges(end-1)*time_grat]);
        % set(gca, 'xtick', [0: 3 : edges(end-1)*time_grat/60]);
        xlabel('Duration of Bursty Box Tickets (hour)'); ylabel('PDF');
        %set(gca, 'xscale', 'log');
        title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'pdf_duration_bursty_box_tickets_',mat2str(ticket_thres(ticket_id)),'_cpu_lag_', mat2str(test_lag)));

        %%%%%%%%%%%%%%%%%%%%%%%% mem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box_original_ticket = BOX_TICKET_INTER_ARRIVAL_MEM{ticket_id};
        box_bursty_ticket = [];
        
        col_num = 1; ticket_num = 1; box_id = box_original_ticket(1,2);
        
        total_col = numel(box_original_ticket(:,1));
        while col_num <= total_col - 1
            if box_original_ticket(col_num + 1, 2) ~= box_id
                if box_original_ticket(col_num, 1) == 0
                    box_bursty_ticket(end+1, 1:2) = [ticket_num + 1, box_id];
                else
                    box_bursty_ticket(end+1:end+2, 1:2) = [ticket_num, box_id; ticket_num, box_id];
                end                
                col_num = col_num + 1;
                ticket_num = 1;
                box_id = box_original_ticket(col_num,2);
            else 
                if box_original_ticket(col_num,1) == 0 && box_original_ticket(col_num + 1,1) == 0
                    ticket_num = ticket_num + 1;
                    col_num = col_num + 1;
                    continue;
                end
                
                if box_original_ticket(col_num, 1) == 0
                    box_bursty_ticket(end+1, 1:2) = [ticket_num + 1, box_id];
                else
                    box_bursty_ticket(end+1:end+2, 1:2) = [ticket_num, box_id; ticket_num, box_id];
                end
                col_num = col_num + 1;
                ticket_num = 1;                
            end
        end
        
        % Consider the last ticket
        if box_original_ticket(col_num, 1) == 0
            box_bursty_ticket(end+1, 1:2) = [ticket_num + 1, box_id];
        else
            box_bursty_ticket(end+1:end+2, 1:2) = [ticket_num, box_id; ticket_num, box_id];
        end 
        
        % Plot the distribution of duration of bursty tickets
        max_inter = max(box_bursty_ticket(:,1));
        min_inter = min(box_bursty_ticket(:,1));
        if min_inter == max_inter
            max_inter = min_inter + 1;
        end
        check_inter = min_inter : (max_inter - min_inter) / 50 : max_inter;
        [N, edges] = histcounts(box_bursty_ticket(:,1), check_inter);
        
        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        plot(edges(1:end-1)*time_grat, N/sum(N), 'k*-', 'linewidth', 1.5)
        set(gca, 'xlim', [0 edges(end-1)*time_grat]);
        set(gca, 'xtick', [0 : 3 : edges(end-1)*time_grat]);
        % set(gca, 'xtick', [0: 3 : edges(end-1)*time_grat/60]);
        xlabel('Duration of Bursty Box Tickets (hour)'); ylabel('PDF');
        %set(gca, 'xscale', 'log');
        title(strcat('RAM, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca, 'fontsize', 24);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'pdf_duration_bursty_box_tickets_',mat2str(ticket_thres(ticket_id)),'_mem_lag_', mat2str(test_lag)));

    end
    
    
    
end




