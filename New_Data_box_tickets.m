% study on the box tickets

close all; clear; clc

load ../New_Data/box_vm_time_series_summary_mem_only.mat
box_vm_time_series_summary_mem = box_vm_time_series_summary;
load ../New_Data/box_vm_time_series_summary_cpu_only.mat

size_box_vm = size(box_vm_time_series_summary);

mkdir('../New_Data_box_tickets_characterization');
path = '../New_Data_box_tickets_characterization/';

ticket_thres = [60, 70, 80];

% Just characterize for no stamp case
BOX_TICKET_CPU = []; BOX_TICKET_MEM = [];
VM_TICKET_CPU = []; VM_TICKET_MEM = [];

BOX_TICKET_VM_TICKET_CPU = {[],[],[]}; BOX_TICKET_VM_TICKET_MEM = {[],[],[]};

possible_related_features = [];

test_lag = 4;

% Consider the time stamp added case
OVERTIME_CPU_BOX_VM = {[],[],[]}; OVERTIME_MEM_BOX_VM = {[],[],[]};

box_num = 1;

for box_id = 1 : size_box_vm(2)
    
    size_box = numel(box_vm_time_series_summary{1, box_id});
    
    % If we don't have time series
    if size_box <= 2
        continue;
    end

    if numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1)) <= 10
        continue;
    end
    
    % First extract the number of tickets for CPU and RAM for different
    % thresholds
    box_cpu = box_vm_time_series_summary{1, box_id}{1,1}(:,4);
    box_mem = box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 4);
    box_ticket_overtime_cpu = [];
    box_ticket_overtime_mem = [];
    for ticket_id = 1 : numel(ticket_thres)
        BOX_TICKET_CPU(box_num, ticket_id) = sum(box_cpu > ticket_thres(ticket_id));
        BOX_TICKET_MEM(box_num, ticket_id) = sum(box_mem > ticket_thres(ticket_id));
        
        box_ticket_overtime_cpu(:, ticket_id) = box_cpu > ticket_thres(ticket_id);
        box_ticket_overtime_mem(:, ticket_id) = box_mem > ticket_thres(ticket_id);
    end  
  
%     % Second look into the possible related features             
%     total_vm_cpu_capacity = 0; used_mean_vm_cpu_capacity = 0;
%     mean_vm_cpu_util_set = []; 
%     
%     total_vm_mem_capacity = 0; used_mean_vm_mem_capacity = 0;
%     mean_vm_mem_util_set = []; 
%     
%     total_box_cpu_capacity = nanmean(box_vm_time_series_summary{1, box_id}{1, 1}(:, 3) ./ ...
%                                      box_vm_time_series_summary{1, box_id}{1, 1}(:, 4) * 100);
%     used_mean_box_cpu_capacity = nanmean(box_vm_time_series_summary{1, box_id}{1, 1}(:, 3));
%     
%     total_box_mem_capacity = nanmean(box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 3));
%     used_mean_box_mem_capacity = nanmean(box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 3) .* ...
%                                      box_vm_time_series_summary_mem{1, box_id}{1, 1}(:, 4) / 100);
    
    vm_ticket_overtime_cpu = {}; vm_ticket_overtime_mem = {};
    for vm_id = 2 : size_box
%         total_vm_cpu_capacity = total_vm_cpu_capacity + nanmean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4) ./ ...
%                                                                 box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5) * 100);
%         used_mean_vm_cpu_capacity = used_mean_vm_cpu_capacity + nanmean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 4));
%         
%         mean_vm_cpu_util_set(end+1) = nanmean(box_vm_time_series_summary{1, box_id}{1, vm_id}(:, 5));
%         
%         total_vm_mem_capacity = total_vm_mem_capacity + nanmean(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 4));
%         used_mean_vm_mem_capacity = used_mean_vm_mem_capacity + nanmean(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 4) .* ...
%                                                                 box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5) / 100);
%         
%         mean_vm_mem_util_set(end+1) = nanmean(box_vm_time_series_summary_mem{1, box_id}{1, vm_id}(:, 5));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%consider VM%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vm_cpu = box_vm_time_series_summary{1, box_id}{1,vm_id}(:,5);
        vm_mem = box_vm_time_series_summary_mem{1, box_id}{1,vm_id}(:, 5);
        VM_TICKET_CPU(box_num, 1:numel(ticket_thres)) = 0;
        VM_TICKET_MEM(box_num, 1:numel(ticket_thres)) = 0;
        for ticket_id = 1 : numel(ticket_thres)
            VM_TICKET_CPU(box_num, ticket_id) = VM_TICKET_CPU(box_num, ticket_id) + sum(vm_cpu > ticket_thres(ticket_id));
            VM_TICKET_MEM(box_num, ticket_id) = VM_TICKET_MEM(box_num, ticket_id) + sum(vm_mem > ticket_thres(ticket_id));

            vm_ticket_overtime_cpu{ticket_id}(:, vm_id - 1) = vm_cpu > ticket_thres(ticket_id);
            vm_ticket_overtime_mem{ticket_id}(:, vm_id - 1) = vm_mem > ticket_thres(ticket_id);          
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % check if box has ticket, how many vms have tickets
    for ticket_id = 1 : numel(ticket_thres)
        % CPU
        box_has_ticket = find(box_ticket_overtime_cpu(:, ticket_id) == 1);
        if sum(box_has_ticket) > 0
            all_lag_sum_vm_has_ticket = [];
            for lag = -test_lag : test_lag
                all_lag_sum_vm_has_ticket(end+1, :) = sum(vm_ticket_overtime_cpu{ticket_id}(max(1,min(box_has_ticket+lag,end)),:)');
            end
            sum_vm_has_ticket = max(all_lag_sum_vm_has_ticket);
            BOX_TICKET_VM_TICKET_CPU{ticket_id}(end+1,:) = [nanmean(sum_vm_has_ticket), nanmedian(sum_vm_has_ticket),...
                                                  prctile(sum_vm_has_ticket, 25), prctile(sum_vm_has_ticket, 75)]; 
        end
        
        % RAM  
        box_has_ticket = find(box_ticket_overtime_mem(:, ticket_id) == 1);
        if sum(box_has_ticket) > 0
            all_lag_sum_vm_has_ticket = [];
            for lag = -test_lag : test_lag
                all_lag_sum_vm_has_ticket(end+1, :) = sum(vm_ticket_overtime_mem{ticket_id}(max(1,min(box_has_ticket+lag,end)),:)');
            end
            sum_vm_has_ticket = max(all_lag_sum_vm_has_ticket);
            BOX_TICKET_VM_TICKET_MEM{ticket_id}(end+1,:) = [nanmean(sum_vm_has_ticket), nanmedian(sum_vm_has_ticket),...
                                                  prctile(sum_vm_has_ticket, 25), prctile(sum_vm_has_ticket, 75)]; 
        end
    end
    
%     mean_vm_cpu_util = nanmean(mean_vm_cpu_util_set);
%     mean_vm_mem_util = nanmean(mean_vm_mem_util_set);
%     
%     vcap_bcap_cpu_ratio = total_vm_cpu_capacity / total_box_cpu_capacity;
%     vusage_vcap_cpu_ratio = min(1, used_mean_vm_cpu_capacity / total_vm_cpu_capacity);
%     vusage_bcap_cpu_ratio = min(1, used_mean_vm_cpu_capacity / total_box_cpu_capacity);
%     vusage_busage_cpu_ratio = min(1, used_mean_vm_cpu_capacity / used_mean_box_cpu_capacity);
%     
%     vcap_bcap_mem_ratio = total_vm_mem_capacity / total_box_mem_capacity;
%     vusage_vcap_mem_ratio = min(1, used_mean_vm_mem_capacity / total_vm_mem_capacity);
%     vusage_bcap_mem_ratio = min(1, used_mean_vm_mem_capacity / total_box_mem_capacity);
%     vusage_busage_mem_ratio = min(1, used_mean_vm_mem_capacity / used_mean_box_mem_capacity);
%     
%     possible_related_features(end+1, :) = [total_vm_cpu_capacity / 1024, used_mean_vm_cpu_capacity / 1024, mean_vm_cpu_util, ...
%                                            total_box_cpu_capacity  / 1024, used_mean_box_cpu_capacity  / 1024, ...
%                                            vcap_bcap_cpu_ratio, vusage_vcap_cpu_ratio, vusage_bcap_cpu_ratio, vusage_busage_cpu_ratio, ...
%                                            VM_TICKET_CPU(box_num, 1:3), ...
%                                            total_vm_mem_capacity  / 1024, used_mean_vm_mem_capacity  / 1024, mean_vm_mem_util, ...
%                                            total_box_mem_capacity  / 1024, used_mean_box_mem_capacity  / 1024, ...
%                                            vcap_bcap_mem_ratio, vusage_vcap_mem_ratio, vusage_bcap_mem_ratio, vusage_busage_mem_ratio, ...
%                                            VM_TICKET_MEM(box_num, 1:3), ...
%                                            size_box - 1];
%                                        
                                       
   box_num = box_num + 1;
    
end

% Check the CDF of box ticket with VM tickets
for ticket_id = 1 : numel(ticket_thres)
    f= {}; x = {}; % mean, median, 25%ile, 75%ile
    fig = figure; set(fig,'Position',[200, 200, 600, 400]);
    for test_id = 1 : numel(BOX_TICKET_VM_TICKET_CPU{ticket_id}(1,:))
        [f{test_id}, x{test_id}] = ecdf(BOX_TICKET_VM_TICKET_CPU{ticket_id}(:, test_id));
        plot(x{test_id}, f{test_id}, 'linewidth',2);
        hold on
    end
    xlabel('Number of VM w/ CPU Tickets','fontsize',18); ylabel('CDF');
    set(gca, 'ylim',[0 1],'fontsize',18); set(gca, 'ytick', [0 : 0.2 : 1]);
    %set(gca, 'yscale', 'log');
    title(strcat('Ticket Threshold = ', {' '}, mat2str(ticket_thres(ticket_id)),'%'), 'fontsize', 18)
    h = legend('Mean', 'Median', '25%ile', '75%ile');
    set(h, 'box', 'on','location','southeast','fontsize',18);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'box_vm_ticket_relation_', mat2str(ticket_id)));
    
    f= {}; x = {}; % mean, median, 25%ile, 75%ile
    fig = figure; set(fig,'Position',[200, 200, 600, 400]);
    for test_id = 1 : numel(BOX_TICKET_VM_TICKET_MEM{ticket_id}(1,:))
        [f{test_id}, x{test_id}] = ecdf(BOX_TICKET_VM_TICKET_MEM{ticket_id}(:, test_id));
        plot(x{test_id}, f{test_id}, 'linewidth',2);
        hold on
    end
    xlabel('Number of VM w/ RAM Tickets','fontsize',18); ylabel('CDF');
    set(gca, 'ylim',[0 1],'fontsize',18); set(gca, 'ytick', [0 : 0.2 : 1]);
    %set(gca, 'yscale', 'log');
    title(strcat('Ticket Threshold = ', {' '}, mat2str(ticket_thres(ticket_id)),'%'), 'fontsize', 18)
    h = legend('Mean', 'Median', '25%ile', '75%ile');
    set(h, 'box', 'on','location','southeast','fontsize',18);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'box_vm_ticket_mem_relation_', mat2str(ticket_id)));
    
end

% boxplot of box ticket (x-axis) and the number of VMs distribution w/
% tickets
test_name = {'mean', 'median', '25_PCTile', '75_PCTile'};
box_size = 10; 
for ticket_id = 1 : numel(ticket_thres)
    box_ticket_idx = BOX_TICKET_CPU(:, ticket_id) > 0;
    box_ticket = BOX_TICKET_CPU(box_ticket_idx, ticket_id);
    box_ticket = ceil(box_ticket / box_size) * box_size;  
    % uniqu box_ticket
    unique_tickets = unique(box_ticket);
    hist_tickets = histcounts(box_ticket, unique_tickets);
    tickets_sum = sum(hist_tickets);
    pdf_tickets = hist_tickets / tickets_sum * 100;
    fig = figure; set(fig,'Position',[200, 200, 600, 400]);
    bar(pdf_tickets);
    xlabel('Percentage of Boxes','fontsize' , 18); ylabel('Number of Boxes tickets');
    set(gca, 'ylim',[0 100],'fontsize' , 18); set(gca, 'ytick', [0 : 10 : 100]);
    %set(gca, 'yscale', 'log');
    title(strcat('Ticket Threshold = ', {' '}, mat2str(ticket_thres(ticket_id)),'%'), 'fontsize', 24)
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'pdf_box_ticket_', mat2str(ticket_id)));
    disp(pdf_tickets);
    for test_id = 1 : numel(BOX_TICKET_VM_TICKET_CPU{ticket_id}(1,:))
        fig = figure; set(fig,'Position',[200, 200, 600, 400]);
        boxplot(BOX_TICKET_VM_TICKET_CPU{ticket_id}(:, test_id), box_ticket);
        xlabel(strcat('Number of BOX CPU Tickets (', test_name{test_id},')'),'fontsize',24); ylabel('Number of VM w/ tickets');
        set(gca, 'ylim',[-1 1],'fontsize',24); %set(gca, 'ytick', [0 : 0.2 : 1]);
        %set(gca, 'yscale', 'log');
        title(strcat('Ticket Threshold = ', {' '}, mat2str(ticket_thres(ticket_id)),'%'), 'fontsize', 24)
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'box_vm_cpu_ticket_boxplot_', mat2str(ticket_id),'_',test_name{test_id}));
    end
    
    
    box_ticket_idx = BOX_TICKET_MEM(:, ticket_id) > 0;
    box_ticket = BOX_TICKET_MEM(box_ticket_idx, ticket_id);
    box_ticket = ceil(box_ticket / box_size) * box_size;    
    unique_tickets = unique(box_ticket);
    hist_tickets = histcounts(box_ticket, unique_tickets);
    tickets_sum = sum(hist_tickets);
    pdf_tickets = hist_tickets / tickets_sum;
    disp(pdf_tickets);
    for test_id = 1 : numel(BOX_TICKET_VM_TICKET_MEM{ticket_id}(1,:))
        fig = figure; set(fig,'Position',[200, 200, 600, 400]);
        boxplot(BOX_TICKET_VM_TICKET_MEM{ticket_id}(:, test_id), box_ticket);
        xlabel(strcat('Number of BOX RAM Tickets (', test_name{test_id},')'),'fontsize',24); ylabel('Number of VM w/ tickets');
        set(gca, 'ylim',[-1 1],'fontsize',24); 
        % set(gca, 'ytick', [0 : 0.2 : 1]);
        %set(gca, 'yscale', 'log');
        title(strcat('Ticket Threshold = ', {' '}, mat2str(ticket_thres(ticket_id)),'%'), 'fontsize', 24)
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'box_vm_mem_ticket_boxplot_', mat2str(ticket_id),'_',test_name{test_id}));
    end
    
end

% labels_name = {'CPU Vcap', 'CPU Vusage', 'CPU Vutil', ...
%                'CPU Bcap', 'CPU Busage', 'CPU Vcap / Bcap', ...
%                'CPU Vusage / Vcap', 'CPU Vuage / Bcap', 'CPU Vusage / Busage', ...
%                'Number of VM Tickets (CPU 60%) ', 'Number of VM Tickets (CPU 70%)', 'Number of VM Tickets (CPU 80%)', ...
%                'MEM Vcap', 'MEM Vusage', 'MEM Vutil', ...
%                'MEM Bcap', 'MEM Busage', 'MEM Vcap / Bcap', ...
%                'MEM Vusage / Vcap', 'MEM Vuage / Bcap', 'MEM Vusage / Busage', ...
%                'Number of VM Tickets (RAM 60%) ', 'Number of VM Tickets (RAM 70%)', 'Number of VM Tickets (RAM 80%)', ...
%                'Number of VMs'};
% 
% all_case_name = {'CPU: BOX TICKETS (60%)','CPU: BOX TICKETS (70%)','CPU: BOX TICKETS (80%)',...
%                  'RAM: BOX TICKETS (60%)','RAM: BOX TICKETS (70%)','RAM: BOX TICKETS (80%)'};
% % Check the correlation among the correaltion value and possible related
% % features
% xcf_cand = [];
% all_case = {BOX_TICKET_CPU(:,1),BOX_TICKET_CPU(:,2),BOX_TICKET_CPU(:,3),...
%             BOX_TICKET_MEM(:,1),BOX_TICKET_MEM(:,2),BOX_TICKET_MEM(:,3)};
% for case_id = 1 : numel(all_case)
%     time_series = all_case{case_id};
%     nan_time_series = find(isnan(time_series));
%     time_series(nan_time_series) = 0;
%     for feature_id = 1 : numel(possible_related_features(1,:))
%         temp = crosscorr(time_series, possible_related_features(:, feature_id), 1);
%         temp(isnan(temp)) = 0;
%         xcf_cand(case_id, feature_id) = temp(2);
%     end 
%     
%     % boxplot the top related features based box plot
%     for cand_id = 1 : numel(possible_related_features(1,:))
%         if cand_id == 10 || cand_id == 11 || cand_id == 12 || cand_id == 22 || cand_id == 23 || cand_id == 24
%             continue;
%         end
%         fig = figure;
%         set(fig, 'Position', [200, 200, 1200, 400]);
%         % [~, cand_id] = max(xcf_cand(case_id, :));
%         % First try: remove the 10%ile and 90%ile
%         lower_bound = prctile(possible_related_features(:, cand_id), 5);
%         upper_bound = prctile(possible_related_features(:, cand_id), 95);
%         
%         used_idx1 = find(possible_related_features(:, cand_id) < upper_bound);
%         used_idx2 = find(possible_related_features(:, cand_id) > lower_bound);
%         used_idx3 = find(time_series > 0);
%         
%         used_idx = intersect(intersect(used_idx1, used_idx2), used_idx3);
%     
%         cand_related_feature = possible_related_features(used_idx, cand_id);
%         cand_time_series = time_series(used_idx);
%         
%         max_feature = max(cand_related_feature);
%         scale = 2;
%         if max_feature > 2
%             cand_related_feature = ceil(log2(cand_related_feature));
%         else
%             interval = max_feature / 10;
%             cand_related_feature = ceil(cand_related_feature / interval);
%             scale = 1;
%         end
%         
%         bucket_range = min(cand_related_feature) : max(cand_related_feature);
% 
%         xticklabel_last = {};
%         if scale == 2
%             for bucket_id = 1 : numel(bucket_range)
%                 xticklabel_last{end+1} = cellstr(num2str(bucket_range(bucket_id), '2^%d'));
%             end
%         else
%             for bucket_id = 1 : numel(bucket_range)
%                 xticklabel_last{end+1} = cellstr(num2str(bucket_range(bucket_id) * interval, '%.2f'));
%             end
%         end
%         
%         xticklabel_name = [xticklabel_last{1}];
%         for name_id = 2 : numel(xticklabel_last)               
%             xticklabel_name(end+1) = strcat(xticklabel_last{name_id-1},'-',xticklabel_last{name_id});
%            
%         end
%         % xticklabel_name = {'1','2','3','4','5'};
%         subplot(1,2,1)
%         boxplot(cand_time_series, cand_related_feature);
%         xlabel(labels_name{cand_id}); ylabel(all_case_name{case_id});
% %         if case_id <= 3
% %             set(gca, 'ylim', [0 20]); set(gca, 'ytick', [0 : 5 : 20], 'fontsize', 15);
% %         else
% %             set(gca, 'ylim', [0 40]); set(gca, 'ytick', [0 : 10 : 40], 'fontsize', 15);
% %         end
%         set(gca, 'xlim', [0 numel(bucket_range)+1], 'fontsize', 15); 
%         set(gca, 'ylim', [0 100]);
%         if scale ~= 2
%             xticklabel_rotate([1: numel(xticklabel_name)],45,xticklabel_name);
%         else
%             set(gca, 'xticklabel', xticklabel_name);
%         end
%         
%         subplot(1,2,2)
%         [X, N] = hist(cand_related_feature, bucket_range);
%         X = 1 : numel(X);
%         bar(X, N/sum(N) * 100);
%         xlabel(labels_name{cand_id}); ylabel('Percentage (%) of Boxes');
%         set(gca, 'xlim', [0 numel(X)+1],'fontsize',15); 
%         if scale ~= 2
%             xticklabel_rotate([1: numel(xticklabel_name)],45,xticklabel_name);
%         else
%             set(gca, 'xticklabel', xticklabel_name);
%         end
%         set(gcf, 'paperpositionmode', 'auto');
%         print('-depsc2','-r300', strcat(path, 'xcf_related_features_boxplot_', mat2str(case_id), '_', mat2str(cand_id)));
%     end
%     close all
% end
% 
% fig = figure;
% grid on
% set(fig,'Position',[200, 200, 1000, 400]);
% plot(xcf_cand(1,:), 'r*-', 'linewidth', 2);
% hold on
% plot(xcf_cand(2,:), 'ko-', 'linewidth', 2);
% hold on
% plot(xcf_cand(3,:), 'm>-', 'linewidth', 2);
% hold on
% plot(xcf_cand(4,:), 'b<-', 'linewidth', 2);
% hold on
% plot(xcf_cand(5,:), 'g^-', 'linewidth', 2);
% hold on
% plot(xcf_cand(6,:), 'y--', 'linewidth', 2);
% hold on
% plot(1:numel(labels_name), ones(numel(labels_name),1)*0.4, 'c-.', 'linewidth', 2)
% hold on
% plot(1:numel(labels_name), ones(numel(labels_name),1)*-0.4, 'c-.', 'linewidth', 2)
% h = legend('CPU: BOX TICKETS (60%)','CPU: BOX TICKETS (70%)','CPU: BOX TICKETS (80%)',...
%            'RAM: BOX TICKETS (60%)','RAM: BOX TICKETS (70%)','RAM: BOX TICKETS (80%)');
% set(h, 'box', 'on', 'location', 'southwest');
% set(gca, 'ylim', [-1 1]); set(gca, 'ytick', [-1 : 0.2 : 1], 'fontsize', 18);
% set(gca, 'xlim', [1 numel(labels_name)]);
% xticklabel_rotate([1 : feature_id], 45, labels_name, 'fontsize', 12);
% ylabel('XCF'); xlabel('Possible Related Feature ID');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path, 'mean_xcf_related_features'));
% 
% % Possibility of box have tickets and ticket distribution
% % Bar plot
% prob_ticket = [numel(find(BOX_TICKET_CPU(:,1) ~=0))/ (box_num -1), numel(find(BOX_TICKET_CPU(:,2) ~=0))/ (box_num -1), ...
%                numel(find(BOX_TICKET_CPU(:,3) ~=0))/ (box_num -1); numel(find(BOX_TICKET_MEM(:,1) ~=0))/ (box_num -1), ...
%                numel(find(BOX_TICKET_MEM(:,2) ~=0))/ (box_num -1), numel(find(BOX_TICKET_MEM(:,3) ~=0))/ (box_num -1)];
% prob_ticket = prob_ticket * 100;
% font_size = 18;
% color = colorgrad(3,'blue_down');
% fig = figure;
% box on
% set(fig, 'Position', [200 200 600 400]);
% hold on
% h1(1,1) = bar(1,prob_ticket(1,1));
% h1(1,2) = bar(1.4,prob_ticket(1,2));
% h1(1,3) = bar(1.8,prob_ticket(1,3));
% h1(2,1) = bar(2.8,prob_ticket(2,1));
% h1(2,2) = bar(3.2,prob_ticket(2,2));
% h1(2,3) = bar(3.6,prob_ticket(2,3));
% set(h1(:,1), 'facecolor', color(1,:));
% set(h1(:,2), 'facecolor', color(2,:));
% set(h1(:,3), 'facecolor', color(3,:));
% set(h1(:,1), 'facecolor', color(1,:), 'barwidth', 0.35); 
% set(h1(:,2), 'facecolor', color(2,:), 'barwidth', 0.35); 
% set(h1(:,3), 'facecolor', color(3,:), 'barwidth', 0.35); 
% 
% set(gca,'ylim', [0 60]); set(gca,'ytick',[0: 10 : 60], 'fontsize', font_size);
% ylabel('Percentage of Boxes (%)', 'fontsize',font_size); % Set the X-axis label
% set(gca, 'xtick', [1.4,3.2]); set(gca, 'xlim', [0.5 4.1]);
% set(gca, 'xticklabel', {'CPU', 'RAM'});
% legen = legend('Ticket Threshold = 60%', 'Ticket Threshold = 70%', 'Ticket Threshold = 80%'); % Add a legend
% set(legen, 'location', 'northeast', 'fontsize',font_size,'box', 'on');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path, 'tickets_probability_different_thres'));
% 
% % check distribution
% mean_ticket_all = {nanzeros(BOX_TICKET_CPU(:,1)), nanzeros(BOX_TICKET_CPU(:,2)), nanzeros(BOX_TICKET_CPU(:,3)), ...
%                nanzeros(BOX_TICKET_MEM(:,1)), nanzeros(BOX_TICKET_MEM(:,2)), nanzeros(BOX_TICKET_MEM(:,3))};  
% mean_ticket = [mean(mean_ticket_all{1}), mean(mean_ticket_all{2}), mean(mean_ticket_all{3}), ...
%                mean(mean_ticket_all{4}), mean(mean_ticket_all{5}), mean(mean_ticket_all{6})];
% std_ticket = [std(mean_ticket_all{1}), std(mean_ticket_all{2}), std(mean_ticket_all{3}), ...
%               std(mean_ticket_all{4}), std(mean_ticket_all{5}), std(mean_ticket_all{6})];
%            
% x_tick = [1,1.4,1.8,2.8,3.2,3.6];
% 
% fig = figure;
% box on
% set(fig, 'Position', [200 200 600 400]);
% hold on
% h1(1,1) = bar(1,mean_ticket(1,1));
% h1(1,2) = bar(1.4,mean_ticket(1,2));
% h1(1,3) = bar(1.8,mean_ticket(1,3));
% h1(2,1) = bar(2.8,mean_ticket(1,4));
% h1(2,2) = bar(3.2,mean_ticket(1,5));
% h1(2,3) = bar(3.6,mean_ticket(1,6));
% set(h1(:,1), 'facecolor', color(1,:));
% set(h1(:,2), 'facecolor', color(2,:));
% set(h1(:,3), 'facecolor', color(3,:));
% set(h1(:,1), 'facecolor', color(1,:), 'barwidth', 0.35); 
% set(h1(:,2), 'facecolor', color(2,:), 'barwidth', 0.35); 
% set(h1(:,3), 'facecolor', color(3,:), 'barwidth', 0.35); 
% 
% errorbar(x_tick, mean_ticket, std_ticket, 'r.', 'linewidth',1);
% 
% set(gca,'ylim', [0 120]); set(gca,'ytick',[0: 20 : 120]);
% ylabel('Number of Tickets', 'fontsize',font_size); % Set the X-axis label
% set(gca, 'xtick', [1.4,3.2]); set(gca, 'xlim', [0.5 4.1], 'fontsize',font_size);
% set(gca, 'xticklabel', {'CPU', 'RAM'}, 'fontsize',font_size);
% legen = legend('Ticket Threshold = 60%', 'Ticket Threshold = 70%', 'Ticket Threshold = 80%'); % Add a legend
% set(legen, 'location', 'northwest', 'fontsize',font_size,'box', 'on');
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path, 'tickets_distribution_different_thres'));