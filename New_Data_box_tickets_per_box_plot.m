% plot the box ticket characterization
close all; clear; clc

path = '../New_Data_box_vm_tickets_characterization_per_box/';

load (strcat(path, 'USED_BOX'));
load (strcat(path, 'BOX_TICKET_ALL'));
load (strcat(path, 'BOX_TICKET_ALL_GROUPED'));
load (strcat(path, 'BOX_TICKET_INTER_ARRIVAL'));
load (strcat(path, 'BOX_TICKET_INTER_ARRIVAL_ACF'));

load (strcat(path, 'BOX_TICKET_INTER_ARRIVAL_GROUPED'));
load (strcat(path, 'BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED'));
load (strcat(path, 'BOX_TICKET_DURATION'));

load (strcat(path, 'PROB_BOX_VM_TICKET'));
load (strcat(path, 'PROB_BOX_VM_TICKET_DIFF_NUM'));
load (strcat(path, 'BOX_TICKET_VM_NO_TICKET'));
load (strcat(path, 'BOX_TICKET_VM_TICKET'));

mkdir(strcat(path, 'Fig/'));

ticket_thres = [60, 70, 80];
box_size = {};
ticket_counts = {}; ticket_counts_grouped = {};
box_w_ticket_prob = [];
for ticket_id = 1 : numel(ticket_thres)
    used_idx = BOX_TICKET_ALL{ticket_id}(:,3) ~= 0;  
    ticket_counts{ticket_id} = BOX_TICKET_ALL{ticket_id}(used_idx,3);
    ticket_counts_grouped{ticket_id}= BOX_TICKET_ALL_GROUPED{ticket_id}(:,3);
    box_size{ticket_id} = BOX_TICKET_ALL{ticket_id}(used_idx,4);
    box_w_ticket_prob(ticket_id) = sum(used_idx) / USED_BOX;
end

%%%%%%%%%% Fig 0: bubble plot of probabity that box has tickets %%%%%%%%%%
first_try = true;
if first_try
    for ticket_id = 1 : numel(ticket_thres)
        prob_case = PROB_BOX_VM_TICKET_DIFF_NUM{1, ticket_id};
        prob_size = 0 : 0.02 : 1.02;
        vm_ratio = 10:10:100;
        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        for vm_num = 1 : numel(prob_case(1,:))
            [N, edges] = histcounts(prob_case(:,vm_num), prob_size);
            N = max(1, N);
            scatter(ones(numel(N),1)*vm_ratio(vm_num), edges(1:end-1), N, 'linewidth', 2);
            hold on
        end
        set(gca, 'ylim', [0 1]); set(gca, 'ytick', [0:0.1:1]); 
        set(gca, 'xlim', [0 110]); set(gca,'xtick', [0:10:100]);
        xlabel('PCT of VMs w/ Tickets (%)'); ylabel('Prob(Box has tickets)');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/box_ticket_w_different_vm_thres_', mat2str(ticket_thres(ticket_id))));
    end
end

%%%%%%%%%%%%%% Fig 1: Number of Tickets per day per box%%%%%%%%%%%%%%%%%%%%
% Fig 1.1 Overall Distribution
second_try = false;
if second_try 

    fig = figure;
    set(fig, 'Position', [200 200 600 400]);
    aboxplot(ticket_counts);
    set(gca,'xticklabel',{' '});
    set(gca,'ylim', [0 20]); set(gca,'ytick',[0: 4 : 20]); 
    xlabel('CPU'); ylabel('No. of BOX Tickets per BOX per Day');
    %set(gca, 'xscale', 'log');
    %title('No Grouping');
    h = legend('Ticket Threshold = 60%', 'Ticket Threshold = 70%','Ticket Threshold = 80%');
    set(gca, 'fontsize', 20);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'Fig/box_counts_overall_distribution'));

    fig = figure;
    set(fig, 'Position', [200 200 600 400]);
    aboxplot(ticket_counts_grouped);
    set(gca,'xticklabel',{' '});
    set(gca,'ylim', [0 5]); set(gca,'ytick',[0: 1 : 5]); 
    xlabel('CPU'); ylabel('No. of BOX Tickets per BOX per Day');
    %set(gca, 'xscale', 'log');
    %title('Bursty Grouped');
    h = legend('Ticket Threshold = 60%', 'Ticket Threshold = 70%','Ticket Threshold = 80%');
    set(gca, 'fontsize', 20);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'Fig/box_counts_grouped_overall_distribution'));

    fig = figure;
    box on
    color = [[0.5 0.5 1];[0.8 0 0];[0 0.7 0]];
    set(fig, 'Position', [200 200 600 400]);
    hold on
    for ticket_id = 1 : numel(ticket_thres)
        h = bar(ticket_id, box_w_ticket_prob(ticket_id));
        set(h, 'facecolor', color(ticket_id, :), 'barwidth', 0.5);
    end
    set(gca,'xticklabel',{' '});
    set(gca,'ylim', [0 0.3]); set(gca,'ytick',[0: 0.1 : 0.3]); 
    xlabel('CPU'); ylabel('PCT. of BOX w/ Tickets');
    %set(gca, 'xscale', 'log');
    %%title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
    h = legend('Ticket Threshold = 60%', 'Ticket Threshold = 70%','Ticket Threshold = 80%');
    set(gca, 'fontsize', 20);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'Fig/box_w_ticket_prob'));

    % plot mean duration for each ticket 
    for ticket_id = 1 : numel(ticket_thres)
        fig = figure;
        box on
        set(fig, 'Position', [200 200 600 400]);
        [f1,x1] = ecdf(BOX_TICKET_DURATION{ticket_id}(:,1));
        [f2,x2] = ecdf(BOX_TICKET_DURATION{ticket_id}(:,3));
        [f3,x3] = ecdf(BOX_TICKET_DURATION{ticket_id}(:,4));
        [f4,x4] = ecdf(BOX_TICKET_DURATION{ticket_id}(:,5));
        [f5,x5] = ecdf(BOX_TICKET_DURATION{ticket_id}(:,6));
        plot(x1, f1, 'k-', 'linewidth', 2);
        hold on
        plot(x2, f2, 'r--', 'linewidth', 2);
        hold on
        plot(x3, f3, 'b--', 'linewidth', 2);
        hold on
        plot(x4, f4, 'm-.', 'linewidth', 2);
        hold on
        plot(x5, f5, 'g-.', 'linewidth', 2);
        h = legend('Mean','5%ile','95%ile','25%ile','75%ile');
        set(h, 'location', 'southeast');
        xlabel('Bursty Ticket Duration * 15 (min)');
        ylabel('CDF');
        set(gca, 'xlim', [0 100]);
        set(gca, 'ytick',[0:0.1:1]); set(gca, 'ylim', [0 1]);
        set(gca, 'xscale', 'log');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/box_ticket_duration_thres_', mat2str(ticket_thres(ticket_id))));
    end
    
    % Fig 1.2 For each ticket threshold, get PDF and box plots for size 
    pre_edge = [2, 8, 16, 32, 64, 128];

    for ticket_id = 1 : numel(ticket_thres)
        [N, edges, bin] = histcounts(box_size{ticket_id}, pre_edge);

        fig = figure;
        set(fig, 'Position', [200 200 600 400]);   
        bar(N/sum(N), 'grouped');
        set(gca, 'xlim',[0 numel(edges)-0.1]); set(gca, 'xtick',[0 : 1 : numel(edges)-0.1]);
        set(gca, 'xticklabel', [0, edges(1:end-1)]);
        set(gca, 'ylim', [0 0.35]); set(gca, 'ytick',[0:0.05:0.35]);
        ylabel('PDF'); xlabel('Number of VMs');
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/size_pdf_thres_', mat2str(ticket_thres(ticket_id))));

        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        temp_box = pre_edge(bin);
        unique_box = unique(temp_box);
        mean_num = [];
        for i = 1 : numel(unique_box)
            idx = find(temp_box == unique_box(i));
            test_box = ticket_counts{ticket_id}(idx);
            % remove outliers
            test_box_upper = prctile(test_box, 95); upper_idx = find(test_box < test_box_upper);
            test_box_lower = prctile(test_box, 5); lower_idx = find(test_box > test_box_lower);
            test_idx = intersect(upper_idx, lower_idx);
            mean_num(i) = mean(test_box(test_idx));
        end

        boxplot(ticket_counts{ticket_id}, temp_box);
        hold on
        plot(1:numel(unique_box), mean_num, 'ko', 'markersize', 5, 'markerfacecolor','m');
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        if ticket_id == 1
            set(gca,'ylim', [0 35]); set(gca,'ytick',[0: 5 : 35]); 
        else
            set(gca,'ylim', [0 25]); set(gca,'ytick',[0: 5 : 25]); 
        end
        ylabel('Number of Tickets per BOX per Day'); xlabel('Number of VMs');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/box_ticket_dist_thres_', mat2str(ticket_thres(ticket_id))));
    end

    %%%%%%%%%%%%%%%%%%%%%%%Fig 2: Conditional Probability%%%%%%%%%%%%%%%%%%%%%%
    for ticket_id = 1 : numel(ticket_thres)
        [N, edges, bin] = histcounts(PROB_BOX_VM_TICKET{ticket_id}(:,3), pre_edge);

        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        temp_box = pre_edge(bin);
        unique_box = unique(temp_box);
        mean_num = [];
        for i = 1 : numel(unique_box)
            idx = find(temp_box == unique_box(i));
            test_box = PROB_BOX_VM_TICKET{ticket_id}(idx,1);
            % remove outliers
            test_box_upper = prctile(test_box, 95); upper_idx = find(test_box < test_box_upper);
            test_box_lower = prctile(test_box, 5); lower_idx = find(test_box > test_box_lower);
            test_idx = intersect(upper_idx, lower_idx);
            mean_num(i) = mean(test_box(test_idx));
        end
        boxplot(PROB_BOX_VM_TICKET{ticket_id}(:,1), temp_box);
        hold on
        plot(1:numel(unique_box), mean_num, 'ko', 'markersize', 5, 'markerfacecolor','m');
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca,'ylim', [0 1]); set(gca,'ytick',[0: 0.2 : 1]); 
        ylabel('Prob(vm | box)'); xlabel('Number of VMs');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/prob_vm_box_thres_', mat2str(ticket_thres(ticket_id))));

        fig = figure;
        set(fig, 'Position', [200 200 600 400]);
        temp_box = pre_edge(bin);
        unique_box = unique(temp_box);
        mean_num = [];
        for i = 1 : numel(unique_box)
            idx = find(temp_box == unique_box(i));
            test_box = PROB_BOX_VM_TICKET{ticket_id}(idx,2);
            % remove outliers
            test_box_upper = prctile(test_box, 95); upper_idx = find(test_box < test_box_upper);
            test_box_lower = prctile(test_box, 5); lower_idx = find(test_box > test_box_lower);
            test_idx = intersect(upper_idx, lower_idx);
            mean_num(i) = mean(test_box(test_idx));
        end
        boxplot(PROB_BOX_VM_TICKET{ticket_id}(:,2), temp_box);
        hold on
        plot(1:numel(unique_box), mean_num, 'ko', 'markersize', 5, 'markerfacecolor','m');
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca,'ylim', [0 0.5]); set(gca,'ytick',[0: 0.1 : 0.5]); 
        ylabel('Prob(box | vm)'); xlabel('Number of VMs');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/prob_box_vm_thres_', mat2str(ticket_thres(ticket_id))));

    end
end

%%%%%%%%%%%%%%%%%%%%%% Fig 3: bubble plot of %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Number of VM has ticket v.s. Prob(box has ticket)%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% Fig 4: Box has tickets and No VM has %%%%%%%%%%%%%%%%%%
% Fig 4.1: Check the overal utilization among these VM
third_try = false;
if third_try 
    for ticket_id = 1 : numel(ticket_thres)
        test_box = BOX_TICKET_VM_NO_TICKET{1, ticket_id};
        test_box_vm = BOX_TICKET_VM_TICKET{1, ticket_id};

        [f, x] = ecdf(test_box(:,end-2));
        [f_vm, x_vm] = ecdf(test_box_vm(:,end-2));

        % mean util and std util across vm 
        fig = figure;
        set(fig, 'Position', [200 200 1200 400])
        subplot(1,2,1)
        plot(x, f, 'k-', 'linewidth',2)
        hold on
        plot(x_vm, f_vm, 'r--', 'linewidth',2)
        % title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca,'xlim', [0 200]); set(gca,'xtick',[0: 20 : 200]); 
        set(gca,'ylim', [0 1]); set(gca,'ytick',[0: .1 : 1]); 
        h = legend('Box Ticket w/o VM Tickets','Box Ticket w/ VM Tickets');
        set(h, 'location', 'southeast');
        xlabel('Max/Min VM Usage'); ylabel('CDF');
        set(gca, 'fontsize', 20);

        [f, x] = ecdf(test_box(:,end));
        [f_vm, x_vm] = ecdf(test_box_vm(:,end));
        subplot(1,2,2)
        plot(x, f, 'k-', 'linewidth',2)
        hold on
        plot(x_vm, f_vm, 'r--', 'linewidth',2)
        % title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca,'xlim', [0 100]); set(gca,'xtick',[0: 10 : 100]); 
        set(gca,'ylim', [0 1]); set(gca,'ytick',[0: .1 : 1]); 
        h = legend('Box Ticket w/o VM Tickets','Box Ticket w/ VM Tickets');
        set(h, 'location', 'southeast');
        xlabel('95%ile/5%ile VM Usage'); ylabel('CDF');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/maxmin_ratio_thres_', mat2str(ticket_thres(ticket_id))));  

        % demand ratio figure
        bin_size = 0.01:0.01:1.01;
        fig = figure;
        set(fig, 'Position', [200 200 1200 400])
        subplot(1,2,1);
        [N, edges, bin] = histcounts(test_box(:,3), bin_size);
        [N_vm, edges_vm, bin_vm] = histcounts(test_box_vm(:,3), bin_size);
        plot(edges(1:end-1),N/sum(N), 'k-', 'linewidth',2)
        hold on
        plot(edges_vm(1:end-1),N_vm/sum(N_vm), 'r--', 'linewidth',2)
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca,'xlim', [0 1]); set(gca,'xtick',[0: 0.1 : 1]); 
        h = legend('Box Ticket w/o VM Tickets','Box Ticket w/ VM Tickets');
        set(h, 'location', 'northwest');
        % set(gca,'ylim', [0 0.08]); set(gca,'ytick',[0: .02 : .08]); 
        xlabel('Demands Ratio'); ylabel('PDF');
        set(gca, 'fontsize', 20);

        subplot(1,2,2);
        [f,x] = ecdf(test_box(:,3));
        [f_vm,x_vm] = ecdf(test_box_vm(:,3));
        plot(x, f, 'k-', 'linewidth',2);
        hold on
        plot(x_vm, f_vm, 'r--', 'linewidth',2);
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca,'xlim', [0 1]); set(gca,'xtick',[0: 0.1 : 1]); 
        set(gca,'ylim', [0 1]); set(gca,'ytick',[0: .1 : 1]); 
        xlabel('Demands Ratio'); ylabel('CDF');
        h = legend('Box Ticket w/o VM Tickets','Box Ticket w/ VM Tickets');
        set(h, 'location', 'northwest');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/interest_demand_ratio_thres_', mat2str(ticket_thres(ticket_id))));  

        % dominant VM ratio figure
        vm_ratio = test_box(:,5) ./ test_box(:,6) * 100;
        vm_ratio_vm = test_box_vm(:,5) ./ test_box_vm(:,6) * 100;
        fig = figure;
        set(fig, 'Position', [200 200 1200 400])
        subplot(1,2,1)
        bin_size = 10:10:110;
        [N, edges, bin] = histcounts(vm_ratio, bin_size);
        [N_vm, edges_vm, bin_vm] = histcounts(vm_ratio_vm, bin_size);
        h = bar([(N/sum(N))', (N_vm/sum(N_vm))']);
        set(h(1), 'facecolor', 'r');
        set(h(2), 'facecolor', 'b');
        h = legend('Box Ticket w/o VM Tickets','Box Ticket w/ VM Tickets');
        set(h, 'location', 'northwest');
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        % set(gca,'xlim', [0 100]); set(gca,'xtick',[0: 10 : 100]); 
        set(gca,'xticklabel', bin_size(1:end-1)); 
        set(gca,'ylim', [0 0.5]); set(gca,'ytick',[0: .05 : .5]);
        xlabel('Dominant VM Ratio (%)'); ylabel('PDF');
        set(gca, 'fontsize', 20);

        % dominant VM figure
        subplot(1,2,2)
        bin_size = [2,3,4,5,6,8,16,32,64,128];
        [N, edges, bin] = histcounts(test_box(:,5), bin_size);
        [N_vm, edges_vm, bin_vm] = histcounts(test_box_vm(:,5), bin_size);
        unique_edge = union(unique(edges),unique(edges_vm));
        h = bar([(N/sum(N))', (N_vm/sum(N_vm))']);
        set(h(1), 'facecolor', 'r');
        set(h(2), 'facecolor', 'b');

        h = legend('Box Ticket w/o VM Tickets','Box Ticket w/ VM Tickets');
        set(h, 'location', 'northwest');
        %title(strcat('CPU, Thres = ', {' '}, mat2str(ticket_thres(ticket_id)), '%'));
        set(gca,'xticklabel', unique_edge(1:end-1)); 
        set(gca,'ylim', [0 0.5]); set(gca,'ytick',[0: .05 : .5]); 
        xlabel('Dominant VM Number'); ylabel('PDF');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/interest_dominant_vm_thres_', mat2str(ticket_thres(ticket_id))));  
    end
end

%%%%%%%%%%%%%%%%%%% Fig 5: plot the inter-arrival time %%%%%%%%%%%%%%%%%%%%
% Fig 5.1: plot the maximum autocorrelation of the inter-arrival time
forth_try = false;
if forth_try
    for ticket_id = 1 : numel(ticket_thres)
        acf_inter_arrival = BOX_TICKET_INTER_ARRIVAL_ACF{ticket_id};
        acf_inter_arrival_grouped = BOX_TICKET_INTER_ARRIVAL_ACF_GROUPED{ticket_id};
        [f, x] = ecdf(acf_inter_arrival(:,1));
        [f_grouped, x_grouped] = ecdf(acf_inter_arrival_grouped(:,1));
        fig = figure; 
        set(fig, 'Position', [200 200 600 400]);
        plot(x, f, 'k-', 'linewidth',2);
        hold on
        plot(x_grouped, f_grouped, 'r--', 'linewidth',2);
        h = legend('No Group','Grouped');
        set(h, 'location', 'northeast');
        set(gca,'xlim', [0 1]); set(gca,'xtick',[0: 0.1 : 1]); 
        set(gca,'ylim', [0 1]); set(gca,'ytick',[0: 0.1 : 1]);
        xlabel('Max ACF of Inter-arrival Time'); ylabel('CDF');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/max_acf_inter_arrival_cdf_', mat2str(ticket_thres(ticket_id))));  

        fig = figure; 
        set(fig, 'Position', [200 200 600 400]);
        boxplot(acf_inter_arrival(:,1), acf_inter_arrival(:,2));
        set(gca,'ylim', [0 1]); set(gca,'ytick',[0: 0.1 : 1]);
        xlabel('Lag (Number of Box Tickets)'); ylabel('ACF');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/acf_no_group_inter_arrival_', mat2str(ticket_thres(ticket_id))));  

        fig = figure; 
        set(fig, 'Position', [200 200 600 400]);
        boxplot(acf_inter_arrival_grouped(:,1), acf_inter_arrival_grouped(:,2));
        set(gca,'ylim', [0 1]); set(gca,'ytick',[0: 0.1 : 1]);
        xlabel('Lag (Number of Box Bursty Tickets)'); ylabel('ACF');
        set(gca, 'fontsize', 20);
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path, 'Fig/acf_inter_arrival_', mat2str(ticket_thres(ticket_id))));  
    end
end

close all




