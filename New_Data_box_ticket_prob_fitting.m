% this script aims at charactering BOX tickets per box

close all; clear; clc

load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros

path = '../New_Data_box_tickets_fit/';
mkdir(path);

mkdir(strcat(path, 'Fig/'));

ticket_thres = [60, 70, 80];

time_grat = 900;
day_point = 96;

PROB_BOX_TICKET = [];
BOX_TICKET_FEATURE = [];

box_num = 1;
size_box_vm = size(box_vm_time_series_summary);

test_lag = 0;

fuzzy_bound = 0;

sig_vm_thres = 0.8;

USED_BOX = 0;

for box_id = 1 : size_box_vm(2)
    
    % feature 0: box size 
    size_box = numel(box_vm_time_series_summary{1, box_id});

    % If we don't have time series
    if size_box < 2
        continue;
    end
    
%     if box_num >= 1000
%         break;
%     end

    % at least we should have 3 days data
    total_points = numel(box_vm_time_series_summary{1, box_id}{1,1}(:,1));
    if total_points < day_point * 3
        continue;
    end

    % First extract the number of tickets for CPU and RAM for different
    % thresholds
    USED_BOX = USED_BOX + 1;
    
    box = box_vm_time_series_summary{1, box_id}{1,1}(:,4);
    box_demands = box_vm_time_series_summary{1, box_id}{1,1}(:,3);
    
    % feature 1: capaciy of box
    box_capacity = nanmean(box_demands ./ box * 100);
    
    % feature set 2: mean usage and std of usage
    mean_box_usage = nanmean(box); std_box_usage = nanstd(box);
    median_box_usage = nanmedian(box); 
    box_usage_25 = prctile(box, 25); box_usage_75 = prctile(box, 75);
    
    % feature set 3: the vm total v_cap
    vm_cap = []; vm_on = []; vm_cap_on = [];
    for vm_id = 2 : size_box
        non_zero_idx = box_vm_time_series_summary{1,box_id}{vm_id}(:, 4) ~= 0;
        vm_on(:, end+1) = non_zero_idx;
        vm_demands = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx, 4);
        vm_usage = box_vm_time_series_summary{1, box_id}{vm_id}(non_zero_idx,5);
        vm_cap(end+1) = nanmean(vm_demands ./ vm_usage * 100);
        
        vm_cap_on(:,end+1) = vm_on(:, end) * vm_cap(end); 
    end
    sum_vcap = sum(vm_cap);
    
    vm_on_sum = sum(vm_on');
    vm_on_cap_sum = sum(vm_cap_on');
    
    mean_on_vm = nanmean(vm_on_sum); std_on_vm = nanstd(vm_on_sum);
    median_on_vm = nanmedian(vm_on_sum); 
    on_vm_25 = prctile(vm_on_sum, 25); on_vm_75 = prctile(vm_on_sum, 75);
    
    mean_on_cap_vm = nanmean(vm_on_cap_sum); std_on_cap_vm = nanstd(vm_on_cap_sum);
    median_on_cap_vm = nanmedian(vm_on_cap_sum); 
    on_cap_vm_25 = prctile(vm_on_cap_sum, 25); on_cap_vm_75 = prctile(vm_on_cap_sum, 75);
      
    for ticket_id = 1 : numel(ticket_thres)
        box_ticket = box > ticket_thres(ticket_id);
        sum_ticket = sum(box_ticket);
        % RECORD THE PROBABILITY OF BOX
        PROB_BOX_TICKET(end+1) = sum_ticket / total_points;
        BOX_TICKET_FEATURE(:, end+1) = [ticket_thres(ticket_id); ...
                                       size_box - 1; box_capacity; sum_vcap;...
                                       mean_box_usage; std_box_usage; ....
                                       median_box_usage; box_usage_25; box_usage_75; ...
                                       mean_on_vm; std_on_vm; ...
                                       median_on_vm; on_vm_25; on_vm_75; ...
                                       mean_on_cap_vm; std_on_cap_vm; ...
                                       median_on_cap_vm; on_cap_vm_25; on_cap_vm_75];
    end
end

save(strcat(path, 'PROB_BOX_TICKET'), 'PROB_BOX_TICKET');
save(strcat(path, 'BOX_TICKET_FEATURE'), 'BOX_TICKET_FEATURE');

disp('The training begins');

% choose the interested freatures
interet_idx = [1, 5, 6, 7, 8, 9];
BOX_TICKET_FEATURE = BOX_TICKET_FEATURE(interet_idx, :);

% First, generate the scatter plot among features;
N = numel(BOX_TICKET_FEATURE(1,:)) / 3; 
mean_usage = BOX_TICKET_FEATURE(2,1:N);
y_name = {'std','median','25_pct','75_pct'};
Y_NAME = {'STD of Box Usage (%)','Median Box Usage (%)','25%ile of Box Usage (%)',...
          '75%ile of Box Usage (%)'};
for feature_id = 3 : numel(interet_idx)
    % first do the polynomial fitting
    if feature_id > 3
        degree = 1;
    else
        degree = 2;
%         fcmdata = [mean_usage', BOX_TICKET_FEATURE(feature_id, 1:N)'];
%         [centers,U] = fcm(fcmdata,2);
%         maxU = max(U);
%         index1 = find(U(1,:) == maxU); index2 = find(U(2,:) == maxU);
%         fig = figure;
%         font_size = 12;
%         set(fig, 'Position', [200 200 300 200]);
%         plot(fcmdata(index1,1),fcmdata(index1,2),'ob')
%         hold on
%         plot(fcmdata(index2,1),fcmdata(index2,2),'or')
%         plot(centers(1,1),centers(1,2),'xb','MarkerSize',15,'LineWidth',3)
%         plot(centers(2,1),centers(2,2),'xr','MarkerSize',15,'LineWidth',3)
%         
%         xlabel('Mean Box Usage (%)'); ylabel(Y_NAME{feature_id - 2});
%         h = legend('Cluster1','Cluster2','Centroid1','Centroid2'); 
%         set(h, 'location', 'northwest');
%         set(gca, 'fontsize', font_size);
%         set(gcf, 'paperpositionmode', 'auto');
%         print('-depsc2','-r300', strcat(path, 'Fig/fuzzy_cluster_', y_name{feature_id-2}));
    end
    p = polyfit(mean_usage, BOX_TICKET_FEATURE(feature_id, 1:N), degree);
    fit_val = polyval(p, mean_usage);
        
    fig = figure;
    font_size = 12;
    set(fig, 'Position', [200 200 300 200]);
    scatter(mean_usage, BOX_TICKET_FEATURE(feature_id, 1:N))
    hold on
    plot(mean_usage, fit_val, 'r*', 'linewidth', 1.5)
    set(gca, 'xlim', [0 100]); set(gca, 'xtick', [0 : 10 : 100]);
    if feature_id == 3
        set(gca, 'ylim', [0 30]); set(gca, 'ytick', [0 : 3 : 30]);
    else
        set(gca, 'ylim', [0 100]); set(gca, 'ytick', [0 : 10 : 100]);
    end
    xlabel('Mean Box Usage (%)'); ylabel(Y_NAME{feature_id - 2});
    h = legend('Original','Fitting'); set(h, 'location', 'northwest');
    set(gca, 'fontsize', font_size);
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'Fig/relation_mean_usage_with_', y_name{feature_id-2}));
end

% decision tree
rtree = fitrtree(BOX_TICKET_FEATURE', PROB_BOX_TICKET);
imp = predictorImportance(rtree);
imp = imp ./ sum(imp);
resuberror = resubLoss(rtree);
prob_fit = predict(rtree, BOX_TICKET_FEATURE');
error = PROB_BOX_TICKET - prob_fit';
nanzeros_idx = PROB_BOX_TICKET ~= 0;
abs_pct_error = abs(error(nanzeros_idx)) ./ PROB_BOX_TICKET(nanzeros_idx);

disp({'resuberror is ', resuberror});
disp({'Fitting error is ', nanmean(error)});
disp({'Fitting percentage error is ', nanmean(abs_pct_error)});

% use the 10-fold cross-validation to test the prediction errors
cvrtree = crossval(rtree);
cvloss = kfoldLoss(cvrtree);

N = numel(PROB_BOX_TICKET); K = 10;
Indices = crossvalind('Kfold', N, K);

split_test = {}; split_test_feature = {};
test_error = {}; test_abs_pct_error = {};
fit_error = {}; fit_abs_pct_error = {};

cross_imp = []; 
for fold_id = 1 : K
    idx = Indices == fold_id;
    split_test{fold_id} = PROB_BOX_TICKET(idx);
    split_test_feature{fold_id} = BOX_TICKET_FEATURE(:,idx);
    
    % fit error summary
    rtree = fitrtree(BOX_TICKET_FEATURE(:, ~idx)', PROB_BOX_TICKET(~idx));
    cross_imp(fold_id, :) = predictorImportance(rtree);
    cross_imp(fold_id, :) = cross_imp(fold_id, :) ./ sum(cross_imp(fold_id, :));
    prob_fit = predict(rtree, BOX_TICKET_FEATURE(:, ~idx)');
    fit_error{end+1} = PROB_BOX_TICKET(~idx) - prob_fit';
    traing_prob = PROB_BOX_TICKET(~idx);
    nanzeros_idx = traing_prob ~= 0;
    fit_abs_pct_error{end+1} = abs(fit_error{end}(nanzeros_idx)) ./ traing_prob(nanzeros_idx);
    
    % test the one-fold data
    prob_predict = predict(rtree, BOX_TICKET_FEATURE(:, idx)');
    test_error{end+1} = PROB_BOX_TICKET(idx) - prob_predict';
    nanzeros_idx = split_test{fold_id} ~= 0;
    test_abs_pct_error{end+1} = abs(test_error{end}(nanzeros_idx)) ./ split_test{fold_id}(nanzeros_idx); 
    
    disp({'Fitting error is ', nanmean(fit_error{end})});
    disp({'Fitting percentage error is ', nanmean(fit_abs_pct_error{end})});
    disp({'Test error is ', nanmean(test_error{end})});
    disp({'Test percentage error is ', nanmean(test_abs_pct_error{end})});
    
    fig = figure;
    set(fig, 'Position', [200 200 900 200])
    subplot(1,3,1)
    prob_size = [0 : 0.01 : 1];
    font_size = 10;
    [N, edges] = histcounts(PROB_BOX_TICKET(~idx), prob_size);
    [N_test, edges_test] = histcounts(PROB_BOX_TICKET(idx), prob_size);
    % h = bar(N ./ sum(N));
    % set(gca, 'xtick', [1:numel(edges)]);
    % set(gca, 'xlim', [1 numel(edges)]);
    % set(gca, 'xticklabel', edges(1:end));
    plot(edges(1:min(end, numel(N))), N/sum(N), 'k-', 'linewidth',2);
    hold on
    plot(edges_test(1:min(end, numel(N_test))), N_test/sum(N_test), 'm-.', 'linewidth',2);
    set(gca, 'xlim', [edges(1) edges(end)]);
    % set(gca, 'xscale', 'log');
    ylabel('PDF'); xlabel('Prob(box w/ tickets)');
    h = legend('Training Set', 'Validation Set');
    set(h, 'location', 'northwest');
    set(gca, 'fontsize', font_size)

    % plot the fitting error
    font_size = 10;
    subplot(1,3,2)
%     error_size = [-0.001 : 0.002/ 1000 : 0.001];
%     [N, edges] = histcounts(fit_error{end}, error_size);
%     [N_test, edges_test] = histcounts(test_error{end}, error_size);
    [f, x] = ecdf(abs(fit_error{end}));
    [f_test, x_test] = ecdf(abs(test_error{end}));
    plot(x * 100, f, 'k-', 'linewidth',2);
    hold on
    plot(x_test * 100, f_test, 'm-.', 'linewidth',2);
    % set(gca, 'xtick', [1:numel(edges)]);
    set(gca, 'xlim', [0 2]);
    ylabel('CDF'); xlabel('Error (%)');
    h = legend(strcat('Fitting: Mean = ', mat2str(round(nanmean(abs(fit_error{end})*100), 4)), '%'), ...
               strcat('Validation: Mean = ', mat2str(round(nanmean(abs(test_error{end})*100), 4)), '%'));
    set(h, 'location', 'southeast');
    set(gca, 'fontsize', font_size)

    subplot(1,3,3)
    [f, x] = ecdf(fit_abs_pct_error{end});
    [f_test, x_test] = ecdf(test_abs_pct_error{end});
    plot(x * 100, f, 'k-', 'linewidth',2);
    hold on 
    plot(x_test * 100, f_test, 'm-.', 'linewidth',2);
    % set(gca, 'xtick', [1:numel(edges)]);
    set(gca, 'xlim', [0 100]);
    ylabel('CDF'); xlabel('Abs PCT Error (%)');
    h = legend(strcat('Fitting: Mean = ', mat2str(round(nanmean(fit_abs_pct_error{end}*100), 2)), '%'), ...
               strcat('Validation: Mean = ', mat2str(round(nanmean(test_abs_pct_error{end}*100), 2)), '%'));
    set(h, 'location', 'southeast');
    set(gcf, 'paperpositionmode', 'auto');
    print('-depsc2','-r300', strcat(path, 'Fig/part_feature_cross_validation_fold_', mat2str(fold_id)));
end

disp({'cvloss is ', cvloss});

% Plot figures
% Figure1: PDF of PROB(BOX have a ticket) 
fig = figure;
set(fig, 'Position', [200 200 900 200])
subplot(1,3,1)
prob_size = [0 : 0.01 : 1];
font_size = 10;
[N, edges] = histcounts(PROB_BOX_TICKET, prob_size);
% h = bar(N ./ sum(N));
% set(gca, 'xtick', [1:numel(edges)]);
% set(gca, 'xlim', [1 numel(edges)]);
% set(gca, 'xticklabel', edges(1:end));
plot(edges(1:min(end, numel(N))), N/sum(N), 'k-', 'linewidth',2);
% set(gca, 'xscale', 'log');
ylabel('PDF'); xlabel('Prob(box w/ tickets)');
h = legend('All as Training Set');
set(h, 'location', 'northwest');
set(gca, 'fontsize', font_size)

% plot the fitting error
font_size = 10;
subplot(1,3,2)
[f, x] = ecdf(error);
plot(x * 100, f, 'k-', 'linewidth',2);
% set(gca, 'xtick', [1:numel(edges)]);
set(gca, 'xlim', [0 2]);
ylabel('CDF'); xlabel('Error (%)');
h = legend(strcat('Fitting: Mean = ', mat2str(round(nanmean(abs(error)*100), 4)), '%'));
set(h, 'location', 'southeast');
set(gca, 'fontsize', font_size)

subplot(1,3,3)
[f, x] = ecdf(abs_pct_error);
plot(x * 100, f, 'k-', 'linewidth',2);
% set(gca, 'xtick', [1:numel(edges)]);
set(gca, 'xlim', [0 100]);
ylabel('PDF'); xlabel('Fitting Abs PCT Error (%)');
h = legend(strcat('Mean = ', {' '}, mat2str(round(nanmean(x*100), 2)), '%'));
set(h, 'location', 'southeast');
set(gcf, 'paperpositionmode', 'auto');
print('-depsc2','-r300', strcat(path, 'Fig/part_feature_pdf_fit_all_error'));


