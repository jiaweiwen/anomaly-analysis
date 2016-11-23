close all; clear; clc

load ../New_Data_7days/box_vm_time_series_summary_cpu_only_with_zeros

path = '../New_Data_box_moments_fit/Fig/';

num_box = numel(box_vm_time_series_summary);

ticket_thres = 70; tail_target = 95;

coeff_all = [];

day_thres = 6;
test_day = 2;

for box_id = 1 : num_box
    
    if numel(box_vm_time_series_summary{box_id}) == 0
        continue;
    end
    
    test_usage = box_vm_time_series_summary{box_id}{1}(:,end);
    
    if numel(test_usage) < 96 * day_thres
        continue;
    end
    
    test_ticket = test_usage >= ticket_thres;
    
    if prctile(test_usage(1:(day_thres-test_day)*96), tail_target) >= ticket_thres % sum(test_ticket(1:(day_thres-test_day)*96)) > 1
        coeff_all(end+1, 1) = 1;
    else
        coeff_all(end+1, 1) = 0;
    end
    
    if prctile(test_usage((day_thres-test_day)*96+1 : day_thres * 96), tail_target) >= ticket_thres % sum(test_ticket((day_thres-test_day)*96+1 : day_thres * 96)) > 1
        coeff_all(end, 2) = 1;
    else
        coeff_all(end, 2) = 0;
    end
%     coeff_usage = autocorr(test_usage, 96 * 4);
%     coeff_ticket = autocorr(test_ticket, 96 * 4);    
%     coeff_ticket(isnan(coeff_ticket)) = 1;
    
    
%     coeff_all(end+1, 1:2) = [max(abs(coeff_usage(96:end))), max(abs(coeff_ticket(96:end)))];  
end

idx_zero = find(coeff_all(:,1) == 0);
idx_zero_zero = find(coeff_all(idx_zero, 2) == 0);
disp(strcat('Prob(NO TICKET | NO TICKET) = ', mat2str(numel(idx_zero_zero) / numel(idx_zero))));
disp(strcat('Prob(TICKET | NO TICKET) = ', mat2str(1 - numel(idx_zero_zero) / numel(idx_zero))));

idx_one = find(coeff_all(:,1) == 1);
idx_one_zero = find(coeff_all(idx_one, 2) == 0);
disp(strcat('Prob(NO TICKET | TICKET) = ', mat2str(numel(idx_one_zero) / numel(idx_one))));
disp(strcat('Prob(TICKET | TICKET) = ', mat2str(1 - numel(idx_one_zero) / numel(idx_one))));

% font_size = 15;
% 
% fig = figure;
% set(fig, 'Position', [200 200 400 300]);
% [f1,x1] = ecdf(coeff_all(:,1));
% [f2,x2] = ecdf(coeff_all(:,2));
% plot(x1, f1, 'r-', 'linewidth', 2);
% hold on
% plot(x2, f2, 'k--', 'linewidth', 2);
% xlabel('Maximum Long-term (> 1 day) ACF'); ylabel('CDF');
% set(gca, 'xlim',[0 1]); set(gca, 'xtick',[0 : 0.1 : 1]);
% set(gca, 'ylim',[0 1]); set(gca, 'ytick',[0 : 0.1 : 1]);
% h = legend('Actual Series', 'Ticket Series'); set(h,'location', 'southeast')
% set(gca, 'fontsize', font_size);
% set(gcf, 'paperpositionmode', 'auto');
% print('-depsc2','-r300', strcat(path, 'cdf_long_term_cpu_usage_acf'));
