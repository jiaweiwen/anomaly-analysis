% check the conditional probability
close all; clear; clc

path = '../New_Data_box_vm_tickets_characterization_per_box_per_tenant_with_zeros/';

load (strcat(path, 'PROB_BOX_VM_TICKET'))

ticket_thres = [60, 70, 80];

box_size_thres = 20;

chosen_tenant_id = [4, 5, 8, 10];

median_prob = {zeros(numel(chosen_tenant_id), 8), zeros(numel(chosen_tenant_id), 8),...
               zeros(numel(chosen_tenant_id), 8)};

for ticket_id = 1 : numel(ticket_thres)
    for tenant_idx = 1 : numel(chosen_tenant_id)       
        tenant_id = chosen_tenant_id(tenant_idx);
        
        fig = figure;
        set(fig, 'Position', [200 200 300 200]);
        f = {}; x = {}; 
        color = {'r', 'b', 'm', 'k'};
        line_style = {'-', '--', '-.', ':'};
        for col_idx = 1 : 4
            [f{col_idx}, x{col_idx}] = ecdf(PROB_BOX_VM_TICKET{ticket_id}{tenant_id}(:, col_idx));
            median_prob{ticket_id}(tenant_idx, (col_idx-1)*2+1 : col_idx*2) = ...
                [nanmedian(PROB_BOX_VM_TICKET{ticket_id}{tenant_id}(:, col_idx)), ...
                 1 - nanmedian(PROB_BOX_VM_TICKET{ticket_id}{tenant_id}(:, col_idx))];
            plot(x{col_idx}, f{col_idx}, 'color', color{col_idx}, 'linestyle', line_style{col_idx})
            hold on
        end
        h = legend('Prob(VM | Box)', 'Prob(VM | Box No)', ...
                   'Prob(Box | VM)', 'Prob(Box | VM No)');
        set(h, 'location', 'southeast');
        ylabel('CDF'); xlabel('Conditional Probability');
        title(strcat('Tenant ID = ', {' '}, mat2str(tenant_id), ...
              ', Ticket Threshold = ', {' '}, mat2str(ticket_thres(ticket_id)),'%'));
        set(gcf, 'paperpositionmode', 'auto');
        print('-depsc2','-r300', strcat(path,'tenant_', mat2str(tenant_id), ...
                               '_thres_', mat2str(ticket_thres(ticket_id))));
    end
end