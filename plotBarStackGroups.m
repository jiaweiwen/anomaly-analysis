function [h] = plotBarStackGroups(stackData, groupLabels)
%% Plot a set of stacked bars, but group them according to labels provided.
%%
%% Params: 
%%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)
%%
%% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu
%%
%% 
NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);
NumStacksElements = size(stackData, 3);

% Ji's setting
cmap = cell(NumStacksElements, NumStacksPerGroup);
pre_defined_color = {'blue_down',  'red_down','red_down','green_down'};
for idx = 1 : NumStacksElements
    temp = colorgrad(NumStacksPerGroup, pre_defined_color{idx});
    for j = 1 : NumStacksPerGroup
        cmap{idx,j} = temp(j, :);
    end
end


% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
legend_name = [{'CPU-Original'};{'RAM-Original'};{'CPU-DTW'};{'RAM-DTW'};{'CPU-CBC'};{'RAM-CBC'}];
fig = figure;
set(fig, 'Position', [200 200 600 400]);
box on;
    hold on; 
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    
    h(i,:) = bar(Y, 'stacked');
    for ele_idx = 1 : NumStacksElements
        set(h(i,ele_idx),'facecolor',cmap{ele_idx, i},'edgecolor','k'); 
    end   
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
    
end
hold off;
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);
end 
