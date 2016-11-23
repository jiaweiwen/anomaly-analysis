function [rtree, error, ape, imp] = DecisionTree(X, Y)

    rtree = fitrtree(X, Y);
    imp = predictorImportance(rtree);
    imp = imp ./ sum(imp);
    resuberror = resubLoss(rtree);
    prob_fit = predict(rtree, X);
    error = Y - prob_fit';
    ape = abs(error) ./ Y;
    ape(isnan(ape)) = 0;
    ape(isinf(ape)) = 0;

    disp({'resuberror is ', resuberror});
    disp({'Fitting error is ', nanmean(error)});
    disp({'Fitting percentage error is ', nanmean(ape)});
end

