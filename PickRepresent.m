function [error_summary, time_series_rank, final_error_summary, ...
          final_coeff_determine, final_new_time_series_set] = PickRepresent(time_series_set, max_lag)
     % Step 1: calculate the XCF with differnet lags    
    num_time_series = numel(time_series_set);
    xcf_set = ones(num_time_series);
    lag_set = zeros(num_time_series);
    for time_idx_1 = 1 : num_time_series
        for time_idx_2 = time_idx_1 + 1 : num_time_series
            time_series1 = time_series_set{time_idx_1};
            time_series2 = time_series_set{time_idx_2};
            [xcf, lags, bounds] = crosscorr(time_series1, time_series2, max_lag);
            [xcf_set(time_idx_1, time_idx_2), lag_idx] = max(abs(xcf));
            lag_set(time_idx_1, time_idx_2) = lags(lag_idx);
            xcf_set(time_idx_2, time_idx_1) = xcf_set(time_idx_1, time_idx_2);
            lag_set(time_idx_2, time_idx_1) = -lags(lag_idx);
        end
    end

    % Step 2: ranking time series based on the xcf
    [max_val, rank_idx] = max(mean(xcf_set));
    % store the original time series index
    time_series_rank = [rank_idx];
    new_time_series = {time_series_set{rank_idx}};
    original_xcf_set = xcf_set;
    for time_idx = 2 : num_time_series
        xcf_set(rank_idx,:) = 0;
        xcf_set(:,rank_idx) = 0;

        [max_val, rank_idx] = max(mean(xcf_set));
        time_series_rank(end+1) = rank_idx;
        new_time_series{end+1} = time_series_set{rank_idx};
    end

    % Step 3: fitting with candidate time series
    error_summary = [];
    for time_idx = 1 : num_time_series - 1
        x_time_idx = time_series_rank(1:time_idx);
        error_current = [];
        for fit_idx = time_idx+1 : num_time_series
            original_time_idx = time_series_rank(fit_idx);
            candidate_lag = xcf_set(x_time_idx, original_time_idx);
            left_lag = min(candidate_lag); right_lag = max(candidate_lag);

            if left_lag >= 0 && right_lag >= 0
                Y = new_time_series{fit_idx}(1 + right_lag : end);
                len_y = numel(Y);
                X = ones(len_y, 1);
                for x_idx = 1 : time_idx
                    begin_idx = 1 + right_lag - candidate_lag(x_idx);
                    end_idx = begin_idx + len_y - 1;
                    X(:, end+1) = new_time_series{x_idx}(begin_idx:end_idx);
                end

            elseif left_lag < 0 && right_lag < 0
                Y = new_time_series{fit_idx}(1 : end + left_lag);
                X = ones(len_y, 1);
                for x_idx = 1 : time_idx
                    begin_idx = 1 - candidate_lag(x_idx);
                    end_idx = begin_idx + len_y - 1;
                    X(:, end+1) = new_time_series{x_idx}(begin_idx:end_idx);
                end

            else
                Y = new_time_series{fit_idx}(1 + right_lag : end + left_lag);
                len_y = numel(Y);
                X = ones(len_y, 1);
                for x_idx = 1 : time_idx
                    begin_idx = 1 + right_lag - candidate_lag(x_idx);
                    end_idx = begin_idx + len_y - 1;                 
                    X(:, end+1) = new_time_series{x_idx}(begin_idx:end_idx);
                end                    
            end

            [b, bint, r, rint, stats] = regress(Y, X);
            % The absolute percentage error is what we care about
            ape = abs(r ./ Y);
            error_current(end+1) = mean(ape);  
        end
        error_summary(end+1,1:2) = [mean(error_current), max(error_current)];
    end
    error_summary(end+1,1:2) = 0;

    % Final Step: Pick up the representatives based on the lower bound of
    % the convex (this is really coarse-try :) later, we could try to find
    % a more robust method)
    if numel(error_summary(:,1)) - 1 <= 2
        k_mean = [1]; k_worse = [1];
    else
        x = [1 : 1 : numel(error_summary(:,1))-1];
        k_mean = convhull(x, error_summary(1:end-1,1)');
        k_worse = convhull(x, error_summary(1:end-1,2)');
    end
    
    for k_mean_idx = 1 : numel(k_mean) - 1
        if k_mean(k_mean_idx+1) < k_mean(k_mean_idx)
            break;
        end
    end
    
    for k_worse_idx = 1 : numel(k_worse) - 1
        if k_worse(k_worse_idx+1) < k_worse(k_worse_idx)
            break;
        end
    end
    
    % Store the index of time series in the time_series_set
    final_time_series_set = union(time_series_rank(k_mean(1:k_mean_idx)), ...
                                  time_series_rank(k_worse(1:k_worse_idx)));
    final_time_series_set = [final_time_series_set, time_series_rank(end)];
    left_time_series_set = setdiff(time_series_rank, final_time_series_set);
%     disp(time_series_rank);
%     disp(k_worse);
    disp(final_time_series_set);
    
    % Store the index of time series in the new_time_series
    final_new_time_series_set = union(k_mean(1:k_mean_idx),k_worse(1:k_worse_idx));
    final_new_time_series_set = [final_new_time_series_set', numel(new_time_series)];
    
    % Fitting: to check the error and the r^2

    error_current = []; coeff_current = [];
    for fit_idx = 1 : numel(left_time_series_set)
        original_time_idx = left_time_series_set(fit_idx);
        candidate_lag = xcf_set(final_time_series_set, original_time_idx);
        left_lag = min(candidate_lag); right_lag = max(candidate_lag);

        if left_lag >= 0 && right_lag >= 0
            Y = time_series_set{original_time_idx}(1 + right_lag : end);
            len_y = numel(Y);
            X = ones(len_y, 1);
            for x_idx = 1 : numel(final_time_series_set)
                begin_idx = 1 + right_lag - candidate_lag(x_idx);
                end_idx = begin_idx + len_y - 1;
                X(:, end+1) = time_series_set{final_time_series_set(x_idx)}(begin_idx:end_idx);
            end

        elseif left_lag < 0 && right_lag < 0
            Y = time_series_set{original_time_idx}(1 : end + left_lag);
            X = ones(len_y, 1);
            for x_idx = 1 : numel(final_time_series_set)
                begin_idx = 1 - candidate_lag(x_idx);
                end_idx = begin_idx + len_y - 1;
                X(:, end+1) = time_series_set{final_time_series_set(x_idx)}(begin_idx:end_idx);
            end

        else
            Y = time_series_set{original_time_idx}(1 + right_lag : end + left_lag);
            len_y = numel(Y);
            X = ones(len_y, 1);
            for x_idx = 1 : numel(final_time_series_set)
                begin_idx = 1 + right_lag - candidate_lag(x_idx);
                end_idx = begin_idx + len_y - 1;                 
                X(:, end+1) = time_series_set{final_time_series_set(x_idx)}(begin_idx:end_idx);
            end                    
        end

        [b, bint, r, rint, stats] = regress(Y, X);
        % The absolute percentage error is what we care about
        ape = abs(r ./ Y);
        error_current(end+1) = mean(ape);  
        
        ssresid = sum(r.^2);
        sstotal = (length(Y)-1) * var(Y);
        coeff_current(end+1) = 1 - ssresid/sstotal;
    end
    final_error_summary = [mean(error_current), max(error_current)];
    final_coeff_determine = [mean(coeff_current), min(coeff_current)];

end

