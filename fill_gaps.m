function [time_series_without_gaps] = fill_gaps(time_series_metric, grat)
% This function is used to fill in the time holes in the time series
    time_series_without_gaps = []; 
    temp = time_series_metric; len = numel(temp(1,:));
    time_begin = temp(1,1); time_end = temp(end,1); t = 1;
    time_lag = 0; total_lag = (time_end-time_begin)/grat;
    while time_lag <= total_lag
        if time_lag*grat+time_begin  == temp(t,1)
            time_series_without_gaps(time_lag+1,:) = temp(t,:);
        else
            while temp(t,1) > time_lag*grat+time_begin
                time_series_without_gaps(time_lag+1,1) = time_lag*grat+time_begin;
                time_series_without_gaps(time_lag+1,2) = temp(1,2);
                time_series_without_gaps(time_lag+1,3:len) = -1;
                time_lag = time_lag + 1;
            end
            time_series_without_gaps(time_lag+1,:) = temp(t,:);
        end
        t = t + 1;
        time_lag = time_lag + 1;
    end
end

