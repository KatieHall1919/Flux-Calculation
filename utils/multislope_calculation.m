function [slope_multi,sigma_multi_std,min_slope,max_slope] = multislope_calculation(x_time,y_data)
%MULTISLOPE_CALCULATION calculates a representative dX/dt for the gas of
%interest by calculating slopes of 1 minute intervals stepping every 10
%seconds, averaging these, and calculating the standard error of the slopes
%to aquire an uncertainty for dX/dt.
%
%[slope_multi,sigma_multi_std,min_slope,max_slope] = multislope_calculation(x_time,y_data)
%   x_time = time data in datetime or duration format (e.g. timetable.TIME)
%   y_data = gas concentration data 
%   slope_multi = representative slope of gas concentration over time of
%       measurement
%   sigma_multi_std = representative uncertainty of the slope ("slope_multi")
%   min_slope,max_slope = minimum and maximum of the individual 1-min slopes
%       calculated
        range = 1; step = 10/60; % in minutes
        % calculate slope over 1 min windows
        tstart = minutes(x_time(1));tstop = tstart + range;
        slope = [];sigma = [];tslope = [];intercept = []; % Preallocation
        % clf(figure(5));figure(5)
        % scatter(minutes(x_time),y_data)
        % hold on
        while tstop < minutes(x_time(end))
            % calculate slope for the current minute
            I = find(minutes(x_time)>=tstart & minutes(x_time)<tstop);
            line = fitlm(minutes(x_time(I)),y_data(I));
            slope = [slope,table2array(line.Coefficients(2,1))];
            sigma = [sigma,table2array(line.Coefficients(2,2))*2];
            intercept = [intercept,table2array(line.Coefficients(1,1))];
            tslope = [tslope,(tstart+tstop)/2];
            % increment forward 10 seconds
            tstart = tstart + step;
            tstop = tstop + step;
            % % can be used to plot individual slopes
            % b = table2array(line.Coefficients(1,1));
            % m = table2array(line.Coefficients(2,1));
            % y_plot = m*minutes(x_time(I))+b;
            % x_plot = minutes(x_time(I));
            % plot(x_plot,y_plot)
            
        end
        % calculate slope
        slope_multi = mean(slope);
        sigma_multi_std = (std(slope)*2)/sqrt(length(slope));
        % maximum and minumum slope calculated
        max_slope = max(slope); min_slope = min(slope);
end

