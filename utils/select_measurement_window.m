function [start,stop] = select_measurement_window(plot_data,f)
%SELECT_MEASUREMENT_WINDOW allows the selection of a subsection of a
%measurement.
%
%[start,stop] = SELECT_MEASUREMENT_WINDOW(plotdata,f)
%It prompts the user to select a start and stop time on figure f where:
%       plot_data = the concentration data in a timetable format.
%       f = the figure where the concentration data is plotted and the
%           measurement window is to be selected from.
%       start = measurement start time selected by the user.
%       stop = measurement stop time selected by the user.

    if isempty(f) == 1;f = gcf;else;end
    figure(f)
    data_frequency_min = median([plot_data.TIME;0]-[0;plot_data.TIME],"omitmissing");
    [Start_min,y1] = ginput; % Asks user to manually select start time from the plot.
    % Finds the couple of data points around this selected time. 
    %(The time selected will rarely correlate to an exact point in the Raw Data. The Picarro takes a data point every .0094 minutes)
    while isempty(Start_min) == 1 % This prompts user to select the start time again because the input didn't work the first time (bad mouse click for example).
        disp("Error: Select Start Time Again")
        [Start_min,y1] = ginput;
    end
    %If the comand window freezes or displays nothing and another time is
    %selected it will produce a vector with two rows instead of a scalar.
    if length(Start_min)>1
        Start_min = Start_min(1);
        disp("First Start chosen was selected")
    end
    K1 = find(plot_data.TIME>Start_min-data_frequency_min/2&plot_data.TIME<Start_min+data_frequency_min/2); % Indexes the points around the selected time.
    s = 0;
    while isempty(K1) == 1
        K1 = find(plot_data.TIME>Start_min-(data_frequency_min+s)&plot_data.TIME<Start_min+(data_frequency_min+s));
        s = s + data_frequency_min;
    end
    K = find(K1,1); % Choses the first of these indexes (out of 1-3 numbers).
    K = K1(K); % Selects index that will be used to select start time from the Raw Data File.
    start = plot_data.time(K); % Selects start time from datetime column in Raw Data
    disp('Start Selected') % Displys the selected start time to the user. 
    disp(start)

    [Stop_min,y2] = ginput; % Same as above but for stop time.
    while isempty(Stop_min) == 1
        disp("Error: Select Stop Time Again")
        [Stop_min,y1] = ginput;
    end
    if length(Stop_min)>1
        Stop_min = Stop_min(1);
        disp("First Start chosen was selected")
    end
    L1 = find(plot_data.TIME>Stop_min-data_frequency_min/2&plot_data.TIME<Stop_min+data_frequency_min/2);
    s = 0;
    while isempty(L1) == 1
        L1 = find(plot_data.TIME>Stop_min-(data_frequency_min+s)&plot_data.TIME<Stop_min+(data_frequency_min+s));
        s = s + data_frequency_min;
    end
    L = find(L1,1);
    L = L1(L);
    stop = plot_data.time(L);
    disp('Stop Selected')
    disp(stop)
end

