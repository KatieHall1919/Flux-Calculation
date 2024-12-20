function [Data,win,Filter_Times,Percent_Filtered] = filter_data(plot_data,filtering,filtering_limit_s)
%FILTER_DATA allows the user to interactively select locations of data
%points which need to be filtered out.
%
%[Data,win,Filter_Times,Percent_Filtered] = filter_data(plot_data,filtering,filtering_limit_s)
%   plot_data = Gas concentration and instrument data in a timetable format
%   filtering = Determines whether filtering will happen ("Y" or "N")
%   filtering_limit_s = Limits the amount of data that can be taken out
%       from a measurement in seconds
%   Data = Filtered gas concentration and instrument data in a timetable format
%   win = number of data points around a selected point which were removed
%       from the data
%   Filter_Times = Times during the measurement that were filtered out
%       given in a list within a string
%   Percent_Filtered = Percent of the data points removed (# of data points
%       after filtering/# of data points before filtering), units = %
    filtering1 = filtering;
    Data = plot_data;
    win = NaN;Filter_Times = NaN;
    Percent_Filtered = NaN;
    while filtering1 == 'Y'
        clf(figure(2))
    %--PLOT UNFILTERED DATA----------------------------------------------------
        f2 = figure(2); set(f2,'color','w');
        hold on
        subplot(311);plot(plot_data.TIME,plot_data.CH4_ppm,'b');ylabel('CH_4 (ppm)');title('Unfiltered');
        subplot(312);plot(plot_data.TIME,plot_data.CavityPressure_torr,'b');ylabel('Cavity Pressure (torr)');
        subplot(313);plot(plot_data.TIME,plot_data.H2O_perc,'b');xlabel('time (min)');ylabel('H2O (%)');
        hold off

        data_frequency_min = median([plot_data.TIME;0]-[0;plot_data.TIME],"omitmissing");
    %--SELECT FILTERING---------------------------------------------          
        % Filter Window: number of points either side of flagged data to be removed (set to 0 to only remove flagged data), 1 sec time interval of data = 0.5
        win1 = input("\nFilter Window (default = 5 for Picarro, 1 for LGR): ");
        if win1 == 0;win1 = 5;end
        win = win1;
        
        Filter_Selected = [];
        done = 'N';
        while done == 'N'
            [Filter_min,y1] = ginput;
            while isempty(Filter_min) == 1 % This prompts user to select the start time again because the input didn't work the first time (bad mouse click for example).
                disp("Error: Select Filter Time Again")
                [Filter_min,y1] = ginput;
            end
            %If the comand window freezes or displays nothing and another time is
            %selected it will produce a vector with two rows instead of a scalar.
            if length(Filter_min)>1
                Filter_min = Filter_min(1);
                disp("First was selected")
            end
            K1 = find(plot_data.TIME>Filter_min-data_frequency_min/2&plot_data.TIME<Filter_min+data_frequency_min/2); % Indexes the points around the selected time.
            while isempty(K1) == 1
                K1 = find(plot_data.TIME>Filter_min-data_frequency_min&plot_data.TIME<Filter_min+data_frequency_min);
            end
            K = find(K1,1); % Choses the first of these indexes (out of 1-3 numbers).
            K = K1(K); % Selects index that will be used to select start time from the Raw Data File.
            Filter_Selected = [Filter_Selected; plot_data.TIME(K)]; % Selects start time from datetime column in Raw Data
            c = input("\nDone? \nEnter 'N' to select another \npoint to filter out.\n(Y/N): ", "s");
            done = c; 
        end
        for j = 1:length(Filter_Selected)
            if j == 1
                Filter_Times1 = string(Filter_Selected(j));
            else
                Filter_Times1 = Filter_Times1 + "," + string(Filter_Selected(j));
            end
        end

        Filter_Times = Filter_Times1;

    %--FILTER DATA-------------------------------------------------------------
        I = [];
        % locations of flagged data
        for k = 1:length(Filter_Selected)
            I = [I;find(Filter_Selected(k) == plot_data.TIME)];
        end

        % Windows around flagged data
        Iminus = I-win1; Iplus = I+win1;
        
        % Locations of data within windows
        Itmp = []; 
        for k = 1:length(I) % Finds locations of indices within windows around flagged data. 
            Itmp = [Itmp Iminus(k):Iplus(k)]; 
        end
        Icut = unique(Itmp); Icut = Icut(Icut>0); % Selects unique positive index values to cut out.
        Icut_cut = find(Icut>length(plot_data.TIME)); % If filtering index is bigger than index of selected window.
        Icut(:,Icut_cut) = []; % Cuts out indices higher that the data in the selected time window.

        % Cut bad data from table
        PlotData_filtered = plot_data;
        PlotData_filtered(Icut,:) = []; % Cuts out selected "bad" data.
        
    %--PLOT FILTERED DATA------------------------------------------------------
        clf(figure(2))
        f2 = figure(2); set(f2,'color','w');
        hold on
        subplot(311);plot(PlotData_filtered.TIME,PlotData_filtered.CH4_ppm, 'r');ylabel('CH_4 (ppm)');title('Filtered');
        subplot(312);plot(PlotData_filtered.TIME,PlotData_filtered.CavityPressure_torr, 'r');ylabel('Cavity Pressure (torr)');
        subplot(313);plot(PlotData_filtered.TIME,PlotData_filtered.H2O_perc, 'r');ylabel('H2O (%)');
        hold off
       
        % Sets a limit of 30s on the amount of continuous data that can
        % be deleted. This will prompt the user to refilter the data.
        filtering_toolittledata = 'Z'; 
        for j = 2:length(PlotData_filtered.TIME)
            time_dif=(PlotData_filtered.TIME(j)-PlotData_filtered.TIME(j-1));
            if seconds(time_dif) >= filtering_limit_s
                disp("Too much data taken out, Refilter")
                filtering_toolittledata = 'X';
            else
            end
        end
        if filtering_toolittledata == 'X'
            filtering1 = 'Y';
        else
            % Ask the user whether filtering is still need:
            a = input("\nNeeds Filtering? \nEnter 'Y' to redo all filtering\nfor this sample.\n(Y/N): ", "s");
            filtering1 = a;   
        end
        clear Data
        Data = PlotData_filtered;
        Percent_Filtered = ((height(plot_data)-(height(PlotData_filtered)))/height(PlotData_filtered))*100;
    end
end

