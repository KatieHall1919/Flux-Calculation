function [ConcentrationData] = pull_data(filepath,foldername,instrument)
%%PULL_DATA loads all concentration data from a particular date into a timetable
%
%[ConcentrationData] = pull_data(filepath,foldername,instrument)
%   foldername = folder where the data is located (must end in '\').
%   filepath = filepath/directory of the folder defined by 'foldername'
%   instrument = instrument from which the data is from, options: "LGR" or "Picarro"
%   ConcentrationData = Gas concentration and remaining data from the analyzer given in a timetable format

    foldername = strcat(filepath,foldername);

    if instrument == "Picarro"      
        a = dir(sprintf('%s/*.dat',foldername)); % Creates a structure based on the directory to read in filenames of the raw Picarro data files. Only takes files ending in .dat in the folder.
        filename = strings; % Creates an empty variable for preallocation (increases speed/accuracy).
        for p=1:length(a);filename(p)=convertCharsToStrings(a(p).name);end % Populates variable filename with names (of data files) from the structure a (based on the names in the file directory)). 
        
        for j= 1:length(filename)
            data_holder = readtable(strcat(foldername,filename(j)));  % Reads in the Picarro data file as a table.
            % Creates a datetime object from the DATE and TIME fields of the raw
            % data table.
            Datestring = string(data_holder.DATE);
            Timestring = string(data_holder.TIME);
            date_and_time = strcat(Datestring, " ", Timestring); % Adds the date and time into one string value.
            %infmt = 'yyyy-MM-dd HH:mm:ss.SSS'; % Input format to convert from string to datetime format.
            infmt = 'dd-MMM-yyyy HH:mm:ss.SSS';
            data_holder.time = datetime(date_and_time,'Inputformat',infmt); % Creates a table called data_holder with the Raw Picarro Data from the jth file.
            
            % Rearrange array to put the string date and time at the end, with the
            % datetime object in the beginning.
            data_holder = [data_holder(:,29) data_holder(:,3:28) data_holder(:,1:2)];
            
            % Appends Raw Picarro data tables from multiple data files.
            if j == 1
                ConcentrationData = data_holder;
            else
                ConcentrationData = vertcat(ConcentrationData,data_holder); % this appends data from the currently read file to the data from previous files using vertical concatenation.
            end
        end 
        ConcentrationData = table2timetable(ConcentrationData); % Converts Picarrodata into a timetable.  
    elseif instrument == "LGR"
        a = dir(sprintf('%s/*.txt',foldername)); % Creates a structure based on the directory to read in filenames. Only takes files ending in .dat in the folder.
        filename = strings; % Creates an empty variable for perallocation (increases speed).
        for i=1:length(a);filename(i)=convertCharsToStrings(a(i).name);end % Populates variable filename with names (of data files) from the structure a (based on the file directory)). 
        
        for j= 1:length(filename)
            % Pull the text files from the folder using the file names from the
            % previous step.
            fileID1 = fopen((strcat(foldername,filename(j))),'r'); % Opens the text file so the code can "read it" ('r'->read).
            c = textscan(fileID1,'%s','delimiter','\n');
            
            startRow = find(contains(c{1,1},'Time,'))+1;
            endRow = find(contains(c{1,1},'-----BEGIN PGP MESSAGE-----')) - 1 ;
            if isempty(endRow) == 1
                disp('PGP Message Missing in Data File')
            end
            fclose(fileID1);
            
            % Specify the format of the data table. All columns '%f' for floating
            % number except for the first column, which is a datetime format. 
            formatSpec = '%{MM/dd/yy HH:mm:ss.SSS}D%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            
            fileID2 = fopen((strcat(foldername,filename(j))),'r');
            textscan(fileID2, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
            dataArray = textscan(fileID2, formatSpec, endRow-startRow+1, 'Delimiter', ',', 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID2);
            
            LGRdata1 = table(dataArray{1:end-1}, 'VariableNames', {'time','CH4_ppm','CH4_ppm_se','H2O_ppm','H2O_ppm_se','CO2_ppm','CO2_ppm_se','GasP_torr','GasP_torr_se','GasT_C','GasT_C_se','AmbT_C','AmbT_C_se','RD0_us','RD0_us_se','RD1_us','RD1_us_se','Fit_Flag'});
            data_holder = LGRdata1;
        
            if j == 1
                ConcentrationData = data_holder;
            else
                ConcentrationData = vertcat(ConcentrationData,data_holder); % this appends data from the currently read file to the data from previous files using vertical concatenation.
            end
        end
        ConcentrationData = table2timetable(ConcentrationData); % Converts LGRdata into a timetable.
        % NOTE: in the timetable, TIME = duration format, time = datetime format
        
        % calculate the derivative, add nmin_avg if changing averaging period
        %LGRdata = LGR_dXdt(LGRdata, 1);
        
        %% Add a time series column to the data table
        % Convert to Fraction of Days since JAN1
        % Months
        M = month(ConcentrationData.time);
        % Days
        d = day(ConcentrationData.time);
        % Hours
        H = hour(ConcentrationData.time);
        % Minutes
        m = minute(ConcentrationData.time);
        % Seconds
        s = second(ConcentrationData.time);
        
        % % Fraction of minutes since beginning of day of data file
        % FRAC_MINS = (H)*60 + m + s/60;
        % ConcentrationData.FRAC_MINS = FRAC_MINS;
        
        ConcentrationData.TIME = duration(H,m,s);
    else
        error("Error: Instrument options are LGR or Picarro.")
    end
end

