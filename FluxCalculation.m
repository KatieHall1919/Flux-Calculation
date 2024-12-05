%% Calculate Flux
% This script is published in conjunction with and should be credited with:
% "Natural Geologic Methane Emissions from Microseepage in the Michigan 
% Basin are Negligible” by Kathleen R. Hall, Thomas Weber, Marika Stock, 
% Marc Buursink, Haoran Piao, Mingzhe Zhu, Katey Walter Anthony, 
% Vasilii Petrenko

clear
clf(figure(1));clf(figure(3));

%% Inputs
foldername = '06Aug2021_Picarro\'; %ddMMMyyyy

instrument = "Picarro";
chamber_type = "Bucket_2";

date = "06/08/2021"; %dd/MM/yyyy
sample = ["1","6"]; % Must be a string: "-your sample-". Make it an empty vector to select all samples in that day.
sample_range = []; % Must be a string. Must be in order of the table. Calls all samples in between. sample_range may not work if sample numbers are not unique for the day.

% load concentration data
filepath = [pwd,'\data\'];addpath(filepath);
ConcentrationData = pull_data(filepath,foldername,instrument);
fprintf('\nUploaded Raw %s Data for %s\n\n',instrument,date);

if instrument == "Picarro"
    gas = ["CH4","C2H6"]; % gases measured
    molar_mass = [16040,30070]; % molar mass of gas
    mm_units = ["mg","mg"]; % units for molar mass
    
    ConcentrationData.(gas(1)) = ConcentrationData.CH4_dry;
    ConcentrationData.(gas(2)) = ConcentrationData.C2H6_dry;
    % ConcentrationData.(gas(3)) = ConcentrationData._____;
    ConcentrationData.CavityPressure_torr = ConcentrationData.CavityPressure;
    ConcentrationData.H2O_perc = ConcentrationData.H2O;
    
    % add a column of smoothed data to because the instrument measures at a higher frequency than the others.
    needs_smoothing = "Y";
    sm_interval = 5; % smoothing averaging window in seconds
elseif instrument == "LGR"
    gas = ["CH4","CO2"]; % gases measured
    molar_mass = [16040,44.010]; % molar mass of gas
    mm_units = ["mg","g"]; % units for molar mass

    ConcentrationData.(gas(1)) = ConcentrationData.CH4_ppm;
    ConcentrationData.(gas(2)) = ConcentrationData.CO2_ppm;
    % ConcentrationData.(gas(3)) = ConcentrationData._____;
    ConcentrationData.CavityPressure_torr = ConcentrationData.GasP_torr;
    ConcentrationData.H2O_perc = ConcentrationData.H2O_ppm./10000; % convert to percent.

    needs_smoothing = "N";
else
    error("Instrument not recognized. Options: LGR or Picarro.")
end

% load fieldnotes and select relevant rows from field notes table
FieldNotes_all = readtable("FieldNotes.xlsx","Sheet",instrument);
FieldNotes_date = FieldNotes_all(string(FieldNotes_all.StartTime,"dd/MM/yyyy") == date,:); % Creates a selection of the table by date using this new index.
if isempty(sample) ~= 1 % specific samples
    [jnk,J] = intersect(FieldNotes_date.SampleNumber,sample);FieldNotes = FieldNotes_date(J,:);
elseif isempty(sample_range) ~= 1 % range of samples
    M = [find(FieldNotes_date.SampleNumber==sample_range(1))];N = [find(FieldNotes_date.SampleNumber==sample_range(2))];
    J=M:N;FieldNotes = FieldNotes_date(J,:);
else % all samples processed
    FieldNotes = FieldNotes_date;
end

depth_in_soil = FieldNotes.ChamberDepth_cm.';seal_type = string(FieldNotes.SealType.'); % Make variables from field notes table
deptherror = 0.2;

% load chamber dimensions
vol_area = readtable("Chamber_Area_Volume.xlsx","Sheet",chamber_type);
    full_volume = vol_area.FullChamberVolume(1); full_volume_unc = vol_area.FullChamberVolumeUncetainty(1); chamber_diameter = vol_area.ChamberDiameter(1);
    full_chamber_area = pi*(chamber_diameter/(100*2))^(2); % assumes circular geometry for the chamber

chamber_diameter_unc = 0.2; % in cm
temp_unc = 4; % in Kelvin 
pressure_unc = 10; % in hPa

% Additional time window (in minutes):
plot_win_min = 2; % determines the number of minutes before and after start and stop times to plot the data

filtering_limit_s = 30; % number of seconds of data that is allowed to be filtered out

filtering = strings([length(FieldNotes.SampleNumber),1]);
for i = 1:length(FieldNotes.SampleNumber)
    %% select start and stop times for measurements
    plot_window = timerange((FieldNotes.StartTime(i)-minutes(plot_win_min)),(FieldNotes.StopTime(i)+minutes(plot_win_min))); % Creates a time range between new Start and new Stop.
    original_window = timerange(FieldNotes.StartTime(i),FieldNotes.StopTime(i));

    % Creates figure of Time vs. Gas 1 to select start and stop times
    clf(figure(1));f1 = figure(1);set(f1,'color','w');subplot(211);hold on
    data_plot_win = plot(ConcentrationData.TIME(plot_window,:),ConcentrationData.(gas(1))(plot_window,:)); % Plots time (min) vs gas 1 for broader time window.
    data_original = plot(ConcentrationData.TIME(original_window,:),ConcentrationData.(gas(1))(original_window,:),"LineWidth",1); % Plots time (min) vs gas 1 for Field Notes time window on top of broader window.
    xlabel('Time');ylabel(sprintf(gas(1)+" (ppm)"));xticklabels(datestr(xticks()',"HH:MM:ss"));title(sprintf('%s',string(FieldNotes.SiteName(i)))+", "+(sprintf('Sample: %s',string(FieldNotes.SampleNumber(i)))));
    subplot(212);hold on;plot(ConcentrationData.TIME(plot_window,:),ConcentrationData.CavityPressure_torr(plot_window,:),"LineWidth",1); plot(ConcentrationData.TIME(original_window,:),ConcentrationData.CavityPressure_torr(original_window,:),"LineWidth",1);
    xlabel('Time');ylabel('Cavity Pressure (torr)');xticklabels(datestr(xticks()',"HH:MM:ss"));
    
    fprintf("----Sample: %s -----------\n",string(FieldNotes.SampleNumber(i))); % Displays Sample Number in Comand Window
    fprintf("Field Description:\n"+string(FieldNotes.FieldNotes_Description(i))+"\n")
    
    % interactively select measurement time window
    [Start_Selected(i),Stop_Selected(i)] = select_measurement_window(ConcentrationData(plot_window,:),f1);
    selected_window = timerange(Start_Selected(i),Stop_Selected(i)); % Creates range of time within selected time window.

    %% calculate volume and area of chamber
    % if FieldNotes.DepthExtraError(i) == 1;deptherror = .5;else;deptherror = .2;end

    if seal_type(i) == 'NA'
        if depth_in_soil(i) == 0;depth_in_soil(i) = 0.01;else;end % Cannot divide by zero in the subsequent equations.
        chamber_area = full_chamber_area; % in m^2; 
        area_unc = chamber_area.*sqrt(2*(chamber_diameter_unc/chamber_diameter).^2); % in m^2
        vol_offset = chamber_area*depth_in_soil(i)/100; % reduction in chamber volume due to chamber being pushed into soil, in m^3; 
        vol_offset_unc = vol_offset*sqrt((deptherror/depth_in_soil(i))^2 + (area_unc/chamber_area)^2); % assuming 2 mm error on depth in soil; this is an array
    else
        vol_offset = (-1)*vol_area.VolumeAddedBySeal(find(vol_area.SealType == seal_type(i))); %Enter volume added by the sandbag here in m^3 % reduction in chamber volume due to chamber being pushed into soil, in m^3; this is negative here because volume is being added
        vol_offset_unc = vol_area.VolumeAddedBySealUncertainty(find(vol_area.SealType == seal_type(i))); %Enter uncertainty surrounding volume here in m^3
        area_offset = vol_area.AreaOffsetFromSeal(find(vol_area.SealType == seal_type(i))); %Enter area of sealant inside the chamber here in m^2: make it positive if it adds area to the footprint of the measurement
        area_offset_unc = area_offset/2;
        chamber_area = full_chamber_area + area_offset/2; % in m^2; see notes above on why area_offset is divided by 2
        area_unc = sqrt((chamber_area.*sqrt(2*(chamber_diameter_unc/chamber_diameter).^2)).^2 + area_offset_unc.^2); % in m^2
    end
    chamber_volume = full_volume - vol_offset; % in m^3
    vol_unc = sqrt(full_volume_unc.^2 + vol_offset_unc.^2); % in m^3

    %% Filter Concentration Data
    % Does this sample's data need to be filtered?
    filtering(i) = input("Needs Filtering? \n(Enter to skip and not filter.)\n(Y/N): ", "s"); % Asks for user to decide whether filtering is needed (Y or N).
    
    % no LGR data was filtered
    if instrument == "LGR"; filtering(i) = 'N';sprintf("No LGR data was filtered");else;end

    % filter data
    clear Data
    [Data,win(i),Filter_Times1(i),Percent_Filtered(i)] = filter_data(ConcentrationData(selected_window,:),filtering(i),filtering_limit_s); % if filtering ~= 'Y', Data is just the data input into the function
    
    %% Flux calculation
    fprintf('\nSample Number %s',FieldNotes.SampleNumber{i})

    % instrument specific adjustments
    if instrument == "Picarro"
        if mean(Data.C2H6_dry) == 0;Sample_type(i) = "CH4Only";else;Sample_type(i) = "C2H6";end % record if sample has an ethane measurement (Some do not)
        % Pull out the coordinates of the sample from the end of the measurement because the more time the GPS has to pin, the more accurate the coordinate is. 
        % Print the lat and long to increase decimal places to increase accuracy.Only the Picarro has a built-in GPS unit.
        lat_string(i) = sprintf("%.9f",Data.GPS_ABS_LAT(end));long_string(i) = sprintf("%.9f",Data.GPS_ABS_LONG(end));
    elseif instrument == "LGR"
        Sample_type(i) = "CO2";
        lat_string(i) = "";long_string(i) = "";
    end
    
    % plot concentration and calculate flux
    clf(figure(3));f3 = figure(3); set(f3,'color','w');str_g = "";
    for g = 1:length(gas) % for each gas
        if needs_smoothing == "N"
            % linearity of concentration trend (R^2 parameter)
            r2(i,g) = (corr2(minutes(Data.TIME),Data.(gas(g))))^2;

            % calculate slope using "multislope method"
            [derivative(i,g),derivative_unc(i,g),min_slope(i,g),max_slope(i,g)] = multislope_calculation(Data.TIME,Data.(gas(g))); % methane
        elseif needs_smoothing == "Y"
            % add a column of smoothed data because the instrument measures at a higher frequency than the others.
            Data.(gas(g)+"_sm") = smoothran((Data.FRAC_DAYS_SINCE_JAN1-Data.FRAC_DAYS_SINCE_JAN1(1))*(24*60),Data.(gas(g)),sm_interval);K = find(~isnan(Data.(gas(g)+"_sm")));Data = Data(K,:); 
            
            % linearity of concentration trend using smoothed data (R^2 parameter)
            r2(i,g) = (corr2(minutes(Data.TIME),Data.(gas(g)+"_sm")))^2;
        
            % calculate slope using "multislope method"
            [derivative(i,g),derivative_unc(i,g),min_slope(i,g),max_slope(i,g)] = multislope_calculation(Data.TIME,Data.(gas(g)+"_sm")); % methane
        end

        % statistics on concentration during measurement
        mu_concentration(i,g) = mean(Data.(gas(g)));
        sigma_concentration(i,g) = std(Data.(gas(g)));

        % calculate flux
        flux(i,g) = derivative(i,g)*(1/((10^(6)*((1/(1440*molar_mass(g)))*chamber_area)*(1013/FieldNotes.AirPressure_hPa(i))*(22.4/273.15)*(FieldNotes.AirTemp_C(i)+273.15))./(1000*chamber_volume)));
        flux_unc(i,g) = abs(flux(i,g)*sqrt(((temp_unc/(FieldNotes.AirTemp_C(i)+273.15))^(2))+((vol_unc/chamber_volume)^(2))+((area_unc/chamber_area)^(2)) + (pressure_unc/FieldNotes.AirPressure_hPa(i))^2 +((derivative_unc(i,g)/derivative(i,g))^(2))));

        str_g =str_g+(sprintf("%s Flux:\n  %4.3f %s%s/m^2d \n",gas(g),flux(i,g),mm_units(g),gas(g)))+(sprintf("%s Flux Uncertainty:\n  %4.3f %s%s/m^2d \n\n",gas(g),flux_unc(i,g),mm_units(g),gas(g)));

        % plot concentration over time
        figure(3);subplot(length(gas),4,([1,3]+((g-1)*4)));hold on;plot(Data.TIME,Data.(gas(g)));
        if needs_smoothing=="Y";plot(Data.TIME,Data.(gas(g)+"_sm"),'linewidth',1);else;end
        xlabel('Time (min)');ylabel(gas(g)+" (ppm)");xlim([min(Data.TIME) max(Data.TIME)]);
    end
    if Sample_type(i) == "CH4Only";r2(i,2) = NaN;derivative(i,2)=NaN;derivative_unc(i,2)=NaN;min_slope(i,2)=NaN;max_slope(i,2)=NaN;mu_concentration(i,2)=NaN;sigma_concentration(i,2)=NaN;flux(i,2)=NaN;flux_unc(i,2)=NaN;else;end

    % Create text box on figure with fluxes
    str = (sprintf("%s \nDate: %s \nSample: %s\n\n",FieldNotes.SiteName{i},date,FieldNotes.SampleNumber{i}))+str_g;annotation('textbox',[.725 .2 .5 .5],'String',str,'FitBoxToText','on','FontSize',8,'BackgroundColor',[0.6 0.8 0.9]);  
    figure_name = "Figures/"+string(datetime(date),"yyyy_MM_dd")+"_Sample_"+FieldNotes.SampleNumber{i}+'_CH4'+'.fig';saveas(figure(3),figure_name); % save figure
   
    % Asks for user to decide whether flux looks good and saves figures of any fluxes that do not make sense for later inspection.
    reinspect_flux = input("Does this flux make sense? \n'N' to save the figure, enter to skip. \n(Y/N): ", "s");disp(fprintf('\n---------------------'));
    if reinspect_flux=='N';figure_name = "ProblemSamples_Figures/"+string(datetime(date),"yyyy_MM_dd")+"_"+instrument+"_Sample_"+FieldNotes.SampleNumber{i}+'.fig';saveas(figure(3),figure_name);else;end
        
    sitename(i) = string(FieldNotes.SiteName{i});sample_number(i) = string(FieldNotes.SampleNumber{i});date1(i) = date;
end

%% Create Table of Results
Header_gas = [];
Data_output_gas = [];
for g = 1:length(gas)
    Data_output_gas = [Data_output_gas,flux(:,g),flux_unc(:,g),mu_concentration(:,g),sigma_concentration(:,g),derivative(:,g),derivative_unc(:,g),min_slope(:,g),max_slope(:,g),r2(:,g)];
    Header_gas = [Header_gas,sprintf(gas(g) + " Flux "+mm_units(g)+gas(g)+"/m2/day"),sprintf(gas(g) + " Flux Error "+mm_units(g)+gas(g)+"/m2/day"),sprintf(gas(g) + " Mean, ppm"),sprintf(gas(g) + " Stdev, ppm"),sprintf("d"+gas(g)+"/dt, ppm/min"),sprintf("d"+gas(g)+"/dt Error, ppm/min"), sprintf("Minimum "+gas(g)+" Slope"),sprintf("Maximum "+gas(g)+" Slope"),sprintf("R2, "+gas(g))];
end
Data_output = [date1', sitename', sample_number', Sample_type', string(FieldNotes.FieldNotes_Description),Data_output_gas,...
    datestr(Start_Selected'),datestr(Stop_Selected'),FieldNotes.AirTemp_C,FieldNotes.ChamberDepth_cm,FieldNotes.AirPressure_hPa,string(FieldNotes.SealType),...
    string(FieldNotes.SoilMoisture_VWC),string(FieldNotes.SoilTemp_C),string(FieldNotes.WindSpeed_ms),lat_string',long_string',...
    string(filtering),Percent_Filtered',Filter_Times1',win',repmat(instrument',length(sample_number),1)]; 
Headers = ["Date","Site Name","Sample Number","Sample Type","Description",Header_gas...
    "Start Time, Data File","Stop Time, Data File","Chamber Temperature, Degrees C","Chamber Depth In Soil, cm","Air Pressure, hPa","Seal Type",...
    "Soil Moisture, % VWC","Soil Temperature, Degrees C","Wind Speed, m/s","Latitude, DD (Analyzer)","Longitude, DD (Analyzer)","Was filtering applied?","Percent Filtered (%)","Times Filtered Out","Filtering Window","Analyzer"];
Results_Table = vertcat(Headers,Data_output);open Results_Table

%% Mapmaking
CH4flux = flux(:,1);
varTypes = ["double","double","double"];varNames = ["CH4 flux, mg CH4/m^2/day", "Latitude (DD)", "Longitude (DD)"];CH_4flux = cell(length(CH4flux(:,1)));CH_4flux = round(CH4flux(:,1));
lat = str2double(lat_string'); long = str2double(long_string');
% plot a map of the samples
f4 = figure(4);set(f4,'color','w');
gb = geobubble(table(CH_4flux,lat,long),"lat","long","SizeVariable","CH_4flux");
min_lat = min(lat)-0.006; max_lat = max(lat)+0.006;min_long = min(long)-0.006; max_long = max(long)+0.006;geolimits ([min_lat,max_lat],[min_long,max_long])
geobasemap streets-light

