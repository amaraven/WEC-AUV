% Ama Hartman

% INPUTS: 
% * Wave resource or power generation data. Possible options:  
%   1 - Sea state selection for built-in data modeling power generation of 
%   a four-meter diameter point absorber WEC at the PacWave South test site 
%   2 - Time series of WEC power generation and time. 
%   3 - Mean value of WEC power generation
%   4 - Time series of significant wave height, wave energy period, peak
%   period and time
%   5 - Mean values of significant wave height, wave energy period, and
%   peak period
%   6 - *.mat data file must contain power matrix and time series vectors
%   for significant wave height, wave energy period, and time. First column 
%   of matrix must include significant wave height values in order from low 
%   to high, and first row of matrix must include wave energy period values 
%   in order from low to high. 
% * Choice of comparison variable. Possible options: 
%   = AUV Model: compares system performance given different AUV models.
%   Model iterates through a list of pre-set AUVs characterized according
%   to market advertizements
%   - WEC Power Generation / Wave Resource: compares system performance
%   given different power generation. Iterates through a list of resource
%   data viles or values
% * WEC properties: 
%   - Characteristic dimension: Largest cross-sectional area of the WEC.
%   Default is 1.5 meters
%   - Hotel load: Constant power-draw by the WEC. Default is 50 W
%   - Battery capacity: Default is 500 Wh
% * Central energy unit hotel load: Default is 10 W
% * AUV dock hotel load: Default is 20 W
%
% OUTPUTS: 
% * Fleet number: Number of AUVs that can be supported by the WEC
% * Central battery size: size of the central battery required 
% * AUV mission schedule: Example schedule of AUV missions and recharge
% periods
% * Component battery level tracker: battery level history for each
%   component (WEC, central battery, and AUV(s)) during simulation time
% * Global comparison variable plots: compares system performance between
%   comparison variable cases


%% Preliminaries
% profile on;
clear; % close ; clc
warning('on')
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% Add component class folder to path
addpath('Component Classes');
addpath('Functions');
addpath('InputData');

auvModels = [{'Iver3-27'}, {'Iver3-38.5'}, {'REMUS 100'}, {'REMUS 300-58.5'}, {'REMUS 300-70.3'}, {'REMUS 620-210'}, {'REMUS 620-279'}, {'REMUS 620-347'}, {'REMUS 6000'}, {'Bluefin-9'}, {'Bluefin-12'}, {'Bluefin-21'}, {'Bluefin-HAUV'}, {'Hugin Superior'}, {'Hugin Endurance'}, {'Hugin 3000'}, {'Hugin 4500'}, {'Boxfish AUV'}, {'Boxfish ARV-i'}, {'Saab Sabertooth Single'}, {'Saab Sabertooth Double'}];

%% User-Defined Inputs / Model Setup

userPrompts = 0;  % 1 - prompts user for inputs, 0 - does not prompt, uses set values

if userPrompts == 0

    % -- Inputs ----------------------------------------------------------- 

    simHrs = 24*14;
    depVar =  'AUV Model';  % 'WEC Power Gen / Wave Resource'; % 'Battery Specs';
    resourceDataType = 4;  % 1) RM3 data, 2) power time series 3) mean power (currently just draws from RM3 data), 4) time series of wave specs, 5) mean wave spec values, 6) Power matrix & Hs, Te time series
    incorpStagger = 1;  % (1) Incorporate stagger into AUV mission scheduling
    maxFleetSize = 0;  % Maximum size for AUV fleet, set to 0 for no maximum

    userDefinedBattery = 0;  % 1 - User inputs total battery capacity, 0 - model outputs battery capacity to support max auv's
    if userDefinedBattery == 1
        energyStorageSize = 4200;
    end

    switch depVar       
        case 'AUV Model'  
            if resourceDataType == 1  % Pick one sea state, try each auv model at that seastate                                                                                             
                seaState = 3; 
            end  
       
            wec = WEC('generic'); 

        case 'WEC Power Gen / Wave Resource'  % Pick one auv model, run through all sea states
            auvModelNum = 18; 
            auv = AUV(string(auvModels(auvModelNum))); 
    end

    % ---------------------------------------------------------------------

else

    % Simulation specs
    simHrs = input('Enter the desired simulation time in hours: ');

    
    % Incorporate AUV deployment stagger?
    staggerAns = questdlg('Would you like to incorporate a stagger into the AUV mission scheduling?','Mission Stagger','Yes','No','Yes');
    switch staggerAns
        case 'Yes'
            incorpStagger = 1;
        case 'No'
            incorpStagger = 0;
    end

    
    % Central Battery Size
    energyStorageSize = input('Enter the desired central battery size [Wh], or leave blank to use an ideal battery size based on the wave resource: ');
    
    if isempty(energyStorageSize) || energyStorageSize == 0
        userDefinedBattery = 0; 
    
    else
        userDefinedBattery = 1;

    end

    
    % Choose dependant variable
    depVar = questdlg('Choose variable to evaluate and change betweeen simulations.','Dependant Variable', ...
        'AUV Model','WEC Power Gen / Wave Resource', 'AUV Model'); 
    
    switch depVar
        case 'AUV Model'
            % Generate WEC & AUV Objects
            wec = WEC('generic');
    
        case 'WEC Power Gen / Wave Resource'
            auvModelNum = listdlg('PromptString',{'Choose AUV model.'}, 'SelectionMode', 'single','ListString', {'Iver3-27','Iver3-38.5', 'REMUS 100', 'REMUS 300-58.5', 'REMUS 300-70.3', 'REMUS 620-210', 'REMUS 620-279', 'REMUS 620-347', 'REMUS 6000', 'Bluefin-9', 'Bluefin-12', 'Bluefin-21', 'Bluefin-HAUV', 'Hugin Superior', 'Hugin Endurance', 'Hugin 3000', 'Hugin 4500', 'Boxfish AUV', 'Boxfish ARV-i'});
            auv = AUV(auvModels{auvModelNum}); 

    end
    
    
    % Specify & load wave resource data type
    resourceDataType = listdlg('PromptString',{'Specify the type of wave resource data you plan to use'}, 'SelectionMode','single', ...
        'ListString', {'Proteus struct containing power generated during different sea states','Time series of power generated','Value of mean power','Time series of wave specs (Hs, Te, Tp)','Mean Wave specs (Hs, Te, Tp)'});

end

% Simulation parameters from user-inputs 
simSec = simHrs * 60 * 60; 
dtSec = 30.0;  % [s]
dt = dtSec/60/60;  % [hr]
simTime = (dtSec:dtSec:simSec)' / 60 / 60;  % [hr]


% Load Resource Data

switch resourceDataType
    case 1  % 'Proteus struct containing power generated during different sea states'
        load('inputData/data.mat', 'RM3'); % Time series output of proteus WEC model for a single point absorber with no dock attached. (Other variables include F3B 'floating third body', WECM 'wec mounted', and SBM 'sea bottom mounted', all referring to auv dock mount location.)

        if strcmp(depVar, 'AUV Model') && userPrompts == 1
            % Prompt for sea state
            seaState = listdlg('PromptString',{'Specify the sea state you would like to use for this simulation'}, 'SelectionMode','single', ...
                'ListString',{'1','2','3','4','5','6','7','8','9','10'});
        end
        
    case 2  % 'Time series of power generated'
        switch depVar
            case 'AUV Model'
                dataFiles = uigetfile('*.mat', 'Select *.mat data file for time series','MultiSelect','off');
                vars = whos('-file', dataFiles);
    
            case 'WEC Power Gen / Wave Resource' 
                dataFiles = uigetfile('*.mat', 'Select *.mat data files','MultiSelect','on');
                if ischar(dataFiles)  % if only one is selected
                    dataFiles = {dataFiles};
                end
            
                vars = whos('-file', dataFiles{1});

        end
        
        tIndx = listdlg('PromptString',{'Select the time series vector'}, 'SelectionMode','single', ...
                'ListString',{vars.name});
        pwrIndx = listdlg('PromptString',{'Select the power gen. time series'}, 'SelectionMode','single', ...
                'ListString',{vars.name});

    case 3  % 'Value of mean power'
        % Generates generic mean power using example data:
        %{
        load('inputData/data.mat', 'RM3'); % Time series output of proteus WEC model for a single point absorber with no dock attached. (Other variables include F3B 'floating third body', WECM 'wec mounted', and SBM 'sea bottom mounted', all referring to auv dock mount location.)

        % ~< input dialog for mean power >~ but for now... 

        switch depVar
            case 'AUV Model' 
                meanPower = RM3(3).Pmean;

            case 'WEC Power Gen / Wave Resource'
                meanPower = [RM3.Pmean];  
        end
        %}
        meanPower = input('Enter mean power generated as a single value or a bracketed, comma-seperated list: ');

    case 4  % Time series of wave specs (Hs, Te, Tp)
        switch depVar
            case 'AUV Model'
                % load('inputData/Hawaii_Wave_Hindcast_Data_NREL.csv');  % generic example
                dataFile = uigetfile('*.csv', 'Select *.csv Data File', 'MultiSelect','off');
                dataTable = readtable(dataFile,'VariableNamingRule','modify');

            case 'WEC Power Gen / Wave Resource'
                dataFiles = uigetfile('*.csv', 'Select *.csv Data Files', 'MultiSelect','on');
                if ischar(dataFiles)
                    dataFiles = {dataFiles};
                end
        end

    case 5  % Mean Wave specs (Hs, Te, Tp)
        if userPrompts == 0
            sigWaveHeight = [1.8774    1.8800    1.9000    2.2000    1.6000]';
            energyPeriod = [8.9367    8.9400    9.0000    9.0000    8.8000];
            peakPeriod = [13.1428   13.1400   13.1400   13.1400   10.1000];

        else
            prompt = {'Enter significant wave height(s) [m]:', 'Enter wave energy period(s) [s]:','Enter peak period(s) [s]:'};
            fieldsize = [1, 45; 1, 45; 1, 45];
            definput = {'1.8774','8.9367','13.1428'};
    
            tryAgain = 1;  % preallocate
            while tryAgain == 1
                waveSpecDlg = inputdlg(prompt, 'Wave Resource Specifications', fieldsize, definput);
                sigWaveHeight = str2num(waveSpecDlg{1});
                energyPeriod = str2num(waveSpecDlg{2});
                peakPeriod = str2num(waveSpecDlg{3});
    
                if isequal( size(sigWaveHeight), size(energyPeriod), size(peakPeriod) )
                    tryAgain = 0;
                else
                    warning('Each specification must have the same number of values. Please re-enter the wave resource specifications.')
                end
    
            end
        end

        waveSpecTable = table(sigWaveHeight', energyPeriod', peakPeriod','VariableNames',{'SignificantWaveHeight','EnergyPeriod','PeakPeriod'});

    case 6  % Power matrix & Hs, Te time series
        % User-Selection of file & variable names:
        switch depVar
            case 'AUV Model'
                dataFiles = uigetfile('*.mat', 'Select *.mat data file for time series','MultiSelect','off');
                vars = whos('-file', dataFiles);

            case 'WEC Power Gen / Wave Resource' 
                dataFiles = uigetfile('*.mat', 'Select *.mat data files','MultiSelect','on');
                if ischar(dataFiles)  % if only one is selected
                    dataFiles = {dataFiles};
                end
                
                vars = whos('-file', dataFiles{1});

        end
        pMatIndx = listdlg('PromptString',{'Select the power matrix variable'}, 'SelectionMode','single', ...
                'ListString',{vars.name});
        tIndx = listdlg('PromptString',{'Select the time series vector'}, 'SelectionMode','single', ...
                'ListString',{vars.name});
        HsIndx = listdlg('PromptString',{'Select the significant wave height time series'}, 'SelectionMode','single', ...
                'ListString',{vars.name});
        TeIndx = listdlg('PromptString',{'Select the wave energy period time series'}, 'SelectionMode','single', ...
                'ListString',{vars.name});

end  % load resource data


%% Loop Preallocation, Constant Calcs, & Output Object Generation

% Generate output object
modOut = ModelOutput(simTime, depVar, resourceDataType);

modOut.incorpStagger = incorpStagger;
modOut.maxFleetSize = maxFleetSize;

switch depVar
    case 'AUV Model'
        
        % Dependent variable loop spec.
        loopLength = length(auvModels); 

        % Calculate power generation...
        switch resourceDataType
            case 1
                wec.reshapePowerGen(RM3(seaState).Power, RM3(seaState).Time, simHrs, dt);  

                modOut.meanPowerGen = wec.meanPowerGen;
                modOut.seaState = seaState;

            case 2
                % Load data into workspace:
                dataTime = load(dataFiles, vars(tIndx).name);
                pGen_dataTime = load(dataFiles, vars(pwrIndx).name);

                % Reshape power gen. to fit simulation timestep
                wec.reshapePowerGen(pGen_dataTime, dataTime, simHrs, dt);
                modOut.meanPowerGen = wec.meanPowerGen; 
                modOut.dataIn = dataFiles;

            case 3
                modOut.meanPowerGen = meanPower;
    
                wec.meanPowerGen = modOut.meanPowerGen;
                wec.lowPowerGen = 0.75 * wec.meanPowerGen; 
                wec.powerGenMeans = ones(size(simTime)) * wec.meanPowerGen; 

            case 4
                wec.calcPowerGen(dataTable,'meanPwr', simTime, 1, 403); %%%%%%%%%%%%%%%%%%%%%%%%%%% reset windowOverrideIndx to 0 to pull mean power from time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                modOut.dataIn = dataFile; 

            case 5  
                wec.calcPowerGen(waveSpecTable(1,:), [], simTime, 0, 0);  % If multiple sets of wave specifications is given, runs simulation with first set only
                wec.lowPowerGen = 0.75*wec.meanPowerGen; 
                wec.powerGenMeans = ones(size(simTime)) * wec.meanPowerGen;
                
                modOut.meanPowerGen = wec.meanPowerGen;
                modOut.dataIn = waveSpecTable;
            
            case 6
                % Load data into workspace:
                pMat = load(dataFiles, vars(pMatIndx).name);  pMat = pMat.pMat;
                dataTime = load(dataFiles, vars(tIndx).name);  dataTime = dataTime.dataTime;
                Hs = load(dataFiles, vars(HsIndx).name);  Hs = Hs.Hs;
                Te = load(dataFiles, vars(TeIndx).name);  Te = Te.Te;
        
                % Calculate power from grid interpolation
                pMatTableFn = griddedInterpolant( repmat(pMat(2:end,1), 1, size(pMat, 2)-1), repmat(pMat(1,2:end), size(pMat, 1)-1, 1), pMat(2:end, 2:end),'linear','nearest');  % X1 grid, X2 grid, valueMatrix, interpolation method, extrapolation method
                pGen_dataTime = pMatTableFn(Hs, Te); 

                % reshape power gen
                wec.reshapePowerGen(pGen_dataTime, dataTime, simHrs, dt);
                modOut.meanPowerGen = wec.meanPowerGen;
                modOut.dataIn = dataFiles;

        end

    case 'WEC Power Gen / Wave Resource'
        modOut.auvModels = auvModels{auvModelNum}; 

        switch resourceDataType
            case 1
                % Dependent variable loop spec.
                loopLength = length(RM3);

            case 2
                loopLength = length(dataFiles);
                modOut.dataIn = dataFiles;
                
            case 3
                loopLength = length(meanPower);
                modOut.meanPowerGen = meanPower;

            case 4
                loopLength = length(dataFiles);

                modOut.dataIn = dataFiles; 

            case 5
                loopLength = length(waveSpecTable.SignificantWaveHeight);
                
                % Calc power gen. for all  cases simultaneously
                wec = WEC('generic');
                wec.calcPowerGen(waveSpecTable, [], simTime, 0);  

                modOut.meanPowerGen = wec.meanPowerGen;
                modOut.dataIn = waveSpecTable;

            case 6
                loopLength = length(dataFiles);
                modOut.dataIn = dataFiles;
        end

end

%% Generate central energy storage object

if userDefinedBattery == 1
    energyStorage = EnergyStorage(energyStorageSize, 10, 20); 

elseif userDefinedBattery == 0
    energyStorage = EnergyStorage([], 10, 20);  % Empty battery storage will be rewritten after fleet number determination...

else
    error('"userDefinedBattery" out of expected range')
end

modOut.userDefinedBattery = userDefinedBattery; 


%% Dependent Variable Loop / Run Simulation

for depVarCount = 1:loopLength

    %% Prep variables (power gen & auv model(s))
    switch depVar
        case 'AUV Model'

            % Clear Variables
            clearvars -except simTime dt depVarCount resourceDataType modOut auvModels depVar wec energyStorage incorpStagger maxFleetSize

            % Generate auv object
            auv = AUV(auvModels{depVarCount}); 


        case 'WEC Power Gen / Wave Resource'
            % Generate WEC & AUV Objects
            wec = WEC('generic');

            switch resourceDataType
                case 1
                    % Save Sea state
                    modOut.seaState(depVarCount) = depVarCount;

                    % Clear Variables
                    clearvars -except simTime dt depVarCount resourceDataType modOut auvModels depVar RM3 auvModelNum auv wec energyStorage incorpStagger maxFleetSize
        
                    % Calc. power gen
                    wec.reshapePowerGen(RM3(depVarCount).Power, RM3(depVarCount).Time, simTime(end), dt);  % Output is wec.powerGenMeans: Timeseries of power generation [W] corresponding to simulation time [hr].

                case 2
                    % Clear Variables
                    clearvars -except simTime dt depVarCount resourceDataType modOut auvModels depVar auvModelNum auv wec energyStorage dataFiles vars tIndx pwrIndx incorpStagger maxFleetSize

                    % Load power data into workspace: 
                    dataTime = load(dataFiles{depVarCount}, vars(tIndx).name);
                    pGen_dataTime = load(dataFiles{depVarCount}, vars(pwrIndx).name);

                    % Reshape power gen
                    wec.reshapePowerGen(pGen_dataTime, dataTime, simHrs, dt);
                    modOut.meanPowerGen = wec.meanPowerGen; 
                    modOut.dataIn = dataFiles; 

                case 3
                    % Clear Variables
                    clearvars -except simTime dt depVarCount resourceDataType modOut auvModels depVar meanPower auvModelNum auv wec energyStorage incorpStagger maxFleetSize

                    % Save power generated (user input)
                    wec.meanPowerGen = modOut.meanPowerGen(depVarCount);
                    wec.lowPowerGen = 0.75 * wec.meanPowerGen; 
                    wec.powerGenMeans = ones(size(simTime)) * wec.meanPowerGen; 

                case 4
                    % Clear Variables
                    clearvars -except simTime dt depVarCount resourceDataType modOut auvModels depVar auvModelNum auv dataFiles wec energyStorage incorpStagger maxFleetSize

                    % Load data & calculate power generation
                    dataTable = readtable(dataFiles{depVarCount},'VariableNamingRule','modify');
                       
                    % Calculate power generation
                    wec.calcPowerGen(dataTable,'meanPwr', simTime, 0);

                case 5
                    % Clear Variables
                    clearvars -except simTime dt depVarCount resourceDataType modOut auvModels depVar auvModelNum auv waveSpecTable wec energyStorage incorpStagger maxFleetSize

                    % Save power gen (already calculated from user-inputs)
                    wec.meanPowerGen = modOut.meanPowerGen(depVarCount); 
                    wec.lowPowerGen = 0.75 * wec.meanPowerGen;
                    wec.powerGenMeans = ones(size(simTime)) * wec.meanPowerGen; 

                case 6
                    % Clear Vars
                    clearvars -except simTime dt depVarCount resourceDataType modOut auvModels depVar auvModelNum auv wec energyStorage dataFiles vars pMatIndx tIndx HsIndx TeIndx incorpStagger maxFleetSize

                    % Load power data into workspace:
                    pMat = load(dataFiles{depVarCount}, vars(pMatIndx).name);
                    dataTime = load(dataFiles{depVarCount}, vars(tIndx).name);
                    Hs = load(dataFiles{depVarCount}, vars(HsIndx).name);
                    Te = load(dataFiles{depVarCount}, vars(TeIndx).name);
            
                    % Calculate power from grid interpolation
                    pMatTableFn = griddedInterpolant( repmat(pMat(2:end,1), 1, size(pMat, 2)-1), repmat(pMat(1,2:end), size(pMat, 1)-1, 1), pMat(2:end, 2:end),'linear','nearest');  % X1 grid, X2 grid, valueMatrix, interpolation method, extrapolation method
                    pGen_dataTime = pMatTableFn(Hs, Te); 

                    % Reshape power gen
                    wec.reshapePowerGen(pGen_dataTime, dataTime, simHrs, dt);
                    modOut.meanPowerGen = wec.meanPowerGen;
                    modOut.dataIn = dataFiles;

            end
    end
    

    %% Fleet number calculation
    runFleetCalc = 1;
    fleetCalcCount = 0;
    while runFleetCalc == 1
        fleetCalcCount = fleetCalcCount + 1;

        if fleetCalcCount > 1
            modOut.fleetNumber(depVarCount) = modOut.fleetNumber(depVarCount) - 1;

            % clear values from prev simulation
            clear modOut.energyStorageBatteryLvl modOut.wecBatteryLvl modOut.auvBatteryLvl modOut.auvSchedule modOut.auvTimeOnMission

        else
            modOut.fleetNumber(depVarCount) = calcFleetNum(wec, auv, energyStorage, maxFleetSize);  % initial calculation

        end
        
        
        %% Update central storage hotel load to account for additional AUV docks
        energyStorage.hotelLoad = energyStorage.baseHotelLoad + energyStorage.dockHotelLoad*(modOut.fleetNumber(depVarCount)); 


        %% Calculate minimum central battery, for battery to not limit fleet
        % Need at least enough energy saved in central battery to fully
        % recharge AUV(s), and support AUV & WEC hotel loads during
        % that time given poor power generation +~ 10% for safety.
        
        wec.calcLowPower(resourceDataType, dt, auv); 

        lowPowerOverflow = wec.lowPowerGen -  wec.hotelLoad/(wec.n_battery^2);  % Assumes wec is already at max battery
        minBattery = ( auv.chargeLoad(auv.mission)*modOut.fleetNumber(depVarCount)/energyStorage.n_battery + auv.chargeTime(auv.mission)*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTransfer - lowPowerOverflow*auv.chargeTime(auv.mission)*energyStorage.n_battery*energyStorage.n_wecPwrTrnsfr  ) *1.05;  % Given lowest possible power generation during recharge OR 5% of WEC battery, if threshold is negative (i.e. if power gen > power draw)

        if minBattery < 0  % Somewhat arbituary minimum battery if pGen > operating loads
            warning('Problem with minimum battery calculation. Using default value for this simulation.')
            minBattery = 0.25 * ( auv.chargeLoad(auv.mission)*modOut.fleetNumber(depVarCount)/energyStorage.n_battery + auv.chargeTime(auv.mission)*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTransfer ); 
        end
 
        % Clear battery storage for new simulation if battery amount is not user-input
        if energyStorage.userDefinedBattery ~= 1
            energyStorage.maxBattery = minBattery;  % set central battery to minimum needed for AUV fleet

        elseif energyStorage.maxBattery < minBattery
            warning('Central energy storage amount is likely too low to support the fleet of %f %s AUV(s) determined by the given wave resource. Consider increasing central battery from %f to %f kWh.', modOut.fleetNumber, auv.model, energyStorage.maxBattery*10e-3, minBattery*10e-3)
        
        elseif minBattery < energyStorage.maxBattery 
            warning('Central energy storage is %f kWh larger than what is required for this configuration.', (energyStorage.maxBattery - minBattery)*10e-3)

        end


        % Build fleet & run model if possible
        if modOut.fleetNumber(depVarCount) > 0
    
            % Build fleet
            switch depVar
                case 'AUV Model'
                    auvFleet = arrayfun(@(x) AUV(auvModels{depVarCount}), 1:modOut.fleetNumber(depVarCount), 'UniformOutput', false);  
                case 'WEC Power Gen / Wave Resource'
                    auvFleet = arrayfun(@(x) AUV(auvModels{auvModelNum}), 1:modOut.fleetNumber(depVarCount), 'UniformOutput', false);
            end

            modOut.auvFleet{depVarCount} = auvFleet;  % save to model output object. 
            
            %% Run Simulation
            % Simulation Goals: 
            % * Calculate battery levels of central storage & AUV(s)
            % * Generate AUV mission schedule given the following rules: 
            %   - AUV must dock with at least 20% of its battery left
            %   - Minimize time when central energy storage & AUV are both at full battery
            %   - energyStorage must have enough battery to charge AUV(s) when they return
            [modOut.energyStorageBatteryLvl(:, depVarCount), modOut.wecBatteryLvl(:, depVarCount), modOut.auvBatteryLvl{depVarCount}, modOut.auvSchedule{depVarCount}] = runPowerModel(simTime, wec, auvFleet, energyStorage, incorpStagger);
           
            %% Post-Simulation Calcs

            % Time spent 'on mission'
            modOut.auvTimeOnMission{depVarCount} = sum((modOut.auvSchedule{depVarCount} == 1) * dt);
            
            % Time 'on-mission' within performance calculation domain in the case of staggered deployments (excluding time AUVs are artificially held at dock)
            if modOut.incorpStagger == 1 
                if modOut.fleetNumber(depVarCount) == 1
                    modOut.auvTimeOnMissionCorrected{depVarCount} = modOut.auvTimeOnMission{depVarCount};
                    
                else
                    staggerHours = (modOut.auvFleet{1,depVarCount}{1,1}.missionSpecs(2) + modOut.auvFleet{1,depVarCount}{1,1}.chargeTime) / modOut.fleetNumber(depVarCount);
                    staggerPreliminaryTime = (modOut.fleetNumber(depVarCount) - 1) * staggerHours;
                    [~, preDomainIndx] = min(abs(modOut.simTime - staggerPreliminaryTime));

                    modOut.auvTimeOnMissionCorrected{depVarCount} = sum((modOut.auvSchedule{depVarCount}(preDomainIndx+1:end, :) == 1) * dt);  % Exclude on-mission time before all AUVs are deployed
                end

            end

            %% Simulation Quality Checks

            % if last auv only went on one mission, and first auv went on more, OR central battery dropped below 0 (even with battery saver) auvFleet is too large by at least 1
            if ( (modOut.auvTimeOnMission{depVarCount}(1,end) - (auvFleet{end}.chargeTime(auvFleet{end}.mission) + auvFleet{end}.missionSpecs(auvFleet{end}.mission, 2)) ) <= 0 ) && ((modOut.auvTimeOnMission{depVarCount}(1,1) - (auvFleet{1}.chargeTime(auvFleet{1}.mission) + auvFleet{1}.missionSpecs(auvFleet{1}.mission, 2)) ) > 0 ) || any(modOut.energyStorageBatteryLvl(:,depVarCount) < 0)
                runFleetCalc = 1;  % re-run fleet number calculation. (Disable to plot results for systems with too-many AUVs)  %%%%%%%%%%%%%%%%%%%%%%%%% RESET TO 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
                warning('Sorry, our initial %s fleet number estimate of %f was too high. Re-running the simulation now...', auv.model, modOut.fleetNumber(depVarCount))
        
            else
                runFleetCalc = 0;

                modOut.centralBatteryCapacity{depVarCount} = [energyStorage.maxBattery];

    
            end 

        else  % Wave resource insufficient to support deployment of AUV here
            % warning('Wave resource for sea state %f insufficient to support the hotel load for any %s AUVs at this site.', depVar, auv.model)
            warning('Wave resource insufficient to support the hotel load for any %s AUVs at this site.', auv.model)
            
            modOut.auvTimeOnMission{depVarCount} = 0;  % No AUV's = no time spent on missions
            modOut.auvTimeOnMissionCorrected{depVarCount} = 0; 
            modOut.auvSchedule{depVarCount} = [];
            modOut.auvBatteryLvl{depVarCount} = [];
            modOut.energyStorageBatteryLvl(:, depVarCount) = zeros(size(simTime));
            modOut.wecBatteryLvl(:, depVarCount) = zeros(size(simTime)); 
            modOut.auvFleet{depVarCount} = [];

            runFleetCalc = 0;
        end

    end  % while fleetCalc == 1 

    %% Values to Track

    % (Energy used [Wh] / mission time [h]) Homogenous auv fleets only!
    modOut.ratePwrUsed(depVarCount) = (auv.missionSpecs(auv.mission, 3)*auv.maxBattery + auv.hotelLoad*auv.chargeTime) / (auv.missionSpecs(auv.mission, 2) + auv.chargeTime);  % AVG rate AUV uses power during a mission + recharge cycle. No efficiencies applied!  
    modOut.auvMissionLength(depVarCount) = auv.missionSpecs(auv.mission, 2);
    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;
    switch depVar
        case 'AUV Model'
            if depVarCount == 1
                modOut.powerGenMeans = wec.powerGenMeans;  % only save power gen once if it isn't recalculated each iteration
            end
        case 'WEC Power Gen / Wave Resource'
            modOut.powerGenMeans{depVarCount} = wec.powerGenMeans; 
    end

end  % depVar 

%% Figures
%{
modOut.plotBatteryTracker

modOut.plotTimeOnMission

modOut.plotPowerGen
%}

%% Save Output
% Edit and run the following to save data to 'outputData' folder: 
%{  
save('outputData/usWestCoast_03Dec25.mat','auv','energyStorage','modOut','wec')
%}

%% Simulation Meta
% profile off; profile viewer;
beep