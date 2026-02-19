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
% * Fleet size: Number of AUVs that can be supported by the WEC
% * Central battery size: size of the central battery required 
% * AUV mission schedule: Example schedule of AUV missions and recharge
% periods
% * Component battery level tracker: battery level history for each
%   component (WEC, central battery, and AUV(s)) during simulation time
% * Global comparison variable plots: compares system performance between
%   comparison variable cases


%% Preliminaries
tic; % profile on;
clear; % close ; clc
warning('on')
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% Add component class folder to path
addpath('Component Classes');
addpath('Functions');
addpath('InputData');

auvModels = [{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}, {'G'}, {'H'}, {'I'}, {'J'}, {'K'}, {'L'}, {'M'}, {'N'}, {'O'}, {'P'}, {'Q'}, {'R'}, {'S'}, {'T'}, {'U'}];

%% User-Defined Inputs / Model Setup

userPrompts = 0;  % 1 - prompts user for inputs, 0 - does not prompt, uses set values

if userPrompts == 0

    % -- Inputs ----------------------------------------------------------- 

    simHrs = 24*14;
    depVar =  'AUV Model';  % 'WEC Power Gen / Wave Resource'; % 'Battery Specs';
    resourceDataType = 4;  % 1) RM3 data, 2) power time series 3) mean power (currently just draws from RM3 data), 4) time series of wave specs, 5) mean wave spec values, 6) Power matrix & Hs, Te time series
    incorpStagger = 1;  % (1) Incorporate stagger into AUV mission scheduling
    maxFleetSize = 0;  % Maximum size for AUV fleet, set to 0 for no maximum

    userDefinedBattery = 0; %4200;  % 0 - model outputs battery capacity to support max auv's, ~=0 - is the user input central battery capacity

    switch depVar       
        case 'AUV Model'  
            if resourceDataType == 1  % Pick one sea state, try each auv model at that seastate                                                                                             
                seaState = 3; 
            end  
       
        case 'WEC Power Gen / Wave Resource'  % Pick one auv model, run through all sea states
            auvModelNum = 18; 
            auv_userIn = AUV(model=auvModels{auvModelNum}); 
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
    udfTryAgain = 1;
    while udfTryAgain == 1
        userDefinedBattery = input('Enter the desired central battery size [Wh], or leave blank to use an ideal battery size based on the wave resource: ');
        
        if isempty(userDefinedBattery)
            userDefinedBattery = 0; 
            udfTryAgain = 0;

        elseif ~isnumeric(userDefinedBattery)  % if user-input is NOT numeric..
            udfTryAgain = 1;
            warning('Central battery size must be empty or numeric.')
        else
            udfTryAgain = 0;
        end
    end

    
    % Choose dependant variable
    depVar = questdlg('Choose variable to evaluate and change betweeen simulations.','Dependant Variable', ...
        'AUV Model','WEC Power Gen / Wave Resource', 'AUV Model'); 
    
    if strcmp(depVar, 'WEC Power Gen / Wave Resource')
            auvModelNum = listdlg('PromptString',{'Choose AUV model.'}, 'SelectionMode', 'single','ListString', {'A','B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S','T','U','Custom'});
            if auvModelNum == 22
                prompt = {'Enter AUV mass [kg]:','Enter mission specs:','Enter battery capacity [Wh]:','Enter recharge rate [W]','Enter charge method (1 for wired OR 2 for wireless):','Enter hotel load [W] (leave empty to estimate from AUV mass):','Enter recharge threshold (0.XX):','Enter battery efficiency (0.XX)','Enter power transfer efficiency (0.XX) (leave empty to use default efficiencies based on power transfer method):'};
                fieldsize = [1 45; 1, 45; 1 45; 1 45; 1 45; 1 45; 1 45; 1 45; 1 45];
                defaults = {'100','[1, 10, 0.8]','1000','200','1','[]','0.20','0.9','[]'};
                auvSpecs = inputdlg(prompt, 'Enter AUV Properties', fieldsize, defaults);
                auv_userIn = AUV(model='Custom', mass=eval(auvSpecs{1}), missionSpecs=eval(auvSpecs{2}), maxBattery=eval(auvSpecs{3}), chargeRate=eval(auvSpecs{4}), chargeMethod=eval(auvSpecs{5}), hotelLoad=eval(auvSpecs{6}), rechargeThreshold=eval(auvSpecs{7}), n_battery=eval(auvSpecs{8}), n_powerTransfer=eval(auvSpecs{9}));
            else
                auv_userIn = AUV(model=auvModels{auvModelNum}); 
            end
    end
    
    
    % Specify & load wave resource data type
    resourceDataType = listdlg('PromptString',{'Specify the type of wave resource data you plan to use'}, 'SelectionMode','single', ...
        'ListString', {'Proteus struct containing power generated during different sea states','Time series of power generated','Value of mean power','Time series of wave specs (Hs, Te, Tp)','Mean Wave specs (Hs, Te, Tp)'});

end

dtSec = 30.0;  % [s]  % no prompt for this (yet) but is a model input

if strcmp(depVar, 'AUV Model')
    modIn = ModelInput(simHrs, depVar, resourceDataType, incorpStagger, maxFleetSize, userDefinedBattery, dtSec, auvModels);
elseif strcmp(depVar, 'WEC Power Gen / Wave Resource')
    modIn = ModelInput(simHrs, depVar, resourceDataType, incorpStagger, maxFleetSize, userDefinedBattery, dtSec, auvModels{auvModelNum});
else
    error('Unknown dependent variable')
end



% Load Resource Data
modIn.extractDataInfo;


%% Loop Preallocation, Constant Calcs, & Output Object Generation

% Generate output & WEC objects
modOut = ModelOutput(modIn.simTime, depVar, resourceDataType);

modOut.incorpStagger = incorpStagger;
modOut.maxFleetSize = maxFleetSize;

wec_userIn = WEC();

% Calculate dependent variable loop length
switch modIn.depVar
    case 'AUV Model'
        loopLength = length(modIn.auvModels); 
    case 'WEC Power Gen / Wave Resource' 
        switch resourceDataType
            case 1
                loopLength = length(RM3); 
            case 2
                loopLength = length(dataFiles);
            case 3
                loopLength = length(meanPower);
            case 4
                loopLength = length(dataFiles);
            case 5
                loopLength = length(waveSpecTable.SignificantWaveHeight);
            case 6
                loopLength = length(dataFiles);
        end
end


% Calculate Power Generation
if strcmp(modIn.depVar, 'AUV Model')
    modIn.calcPowerGen(wec_userIn, modOut)
end


%% Generate central energy storage object

if userDefinedBattery ~= 0
    energyStorage_userIn = EnergyStorage(maxBattery=userDefinedBattery); 

else  % userDefinedBattery == 0 - using code-calculated size
    energyStorage_userIn = EnergyStorage();  % Empty battery storage will be rewritten after fleet size determination...
end

modOut.userDefinedBattery = userDefinedBattery; 


%% Dependent Variable Loop / Run Simulation %% Could probably do this parallel..
% Preallocate wec & auv objects
switch modIn.depVar
    case 'AUV Model'
        auv = arrayfun(@(i) AUV('model', modIn.auvModels{i}), 1:numel(modIn.auvModels));
    case 'WEC Power Gen / Wave Resource'
        auv = repmat(auv_userIn, 1, 21);
end
wec = repmat(wec_userIn, 1, 21);
energyStorage = repmat(energyStorage_userIn, 1, 21);

% Preallocate local variables
fleetSize = zeros(1, loopLength); 
auvFleetLocal = cell(1, loopLength);
energyStorageBatteryLvl = zeros(size(modIn.simTime, 1), loopLength);
wecBatteryLvl = zeros(size(modIn.simTime, 1), loopLength);
auvBatteryLvl = cell(1, loopLength);
auvSchedule = cell(1, loopLength);
auvTimeOnMission = cell(1, loopLength);
auvTimeOnMissionCorrected = cell(1, loopLength);
centralBatteryCapacity = zeros(1, loopLength); 
ratePwrUsed = zeros(1, loopLength); 
auvMissionLength = zeros(1, loopLength); 
meanPowerGen = zeros(1, loopLength); 
switch modIn.depVar
    case 'AUV Model'
        powerGenMeans = zeros(size(modIn.simTime, 1), 1);  % only save power gen once if it isn't recalculated each iteration
    case 'WEC Power Gen / Wave Resource'
        powerGenMeans = zeros(size(modIn.simTime, 1), loopLength);
end

parfor depVarCount = 1:loopLength

    %% Fleet size calculation
    runFleetCalc = 1;
    fleetCalcCount = 0;
    while runFleetCalc == 1
        fleetCalcCount = fleetCalcCount + 1;

        if fleetCalcCount > 1
            fleetSize(depVarCount) = fleetSize(depVarCount) - 1;

            % clear values from prev simulation
            % clear modOut.energyStorageBatteryLvl modOut.wecBatteryLvl modOut.auvBatteryLvl modOut.auvSchedule modOut.auvTimeOnMission

        else
            fleetSize(depVarCount) = calcFleetSize(wec(depVarCount), auv(depVarCount), energyStorage(depVarCount), modIn.maxFleetSize);  % initial calculation

        end
        
        
        %% Update central storage hotel load to account for additional AUV docks
        energyStorage(depVarCount).hotelLoad = energyStorage(depVarCount).baseHotelLoad + energyStorage(depVarCount).dockHotelLoad*(fleetSize(depVarCount)); 


        %% Calculate minimum central battery, for battery to not limit fleet
        % Need at least enough energy saved in central battery to fully
        % recharge AUV(s), and support AUV & WEC hotel loads during
        % that time given poor power generation +~ 10% for safety.
        
        wec(depVarCount).calcLowPower(modIn.resourceDataType, modIn.dt, auv(depVarCount)); 

        lowPowerOverflow = wec(depVarCount).lowPowerGen -  wec(depVarCount).hotelLoad/(wec(depVarCount).n_battery^2);  % Assumes wec is already at max battery
        minBattery = ( auv(depVarCount).chargeLoad(auv(depVarCount).mission)*fleetSize(depVarCount)/energyStorage(depVarCount).n_battery + auv(depVarCount).chargeTime(auv(depVarCount).mission)*energyStorage(depVarCount).hotelLoad/energyStorage(depVarCount).n_battery/energyStorage(depVarCount).n_powerTrnsfr - lowPowerOverflow*auv(depVarCount).chargeTime(auv(depVarCount).mission)*energyStorage(depVarCount).n_battery*energyStorage(depVarCount).n_wecPwrTrnsfr  ) *1.05;  % Given lowest possible power generation during recharge OR 5% of WEC battery, if threshold is negative (i.e. if power gen > power draw)

        if minBattery < 0  % Somewhat arbituary minimum battery if pGen > operating loads
            warning('Problem with minimum battery calculation. Using default value for this simulation.')
            minBattery = 0.25 * ( auv(depVarCount).chargeLoad(auv(depVarCount).mission)*fleetSize(depVarCount)/energyStorage(depVarCount).n_battery + auv(depVarCount).chargeTime(auv(depVarCount).mission)*energyStorage(depVarCount).hotelLoad/energyStorage(depVarCount).n_battery/energyStorage(depVarCount).n_powerTrnsfr ); 
        end
 
        % Clear battery storage for new simulation if battery amount is not user-input
        if energyStorage(depVarCount).userDefinedBattery ~= 1
            energyStorage(depVarCount).maxBattery = minBattery;  % set central battery to minimum needed for AUV fleet

        elseif energyStorage(depVarCount).maxBattery < minBattery
            warning('Central energy storage amount is likely too low to support the fleet of %f %s AUV(s) determined by the given wave resource. Consider increasing central battery from %f to %f kWh.', fleetSize(depVarCount), auv(depVarCount).model, energyStorage(depVarCount).maxBattery*10e-3, minBattery*10e-3)
        
        elseif minBattery < energyStorage(depVarCount).maxBattery 
            warning('Central energy storage is %f kWh larger than what is required for this configuration.', (energyStorage(depVarCount).maxBattery - minBattery)*10e-3)

        end


        % Build fleet & run model if possible
        if fleetSize(depVarCount) > 0
    
            % Build fleet
            switch modIn.depVar
                case 'AUV Model'
                    auvFleet = arrayfun(@(x) AUV(model=modIn.auvModels{depVarCount}), 1:fleetSize(depVarCount));  
                case 'WEC Power Gen / Wave Resource'
                    auvFleet = arrayfun(@(x) AUV(model=modIn.auvModels), 1:fleetSize(depVarCount));
            end

            auvFleetLocal{depVarCount} = auvFleet;  % save to model output object. 
            
            %% Run Simulation
            % Simulation Goals: 
            % * Calculate battery levels of central storage & AUV(s)
            % * Generate AUV mission schedule given the following rules: 
            %   - AUV must dock with at least 20% of its battery left
            %   - Minimize time when central energy storage & AUV are both at full battery
            %   - energyStorage must have enough battery to charge AUV(s) when they return
            [energyStorageBatteryLvl(:, depVarCount), wecBatteryLvl(:, depVarCount), auvBatteryLvl{depVarCount}, auvSchedule{depVarCount}] = runPowerModel(modIn.simTime, wec(depVarCount), auvFleet, energyStorage(depVarCount), modIn.incorpStagger);
           
            %% Post-Simulation Calcs

            % Time spent 'on mission'
            auvTimeOnMission{depVarCount} = sum((auvSchedule{depVarCount} == 1) * modIn.dt);
            
            % Time 'on-mission' within performance calculation domain in the case of staggered deployments (excluding time AUVs are artificially held at dock)
            if modIn.incorpStagger == 1 
                if fleetSize(depVarCount) == 1
                    auvTimeOnMissionCorrected{depVarCount} = auvTimeOnMission{depVarCount};
                    
                else
                    staggerHours = (auvFleet(1).missionSpecs(2) + auvFleet(1).chargeTime) / fleetSize(depVarCount);
                    staggerPreliminaryTime = (fleetSize(depVarCount) - 1) * staggerHours;
                    [~, preDomainIndx] = min(abs(modIn.simTime - staggerPreliminaryTime));

                    auvTimeOnMissionCorrected{depVarCount} = sum((auvSchedule{depVarCount}(preDomainIndx+1:end, :) == 1) * modIn.dt);  % Exclude on-mission time before all AUVs are deployed
                end

            end

            %% Simulation Quality Checks

            % if last auv only went on one mission, and first auv went on more, OR central battery dropped below 0 (even with battery saver) auvFleet is too large by at least 1
            if ( (auvTimeOnMission{depVarCount}(1,end) - (auvFleet(end).chargeTime(auvFleet(end).mission) + auvFleet(end).missionSpecs(auvFleet(end).mission, 2)) ) <= 0 ) && ((auvTimeOnMission{depVarCount}(1,1) - (auvFleet(1).chargeTime(auvFleet(1).mission) + auvFleet(1).missionSpecs(auvFleet(1).mission, 2)) ) > 0 ) || any(energyStorageBatteryLvl(:,depVarCount) < 0)
                runFleetCalc = 1;  % re-run fleet size calculation. (Disable to plot results for systems with too-many AUVs)  
                warning('Sorry, our initial %s fleet size estimate of %f was too high. Re-running the simulation now...', auv(depVarCount).model, fleetSize(depVarCount))
        
            else
                runFleetCalc = 0;

                centralBatteryCapacity(depVarCount) = energyStorage(depVarCount).maxBattery;

    
            end 

        else  % Wave resource insufficient to support deployment of AUV here
            warning('Wave resource insufficient to support the hotel load for any %s AUVs at this site.', auv(depVarCount).model)
            
            auvTimeOnMission{depVarCount} = 0;  % No AUV's = no time spent on missions
            auvTimeOnMissionCorrected{depVarCount} = 0; 
            auvSchedule{depVarCount} = [];
            auvBatteryLvl{depVarCount} = [];
            energyStorageBatteryLvl(:, depVarCount) = zeros(size(modIn.simTime));
            wecBatteryLvl(:, depVarCount) = zeros(size(modIn.simTime)); 
            auvFleetLocal{depVarCount} = [];

            runFleetCalc = 0;
        end

    end  % while fleetCalc == 1 

    %% Values to Track

    % (Energy used [Wh] / mission time [h]) Homogenous auv fleets only!
    ratePwrUsed(depVarCount) = (auv(depVarCount).missionSpecs(auv(depVarCount).mission, 3)*auv(depVarCount).maxBattery + auv(depVarCount).hotelLoad*auv(depVarCount).chargeTime) / (auv(depVarCount).missionSpecs(auv(depVarCount).mission, 2) + auv(depVarCount).chargeTime);  % AVG rate AUV uses power during a mission + recharge cycle. No efficiencies applied!  
    auvMissionLength(depVarCount) = auv(depVarCount).missionSpecs(auv(depVarCount).mission, 2);
    meanPowerGen(depVarCount) = wec(depVarCount).meanPowerGen;
    switch modIn.depVar
        case 'AUV Model'
            if depVarCount == 1
                powerGenMeans(:, depVarCount) = wec(depVarCount).powerGenMeans;  % only save power gen once if it isn't recalculated each iteration
            end
        case 'WEC Power Gen / Wave Resource'
            powerGenMeans(:, depVarCount) = wec(depVarCount).powerGenMeans; 
    end

end  % depVar parfor

%% Save local variables to modOut object
modOut.fleetSize = fleetSize; 
modOut.auvFleet = auvFleetLocal; 
modOut.energyStorageBatteryLvl = energyStorageBatteryLvl; 
modOut.wecBatteryLvl = wecBatteryLvl; 
modOut.auvBatteryLvl = auvBatteryLvl; 
modOut.auvSchedule = auvSchedule; 
modOut.auvTimeOnMission = auvTimeOnMission; 
modOut.auvTimeOnMissionCorrected = auvTimeOnMissionCorrected; 
modOut.centralBatteryCapacity = centralBatteryCapacity; 
modOut.ratePwrUsed = ratePwrUsed; 
modOut.auvMissionLength = auvMissionLength; 
modOut.meanPowerGen = meanPowerGen; 
modOut.powerGenMeans = powerGenMeans; 

%% Figures
%{
modOut.plotBatteryTracker

modOut.plotTimeOnMission

modOut.plotPowerGen
%}

%% Save Output
% Edit and run the following to save data to 'outputData' folder: 
%{  
save('outputData/usWestCoast_03Dec25.mat','auv','energyStorage','modOut','modIn','wec')
%}

%% Simulation Meta
% profile off; profile viewer;
beep
toc