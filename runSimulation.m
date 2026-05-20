% Ama Hartman

% Function compiles and calculates all inputs needed to run simulation,
% calls runPowerModel.m to execute and compiles outputs into a ModelOutput
% class object 'modOut'.
%
% INPUTS: 
% * uiSettings - Struct containing all user-input settings including power
%   generation / resource data, simulation settings, hardware settings, and
%   output plot selections. 
% 
% OUTPUTS: 
% * modOut - object containing outputs of the power model simulation(s)
% * modIn - object containing select inputs of the power model
%   simulation(s) along with preprocessing code
% * wec - object containing wec properties
% * auv - object containing auv properties
% * energyStorage -



%
% OUTPUTS: 
% modOut contains...
% * Fleet size: Number of AUVs that can be supported by the WEC
% * Central battery size: size of the central battery required 
% * AUV mission schedule: Example schedule of AUV missions and recharge
% periods
% * Component battery level tracker: battery level history for each
%   component (WEC, central battery, and AUV(s)) during simulation time
% * Global comparison variable plots: compares system performance between
%   comparison variable cases

function simResults = runSimulation(auv, wec, energyStorage, modIn)

% Initialze results struct
simResults = struct('fleetSize',    0, ...
    'centralBatteryCapacity',       0, ...
    'energyStorageBatteryLvl',      [], ...
    'wecBatteryLvl',                [], ...
    'auvBatteryLvl',                [], ...
    'auvSchedule',                  [], ...
    'auvTimeOnMission',             0, ...
    'auvTimeOnMissionCorrected',    0, ...
    'auvFleet',                     [], ...
    'ratePwrUsed',                  0, ...
    'auvMissionLength',             0, ...
    'meanPowerGen',                 0, ...
    'powerGenMeans',                []);


%% Fleet size calculation
runFleetCalc = 1;
fleetCalcCount = 0;
while runFleetCalc == 1
    fleetCalcCount = fleetCalcCount + 1;

    if fleetCalcCount > 1
        simResults.fleetSize = simResults.fleetSize - 1;

    else
        simResults.fleetSize = calcFleetSize(wec, auv, energyStorage, modIn.maxFleetSize);  % initial calculation

    end
    
    
    %% Update central storage hotel load to account for additional AUV docks
    energyStorage.hotelLoad = energyStorage.baseHotelLoad + energyStorage.dockHotelLoad*(simResults.fleetSize); 


    %% Calculate minimum central battery, for battery to not limit fleet
    % Need at least enough energy saved in central battery to fully
    % recharge AUV(s), and support AUV & WEC hotel loads during
    % that time given poor power generation +~ 10% for safety.
    
    wec.calcLowPower(modIn.resourceDataType, modIn.dt, auv); 

    lowPowerOverflow = wec.lowPowerGen -  wec.hotelLoad/(wec.n_battery^2);  % Assumes wec is already at max battery
    minBattery = ( auv.chargeLoad(auv.mission)*simResults.fleetSize/energyStorage.n_battery + auv.chargeTime(auv.mission)*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTrnsfr - lowPowerOverflow*auv.chargeTime(auv.mission)*energyStorage.n_battery*energyStorage.n_wecPwrTrnsfr  ) *1.05;  % Given lowest possible power generation during recharge OR 5% of WEC battery, if threshold is negative (i.e. if power gen > power draw)

    if minBattery < 0  % Somewhat arbituary minimum battery if pGen > operating loads
        warning('Problem with minimum battery calculation. Using default value for this simulation.')
        minBattery = 0.25 * ( auv.chargeLoad(auv.mission)*simResults.fleetSize/energyStorage.n_battery + auv.chargeTime(auv.mission)*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTrnsfr ); 
    end

    % Clear battery storage for new simulation if battery amount is not user-input
    if modIn.userDefinedBattery ~= 1
        energyStorage.batteryCapacity = minBattery;  % set central battery to minimum needed for AUV fleet

    elseif energyStorage.batteryCapacity < minBattery
        warning('Central energy storage amount is likely too low to support the fleet of %f %s AUV(s) determined by the given wave resource. Consider increasing central battery from %f to %f kWh.', simResults.fleetSize, auv.model, energyStorage.batteryCapacity*10e-3, minBattery*10e-3)
    
    elseif minBattery < energyStorage.batteryCapacity 
        warning('Central energy storage is %f kWh larger than what is required for this configuration.', (energyStorage.batteryCapacity - minBattery)*10e-3)

    end


    % Build fleet & run model if possible
    if simResults.fleetSize > 0

        % Build fleet
        auvFleet = arrayfun(@(x) AUV(auv.model, auv.mass, auv.batteryCapacity, auv.missionSpecs(1,2), auv.missionSpecs(1,3), auv.hotelLoad, auv.chargeRate, auv.chargeMethod, auv.rechargeThreshold, auv.n_battery, auv.n_powerTransfer), 1:simResults.fleetSize);  

        simResults.auvFleet = auvFleet;  % save to output struct. 
        
        %% Run Simulation
        % Simulation Goals: 
        % * Calculate battery levels of central storage & AUV(s)
        % * Generate AUV mission schedule given the following rules: 
        %   - AUV must dock with at least 20% of its battery left
        %   - Minimize time when central energy storage & AUV are both at full battery
        %   - energyStorage must have enough battery to charge AUV(s) when they return
        [simResults.energyStorageBatteryLvl, simResults.wecBatteryLvl, simResults.auvBatteryLvl, simResults.auvSchedule] = runPowerModel(modIn.simTime, wec, auvFleet, energyStorage, modIn.incorpStagger);
       
        %% Post-Simulation Calcs

        % Time spent 'on mission'
        simResults.auvTimeOnMission = sum((simResults.auvSchedule == 1) * modIn.dt);
        
        % Time 'on-mission' within performance calculation domain in the case of staggered deployments (excluding time AUVs are artificially held at dock)
        if modIn.incorpStagger == 1 
            if simResults.fleetSize == 1
                simResults.auvTimeOnMissionCorrected = simResults.auvTimeOnMission;
                
            else
                staggerHours = (simResults.auvFleet(1).missionSpecs(2) + simResults.auvFleet(1).chargeTime) / simResults.fleetSize;
                staggerPreliminaryTime = (simResults.fleetSize - 1) * staggerHours;
                [~, preDomainIndx] = min(abs(modIn.simTime - staggerPreliminaryTime));

                simResults.auvTimeOnMissionCorrected = sum((simResults.auvSchedule(preDomainIndx+1:end, :) == 1) * modIn.dt);  % Exclude on-mission time before all AUVs are deployed
            end

        end

        %% Simulation Quality Checks

        % if last auv only went on one mission, and first auv went on more, OR central battery dropped below 0 (even with battery saver) auvFleet is too large by at least 1
        if ( (simResults.auvTimeOnMission(1,end) - (auvFleet(end).chargeTime(auvFleet(end).mission) + auvFleet(end).missionSpecs(auvFleet(end).mission, 2)) ) <= 0 ) && ((simResults.auvTimeOnMission(1,1) - (auvFleet(1).chargeTime(auvFleet(1).mission) + auvFleet(1).missionSpecs(auvFleet(1).mission, 2)) ) > 0 ) || any(simResults.energyStorageBatteryLvl < 0)
            runFleetCalc = 1;  % re-run fleet size calculation. (Disable to plot results for systems with too-many AUVs)  
            warning('Sorry, our initial %s fleet size estimate of %f was too high. Re-running the simulation now...', auv.model, simResults.fleetSize)
    
        else
            runFleetCalc = 0;

            simResults.centralBatteryCapacity = [energyStorage.batteryCapacity];


        end 

    else  % Wave resource insufficient to support deployment of AUV here
        warning('Wave resource insufficient to support the hotel load for any %s AUVs at this site.', auv.model)
        
        simResults.auvTimeOnMission = 0;  % No AUV's = no time spent on missions
        simResults.auvTimeOnMissionCorrected = 0; 
        simResults.auvSchedule = [];
        simResults.auvBatteryLvl = [];
        simResults.energyStorageBatteryLvl = zeros(size(modIn.simTime));
        simResults.wecBatteryLvl = zeros(size(modIn.simTime)); 
        simResults.auvFleet = [];

        runFleetCalc = 0;
    end

end  % while fleetCalc == 1 

%% Values to Track

% (Energy used [Wh] / mission time [h]) Homogenous auv fleets only!
simResults.ratePwrUsed = (auv.missionSpecs(auv.mission, 3)*auv.batteryCapacity + auv.hotelLoad*auv.chargeTime) / (auv.missionSpecs(auv.mission, 2) + auv.chargeTime);  % AVG rate AUV uses power during a mission + recharge cycle. No efficiencies applied!  
simResults.auvMissionLength = auv.missionSpecs(auv.mission, 2);
simResults.meanPowerGen = wec.meanPowerGen;
simResults.powerGenMeans = wec.powerGenMeans;  % maybe only save power gen once if it isn't recalculated each iteration?


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

% %% Simulation Meta
% profile off; profile viewer;
% beep
% toc

end
