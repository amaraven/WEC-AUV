% runConfiguration compiles and calculates all inputs needed to run the
% model, calls runPowerModel.m to execute, and compiles outputs into a 
% struct named 'simResults'.
%
% INPUTS: 
% - auv: AUV object or array of objects
% - wec: WEC object
% - energyStorage: EnergyStorage object 
% - modIn: ModelInput object containing all user-input parameters
% 
% OUTPUTS: 
% - simResults: Struct containing the following simulation outputs: 
%   - fleetSize: Size of AUV fleet
%   - centralBatteryCapacity: [Wh] total capacity of the central battery
%   - energyStorageBatteryLvl: [Wh] (mx1) array of central battery level
%     with simulation time
%   - wecBatteryLvl: [Wh] (mx1) array of WEC battery level with simulation 
%     time
%   - auvBatteryLvl: [Wh] (mxn, w/ n = fleetSize) array of AUV battery
%     level(s) with simulation time
%   - auvSchedule: (mxn, w/ n = fleetSize) array of AUV operational states
%     with simulation time
%   - auvTimeOnMission: (1xn w/ n = fleetSize) array of the amount of time
%     each AUV spends 'on-mission'
%   - auvTimeOnMissionCorrected: (1xn ...) array of the amount of time each
%     AUV spends 'on-mission' excluding time before all AUVs are deployed
%     at the start of the simulation
%   - auvFleet: (1xn ...) array of AUV objects
%   - meanPowerGen: [W] mean power generation throughout simulation
%   - powerGenMeans: [W]
%
% Ama Hartman

function runConfig(auv, wec, energyStorage, modIn, simResults)

%% Build wec fleet if empty (ocurrs when modIn.simGoal == 1)
if isempty(simResults.wecFleet)
    simResults.numWECs = 1;
    simResults.wecFleet = wec;
    updateWECAggregates(simResults, modIn, auv);
end


%% Fleet size calculation
reRunSim = 1;
simAttemptNum = 0;
while reRunSim == 1
    simAttemptNum = simAttemptNum + 1;

    if simAttemptNum > 10
        error('runConfig: Exceeded maximum retry attempts (10). Central battery repeatedly failed quality checks — check WEC generation estimates and fleet sizing logic.');
    end

    if modIn.simGoal == 1
        if simAttemptNum > 1
            simResults.fleetSize = simResults.fleetSize - 1;
    
            % reset wec & auv timestep parameters
            % wec.reset();
            arrayfun(@(w) w.reset(), simResults.wecFleet);
            arrayfun(@(a) a.reset(), auv); 
            energyStorage.reset();
        else
            simResults.fleetSize = calcFleetSize(simResults, auv, energyStorage, modIn.maxFleetSize);  % initial calculation
        end

    elseif modIn.simGoal == 2
        % Set AUV fleet size (given via user-input)
        simResults.fleetSize = numel(auv);

        % if re-running simulation, use a larger power generation
        if simAttemptNum > 1
            % newMinPwr = wec.meanPowerGen * 1.10;  % Increase power gen. by 10%
            newMinPwr = modIn.resourceDataVars.minPowerRequired * 1.10;  % Increase power gen. requirement by 10%
            modIn.addResourceDataVar('minPowerRequired', newMinPwr);
            modIn.calcPowerGen(wec, simResults)

            % reset wec & auv timestep parameters
            arrayfun(@(w) w.reset(), simResults.wecFleet);
            arrayfun(@(a) a.reset(), auv); 
            energyStorage.reset();
        end
    end
    
    
    %% Update central storage hotel load to account for additional AUV docks
    energyStorage.hotelLoad = energyStorage.baseHotelLoad + energyStorage.dockHotelLoad*(simResults.fleetSize); 

    
    %% Calculate wec fleet aggregate-power values
    % % Pre-simulation wec fleet calculations
    updateWECAggregates(simResults, modIn, auv);
       

    %% Calculate minimum central battery, for battery to not limit fleet
    % Need at least enough energy saved in central battery to fully
    % recharge AUV(s), and support AUV & WEC hotel loads during
    % that time given poor power generation +~ 20% for safety.
    lowPowerOverflow = simResults.lowPowerGen -  simResults.aggWECHotelLoad;  % Assumes wec is already at max battery
    % lowPowerOverflow = wec.lowPowerGen -  wec.hotelLoad/(wec.n_battery^2);  % Assumes wec is already at max battery
    switch modIn.simGoal
        case 1
            minBattery = ( auv.chargeLoad*simResults.fleetSize/energyStorage.n_battery + auv.chargeTime*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTrnsfr - lowPowerOverflow*auv.chargeTime*energyStorage.n_battery*energyStorage.n_wecPwrTrnsfr  )/0.8;% *1.05;  % Given lowest possible power generation during recharge OR 5% of WEC battery, if threshold is negative (i.e. if power gen > power draw)
        case 2
            minBattery = ( sum([auv.chargeLoad])/energyStorage.n_battery + mean([auv.chargeTime])*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTrnsfr - lowPowerOverflow*mean([auv.chargeTime])*energyStorage.n_battery*energyStorage.n_wecPwrTrnsfr  )/0.8;% *1.05; 
    end

    if minBattery < 0  % Somewhat arbituary minimum battery if pGen > operating loads
        warning('Power generation outweighs maximum expected operational loads for this system, using a default central battery capacity for this simulation.')
        switch modIn.simGoal
            case 1
                minBattery = 0.25 * ( auv.chargeLoad*simResults.fleetSize/energyStorage.n_battery + auv.chargeTime*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTrnsfr ); 
            case 2
                minBattery = 0.25 * ( sum([auv.chargeLoad])/energyStorage.n_battery + mean([auv.chargeTime])*energyStorage.hotelLoad/energyStorage.n_battery/energyStorage.n_powerTrnsfr );
        end
    end

    % Set central battery capacity
    if modIn.userDefinedBattery == 0  
        % No user-defined battery, set capacity to calculated minimum needed for AUV fleet
        energyStorage.batteryCapacity = minBattery;

    elseif energyStorage.batteryCapacity < minBattery
        switch modIn.simGoal 
            case 1
                warning('Central energy storage amount is likely too low to support the fleet of %f %s AUV(s) determined by the given wave resource. Consider increasing central battery from %f to %f kWh.', simResults.fleetSize, auv.model, energyStorage.batteryCapacity*10e-3, minBattery*10e-3)
            case 2
                warning('Central energy storage amount is likely too low to support the given AUV fleet. Consider increasing central battery from %f to %f kWh.', energyStorage.batteryCapacity*10e-3, minBattery*10e-3)
        end

    elseif minBattery < energyStorage.batteryCapacity 
        warning('Central energy storage is %f kWh larger than what is required for this configuration.', (energyStorage.batteryCapacity - minBattery)*10e-3)

    end

    % Modify central battery (optional toggle for analysis)
    % energyStorage.batteryCapacity = energyStorage.batteryCapacity*1.75;


    % Build fleet & run model if possible
    if simResults.fleetSize > 0

        % Build fleet
        switch modIn.simGoal
            case 1
                auvFleet = arrayfun(@(~) AUV(auv.model, auv.mass, auv.batteryCapacity, auv.missionTime, auv.missionBattPrcnt, auv.hotelLoad, auv.chargeRate, auv.chargeMethod, auv.rechargeThreshold, auv.n_battery, auv.n_powerTransfer), 1:simResults.fleetSize);  
            case 2
                auvFleet = auv;  % (user-input AUV fleet saved in app.auv then used as input to runConfig.m)
        end

        simResults.auvFleet = auvFleet;  % save to output struct. 
        
        %% Run Simulation
        % Simulation Goals: 
        % * Calculate battery levels of central storage & AUV(s)
        % * Generate AUV mission schedule given the following rules: 
        %   - AUV must dock with at least 20% of its battery left
        %   - Minimize time when central energy storage & AUV are both at full battery
        %   - energyStorage must have enough battery to charge AUV(s) when they return
        [simResults.energyStorageBatteryLvl, simResults.wecBatteryLvl, simResults.auvBatteryLvl, simResults.auvSchedule] = simulateSystemOps(modIn.simTime, simResults, auvFleet, energyStorage, modIn.incorpStagger);
       
        %% Post-Simulation Calcs

        % Time spent 'on mission'
        simResults.auvTimeOnMission = sum((simResults.auvSchedule == 1) * modIn.dt);
        
        % Time 'on-mission' within performance calculation domain in the case of staggered deployments (excluding time AUVs are artificially held at dock)
        if modIn.incorpStagger == 1 
            if simResults.fleetSize == 1
                simResults.auvTimeOnMissionCorrected = simResults.auvTimeOnMission;
                
            else
                staggerHours = mean( ([auvFleet.missionTime]+[auvFleet.chargeTime]) ./ length(auvFleet) );
                staggerPreliminaryTime = (simResults.fleetSize - 1) * staggerHours;
                [~, preDomainIndx] = min(abs(modIn.simTime - staggerPreliminaryTime));

                simResults.auvTimeOnMissionCorrected = sum((simResults.auvSchedule(preDomainIndx+1:end, :) == 1) * modIn.dt);  % Exclude on-mission time before all AUVs are deployed
            end

        end

        %% Simulation Quality Checks (fleet size & central battery level)
        
        switch modIn.simGoal
            case 1
                % if last auv only went on one mission, and first auv went on more, OR central battery dropped below 0 (even with battery saver) auvFleet is too large by at least 1
                if ( (simResults.auvTimeOnMission(1,end) - (auvFleet(end).chargeTime + auvFleet(end).missionTime) ) <= 0 ) && ((simResults.auvTimeOnMission(1,1) - (auvFleet(1).chargeTime + auvFleet(1).missionTime) ) > 0 ) || any(simResults.energyStorageBatteryLvl < 0)
                    reRunSim = 1;  % re-run fleet size calculation. (Disable to plot results for systems with too-many AUVs)  
                    warning('Sorry, our initial %s fleet size estimate of %f was too high. Re-running the simulation now...', auv.model, simResults.fleetSize)
            
                else
                    reRunSim = 0;
        
                    simResults.centralBatteryCapacity = [energyStorage.batteryCapacity];
                end 

            case 2
                % If central battery dropped below 0, or last AUV never
                % went on mission, powerGen insufficient. Run again w/
                % higher power generation value %%%%%%%%%%%%%%% todo
                if any(simResults.energyStorageBatteryLvl < 0)
                    reRunSim = 1;
                    warning('Central battery dropped below zero. Re-running simulation with a 10% higher minimum power generation.')
                else
                    reRunSim = 0; 

                    simResults.centralBatteryCapacity = [energyStorage.batteryCapacity];
                end              
        end

    else  % Wave resource insufficient to support deployment of AUV here
        switch modIn.simGoal
            case 1
                warning('Wave resource insufficient to support the hotel load for any %s AUVs at this site.', auv.model)
                
                simResults.auvTimeOnMission = 0;  % No AUV's = no time spent on missions
                simResults.auvTimeOnMissionCorrected = 0; 
                simResults.auvSchedule = [];
                simResults.auvBatteryLvl = [];
                simResults.energyStorageBatteryLvl = zeros(size(modIn.simTime));
                simResults.wecBatteryLvl = zeros(size(modIn.simTime)); 
                simResults.auvFleet = [];
        
                reRunSim = 0;
        end
    end

end  % while fleetCalc == 1 

end

function lowPowerGen = calcLowPowerGen(simResults, modIn, auv)
% Calculates aggregate low power generation for the wec fleet
% 
% INPUTS: 
% simResults: SimResults object with wecFleet generation values
% modIn: ModelInput object
% auv: AUV object 

switch modIn.resourceDataType
    case {1, 2, 4, 6}
        chargeWindow = ceil([auv.chargeTime]./modIn.dt);
        if numel(chargeWindow) > 1
            chargeWindow = mean(chargeWindow);
        end
        tmpPwrGen = movmean(simResults.powerGenMeans, chargeWindow);
        lowPowerGen = min(tmpPwrGen(floor(0.5*chargeWindow) : end - floor(0.5*chargeWindow)) );  % Exclude truncated means

    otherwise
        lowPowerGen = 0.75*simResults.meanPowerGen;
end

end


function updateWECAggregates(simResults, modIn, auv)
% Updates WEC object parameters in wecFleet

simResults.powerGenMeans = sum([simResults.wecFleet.powerGenMeans], 2);
simResults.aggWECHotelLoad = sum([simResults.wecFleet.hotelLoad]./([simResults.wecFleet.n_battery].^2));
simResults.meanPowerGen = mean(simResults.powerGenMeans);
simResults.lowPowerGen = calcLowPowerGen(simResults, modIn, auv); %%%%%%%%%%%%%%%%%%%%%% move low power calculation to helper function

end