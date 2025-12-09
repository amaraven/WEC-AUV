% Ama Hartman - Su 2025

%% Define Class
classdef EnergyStorage < handle  
    % EnergyStorage defines the properties and methods (functions) saved in
    % EnergyStorage objects. These objects store information on the
    % 'central energy storage unit' properties (max battery, hotel load,
    % efficiencies, ...) and can estimate future battery levels of the unit. 
    %
    % To crate an EnergyStorage object named 'energyStorage', use the
    % following syntax: 
    % energyStorage = EnergyStorage(maxBattery, baseHotelLoad, dockHotelLoad)

    properties
        maxBattery      % [Wh] Total energy storage onboard 
        battery         % [Wh] Current battery level 
        userDefinedBattery  % 1/0 - User input max battery / model outputs battery capacity
        baseHotelLoad   % [W] Hotel load of unit, excluding docks
        dockHotelLoad   % [W] Hotel load of a single AUV dock
        hotelLoad       % [W] Total baseline power usage of hardware (unit + AUV dock(s)). Defaults to base + 1 dock
        n_battery       % Battery efficiency [0.XX]
        n_wecPwrTrnsfr  % Efficiency of power transfer between WEC & central energy storage [0.XX]
        n_powerTransfer % Efficiency of power transfer between central battery and docks [0.XX]
    end
    properties (Dependent)
        lowBatteryLvl;  % [Wh] Battery level that triggers 'low power' mode (stops charging AUV(s), continues to supply hotel loads)
    end
    
    %% Instance Methods (need object as input)
    methods

        %% Constructor: Creates & returns an EnergyStorage object
        function energyStorage = EnergyStorage(maxBattery, baseHotelLoad, dockHotelLoad)
            % Constructor function user-inputs
            energyStorage.maxBattery = maxBattery; 
            energyStorage.baseHotelLoad = baseHotelLoad; 
            energyStorage.dockHotelLoad = dockHotelLoad;

            % Default efficiency values
            energyStorage.n_battery = 0.90;  % 0.90 standard used in Chen '24
            energyStorage.n_wecPwrTrnsfr = 0.9;  % https://www.researchgate.net/publication/264124522_Efficiency_Comparison_of_Wire_and_Wireless_Battery_Charging_Based_on_Connection_Probability_Analysis   
            energyStorage.n_powerTransfer = 0.9; % ^" assumes cabled connection. 

            % Default hotel load = baseline + 1 dock
            energyStorage.hotelLoad = baseHotelLoad + dockHotelLoad;  % Default value assumes 1 AUV

            % Set userDefinedBattery flag
            if isempty(energyStorage.maxBattery)
                energyStorage.userDefinedBattery = 0; 
            else
                energyStorage.userDefinedBattery = 1; 
            end

        end  % constructor fn


        %% Dependent property calculations...
        function lowBatteryLvl = get.lowBatteryLvl(energyStorage)  
            % Sets lowBatteryLevel to 20% of the max capacity
            lowBatteryLvl = energyStorage.maxBattery * 0.20; 

        end  % low battery fn
        

        %% Battery Estimator
        function [eStorageFutureBattery, eStorageNoBatteryFlag] = estBattery(energyStorage, auvFleet, simTime, dt, auvToDeploy, wec)
            % Estimates the central battery level at a given time
            % Inputs: 
            % - auvFleet: Cell array of AUV objects
            % - simTime: [hr] Current simulation time
            % - dt: [hr] Simulation timestep
            % - auvToDeploy: index of AUV to deploy in auvFleet
            % - wec: WEC object
            % Outputs: 
            % - eStorageFutureBattery: Estimated energy storage battery
            %   level at the time auvToDeploy finishes charging after
            %   returning from its mission
            % - eStorageNoBatteryFlag: Battery level flag. (1) 'on' if 
            %   battery is expected to drop at or below 5% of max

            numAUVs = numel(auvFleet); 
            auvChargeRate = zeros(1, numAUVs); 
            auvHotelLoad = zeros(1, numAUVs); 
            auv_n_powerTransfer = zeros(1, numAUVs); 
            auv_n_battery = zeros(1, numAUVs); 
            auvOpState =zeros(1, numAUVs); 
            auvOpTimeComplete = zeros(1, numAUVs); 
            auvMission = zeros(1, numAUVs); 
            auvMissionBattUsed = zeros(1, numAUVs); 
            auvMissionTime = zeros(1, numAUVs); 
            auvMaxBatt = zeros(1, numAUVs);
            auvChargeTime = zeros(1, numAUVs);
            eStorageNoBatteryFlag = 0; 

            % Extract auvFleet values
            for auvNum = 1: numAUVs
                auv = auvFleet{auvNum};
                auvChargeRate(auvNum) = auv.chargeRate; 
                auvHotelLoad(auvNum) = auv.hotelLoad;
                auv_n_powerTransfer(auvNum) = auv.n_powerTransfer; 
                auv_n_battery(auvNum) = auv.n_battery; 
                auvOpState(auvNum) = auv.opState(1);
                auvOpTimeComplete(auvNum) = auv.opTimeComplete;
                auvMission(auvNum) = auv.mission;
                auvMissionBattUsed(auvNum) = auv.missionSpecs(auv.mission, 3) * auv.maxBattery;
                auvMissionTime(auvNum) = auv.missionSpecs(auv.mission, 2); 
                auvMaxBatt(auvNum) = auv.maxBattery;
                auvChargeTime(auvNum) = auv.chargeTime;
            end  

            % Calculate a linear recharge rate which accounts for charging slowdown. (Total battery / (time to charge to 80% + time to charge last 20%))
            adjustedRechargeRate = auvChargeRate / 1.2;  % Simplified from: maxB / ( (0.8*maxB / chargeRate) + (0.4*maxB / chargeRate) ) = chargeRate / 1.2
           
            % Calculate time auv auvToDeploy will complete charging after return from mission
            t_return = simTime + auvFleet{auvToDeploy}.missionSpecs(auvFleet{auvToDeploy}.mission, 2) + auvFleet{auvToDeploy}.chargeTime(auvFleet{auvToDeploy}.mission); 

            % Power draw rates[W]
            hotelDrawRate = auvHotelLoad ./ auv_n_powerTransfer ./ auv_n_battery; 
            chargingDrawRate = (adjustedRechargeRate ./ auv_n_powerTransfer ./ auv_n_battery + hotelDrawRate);  % [W] - All auv efficiencies applied
            
            % Calculate power generation overflow during mission + charge
            wecPwrGenFx = mean(wec.powerGenMeans(int32(simTime/dt) : int32(min(t_return/dt, length(wec.powerGenMeans))) ));
            pwrOverflowRate = ( wecPwrGenFx - ( (wec.maxBattery - wec.battery)/(dt*wec.n_battery) + wec.hotelLoad/(wec.n_battery^2) ) ) * energyStorage.n_wecPwrTrnsfr * energyStorage.n_battery;
            totPwrOverflow = pwrOverflowRate * (t_return - simTime);


            % Account for power dumping

            % AUV opState logic switches
            onMissionChargeLogic = (auvOpState == 1) .* (auvOpTimeComplete < t_return) == 1;  % If currently on a mission but will return and charge before auvToDeploy returns
            chargingLogic = auvOpState == 2;
            chargingHotelLogic = logical(chargingLogic .* ( auvOpTimeComplete < t_return ));  % Charging & will be done before auvToDploy returns
            onMissionHotelLogic = logical( onMissionChargeLogic .* (auvChargeTime + auvOpTimeComplete < t_return) );  % If on mission, will return before auvToDeploy returns, and will fully recharge before auvToDeploy returns
            
            % For an auv on mission that will recharge before auvToDeploy returns, this yields the timestamp it will be back and fully recharged
            timeDoneChargingAfterMission = auvChargeTime(onMissionHotelLogic) + auvOpTimeComplete(onMissionHotelLogic);

            % Track timestamps for changes in fleet op state
            missionOpChangeTimes = [auvOpTimeComplete(onMissionChargeLogic), timeDoneChargingAfterMission];  % [times auv(s) will go from on mission to charging, times auv on mission will go from charging to docked] (between simTime & t_return)
            chargingOpChangeTimes = auvOpTimeComplete(chargingHotelLogic);  
            auvToDeployOpChangeTimes = [simTime + auvMissionTime(auvToDeploy), simTime + auvMissionTime(auvToDeploy) + auvChargeTime(auvToDeploy)];  % [time back from mission, time done charging]
            fleetOpChangeTimes = sort([missionOpChangeTimes, chargingOpChangeTimes, auvToDeployOpChangeTimes]);  % all auv operation change times sorted from low to high. 
            
            if any(fleetOpChangeTimes < simTime) || any(fleetOpChangeTimes > t_return)
                error('Power dump estimation fault. One of the fleet change times is outside of expected range.')
            end
                
            fleetOpChangeTimes = [simTime, fleetOpChangeTimes];
            timeIntervals = diff(fleetOpChangeTimes);  % durations between changes in opState
            
            % Track AUV opStates (within intervals between changes in fleet opState)

            chargingMatrix = zeros(length(auvFleet), length(timeIntervals));  
            hotelMatrix = zeros(length(auvFleet), length(timeIntervals));
            powerDumps = zeros(1, length(timeIntervals));
            
            for j = 1:length(timeIntervals)  
                if any(onMissionChargeLogic)
                    chargingMatrix(logical((auvOpTimeComplete.*onMissionChargeLogic <= fleetOpChangeTimes(j)) .* ((auvChargeTime + auvOpTimeComplete).*onMissionChargeLogic > fleetOpChangeTimes(j))) , j) = 1;  % was on mission, now charging (for interval t_1-t_2: opTimeComplete <= t_1 && rechargeTimeComplete > t_1 )
                end

                chargingMatrix( (auvOpTimeComplete.*chargingLogic > fleetOpChangeTimes(j)), j) = 1;  % started charging & current time < time charging will be complete
                
                if (simTime + auvMissionTime(auvToDeploy) ) <= fleetOpChangeTimes(j) && fleetOpChangeTimes(j) < t_return  % if done w/ mission but not done charging...
                    chargingMatrix(auvToDeploy, j) = 1;
                end 

                hotelMatrix( (auvChargeTime + auvOpTimeComplete).*onMissionHotelLogic <= fleetOpChangeTimes(j), j) = 1;  % started on mission, finished mission, finished charging, now on hotel.
                hotelMatrix( logical((auvOpTimeComplete <= fleetOpChangeTimes(j)).*chargingLogic) , j) = 1;  % started charging, finished charging
                hotelMatrix( auvOpState == 3, j) = 1;  % started docked, stays docked
                hotelMatrix(auvToDeploy, j) = 0;  

                % calc power dumped during interval of constant fleet opState

                if j == 1
                    eStorageBattInit = energyStorage.battery;  % initialize battery level
                end

                % Power dumped = PowerGenOverflow - eStorageHotelLoad - eStorageMaxBattery + eStorageCurrentBattery - auvTotalDraw
                powerDumps(j) = max(pwrOverflowRate.*timeIntervals(j) - energyStorage.hotelLoad.*timeIntervals(j)/energyStorage.n_battery/energyStorage.n_powerTransfer - energyStorage.maxBattery + eStorageBattInit - sum(chargingMatrix(:,j).*chargingDrawRate').*timeIntervals(j)/energyStorage.n_battery - sum(hotelMatrix(:,j).*hotelDrawRate')*timeIntervals(j)/energyStorage.n_battery, 0);

                % Calc new 'current' battery for start of next interval
                eStorageBattInit = eStorageBattInit + pwrOverflowRate.*timeIntervals(j) - energyStorage.hotelLoad.*timeIntervals(j)/energyStorage.n_battery/energyStorage.n_powerTransfer - sum(chargingMatrix(:,j).*chargingDrawRate')*timeIntervals(j)/energyStorage.n_battery - sum(hotelMatrix(:,j).*hotelDrawRate')*timeIntervals(j)/energyStorage.n_battery - powerDumps(j);  % est. wec battery at beginning of next interval (prev batt + power gen - power draw - power dumped)
            
                % Intermediate battery level flag
                if eStorageBattInit < 0
                    eStorageNoBatteryFlag = 1; 
                end

            end

            totalDraw = sum(sum(chargingMatrix.*(chargingDrawRate').*timeIntervals/energyStorage.n_battery) + sum(hotelMatrix.*(hotelDrawRate').*timeIntervals/energyStorage.n_battery));  % sum of auv draw for each interval - UNDERESTIMATES for any auv not completing their charge

            totPwrDump = sum(powerDumps); 
            
            eStorageFutureBattery = energyStorage.battery + totPwrOverflow - totPwrDump - energyStorage.hotelLoad*(t_return - simTime)/energyStorage.n_battery/energyStorage.n_powerTransfer - totalDraw;

            if eStorageFutureBattery < 0 
                eStorageNoBatteryFlag = 1; 

            end

        end  % battery estimator fn


    end  % instance methods

end