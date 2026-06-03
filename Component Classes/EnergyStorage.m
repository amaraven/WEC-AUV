% Ama Hartman

%% Define Class
classdef EnergyStorage < handle  
    % EnergyStorage defines the properties and methods (functions) saved in
    % EnergyStorage objects. These objects store information on the
    % 'central energy storage unit' properties (max battery, hotel load,
    % efficiencies, ...) and can estimate future battery levels of the unit. 
    %
    % To crate an EnergyStorage object named 'energyStorage', use the
    % following syntax: 
    % energyStorage = EnergyStorage(batteryCapacity, baseHotelLoad,
    %                 dockHotelLoad, n_battery, n_wecPowerTransfer, n_dockPowerTransfer);

    properties (GetAccess = public, SetAccess = private)
        baseHotelLoad   % [W] Hotel load of unit, excluding docks
        dockHotelLoad   % [W] Hotel load of a single AUV dock
        n_battery       % Battery efficiency [0.XX]
        n_wecPwrTrnsfr  % Efficiency of power transfer between WEC & central energy storage [0.XX]
        n_powerTrnsfr   % Efficiency of power transfer between central battery and docks [0.XX]
    end

    properties (Access = public)
        battery         % [Wh] Current battery level 
        batteryCapacity % [Wh] Total energy storage onboard 
        hotelLoad       % [W] Total baseline power usage of hardware (unit + AUV dock(s)). Defaults to base + 1 dock
    end

    properties (Dependent)
        lowBatteryLvl;  % [Wh] Battery level that triggers 'low power' mode (stops charging AUV(s), continues to supply hotel loads)
    end
    
    %% Instance Methods (need object as input)
    methods

        %% Constructor: Creates & returns an EnergyStorage object
        function energyStorage = EnergyStorage(batteryCapacity, baseHotelLoad, dockHotelLoad, n_battery, n_wecPowerTransfer, n_dockPowerTransfer)
            % EnergyStorage object constructor generates an object with the given
            % properties. Default values are used if no inputs are given.
            arguments
                batteryCapacity (1,1) {mustBeNumeric, mustBeNonnegative} = 0  % Total energy storage
                baseHotelLoad (1,1) {mustBeNumeric, mustBeNonnegative} = 10  % [W] Hotel load of unit, excluding AUV docks
                dockHotelLoad (1,1) {mustBeNumeric, mustBeNonnegative} = 20  % [W] Hotel load of a single AUV dock
                n_battery (1,1) {mustBeNumeric,  mustBeInRange(n_battery, 0, 1, 'exclude-lower')} = 0.9  % [0.XX] Battery efficiency
                n_wecPowerTransfer (1,1) {mustBeNumeric,  mustBeInRange(n_wecPowerTransfer, 0, 1, 'exclude-lower')} = 0.9  % [0.XX] Power transfer efficiency between the WEC & central battery
                n_dockPowerTransfer (1,1) {mustBeNumeric,  mustBeInRange(n_dockPowerTransfer, 0, 1, 'exclude-lower')} = 0.9  % [0.XX] Power transfer efficiency between the central battery and docks
            end

            energyStorage.batteryCapacity = batteryCapacity; 
            energyStorage.baseHotelLoad = baseHotelLoad;
            energyStorage.dockHotelLoad = dockHotelLoad;
            energyStorage.n_battery = n_battery;
            energyStorage.n_wecPwrTrnsfr = n_wecPowerTransfer;
            energyStorage.n_powerTrnsfr = n_dockPowerTransfer;

            % Default hotel load = base + 1 dock
            energyStorage.hotelLoad = energyStorage.baseHotelLoad + energyStorage.dockHotelLoad;  % Default value assumes 1 AUV

            % Initialize energyStorage battery at 85% of max, if given
            if ~( batteryCapacity == 0 || isempty(batteryCapacity) )
                energyStorage.battery = energyStorage.batteryCapacity * 0.85;
            else
                energyStorage.battery = [];
            end

        end  % constructor fn


        %% Dependent property calculations...
        function lowBatteryLvl = get.lowBatteryLvl(energyStorage)  
            % Sets lowBatteryLevel to 20% of the max capacity
            lowBatteryLvl = energyStorage.batteryCapacity * 0.20; 
        end  % low battery fn
        

        %% Re-Initialize energyStorage Object
        function reset(energyStorage)
            % Initialize energyStorage battery at 85% of max, if given
            if ~( energyStorage.batteryCapacity == 0 || isempty(energyStorage.batteryCapacity) )
                energyStorage.battery = energyStorage.batteryCapacity * 0.85;
            else
                energyStorage.battery = [];
            end

            % Default hotel load = base + 1 dock
            energyStorage.hotelLoad = energyStorage.baseHotelLoad + energyStorage.dockHotelLoad;  % Default value assumes 1 AUV
        end


        %% Battery Estimator
        function [eStorageFutureBattery, eStorageNoBatteryFlag] = estBattery(energyStorage, auvFleet, simTime, dt, auvToDeploy, wec)
            % Estimates the central battery level at a given time
            % Inputs: 
            % - auvFleet: Array of AUV objects
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
            auvMissionTime = zeros(1, numAUVs); 
            auvChargeTime = zeros(1, numAUVs);
            eStorageNoBatteryFlag = 0; 

            % Extract auvFleet values
            for auvNum = 1: numAUVs
                auv = auvFleet(auvNum);
                auvChargeRate(auvNum) = auv.chargeRate; 
                auvHotelLoad(auvNum) = auv.hotelLoad;
                auv_n_powerTransfer(auvNum) = auv.n_powerTransfer; 
                auv_n_battery(auvNum) = auv.n_battery; 
                auvOpState(auvNum) = auv.opState(1);
                auvOpTimeComplete(auvNum) = auv.opTimeComplete;
                auvMissionTime(auvNum) = auv.missionSpecs(auv.mission, 2); 
                auvChargeTime(auvNum) = auv.chargeTime;
            end  

            % Calculate a linear recharge rate which accounts for charging slowdown. (Total battery / (time to charge to 80% + time to charge last 20%))
            adjustedRechargeRate = auvChargeRate / 1.2;  % Simplified from: maxB / ( (0.8*maxB / chargeRate) + (0.4*maxB / chargeRate) ) = chargeRate / 1.2
           
            % Calculate time auv auvToDeploy will complete charging after return from mission
            t_return = simTime + auvFleet(auvToDeploy).missionSpecs(auvFleet(auvToDeploy).mission, 2) + auvFleet(auvToDeploy).chargeTime(auvFleet(auvToDeploy).mission); 

            % Power draw rates[W]
            hotelDrawRate = auvHotelLoad ./ auv_n_powerTransfer ./ auv_n_battery; 
            chargingDrawRate = (adjustedRechargeRate ./ auv_n_powerTransfer ./ auv_n_battery + hotelDrawRate);  % [W] - All auv efficiencies applied
            
            % Calculate power generation overflow during mission + charge
            wecPwrGenFx = mean(wec.powerGenMeans(int32(simTime/dt) : int32(min(t_return/dt, length(wec.powerGenMeans))) ));
            % Slightly faster (?) than using mean() each timestep this fn is called...
            % i1 = int32(simTime/dt);
            % i2 = int32(min(t_return/dt, length(wec.powerGenMeans)));
            % wecPwrGenFx = (wec.cumPwrGen(i2+1) - wec.cumPwrGen(i1)) / (i2 - i1 +1);  

            pwrOverflowRate = ( wecPwrGenFx - ( (wec.batteryCapacity - wec.battery)/(dt*wec.n_battery) + wec.hotelLoad/(wec.n_battery^2) ) ) * energyStorage.n_wecPwrTrnsfr * energyStorage.n_battery;
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
            
            chargingMatrix = zeros(numel(auvFleet), length(timeIntervals));  
            hotelMatrix = zeros(numel(auvFleet), length(timeIntervals));
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

                % Power dumped = PowerGenOverflow - eStorageHotelLoad - eStoragebatteryCapacity + eStorageCurrentBattery - auvTotalDraw
                powerDumps(j) = max(pwrOverflowRate.*timeIntervals(j) - energyStorage.hotelLoad.*timeIntervals(j)/energyStorage.n_battery/energyStorage.n_powerTrnsfr - energyStorage.batteryCapacity + eStorageBattInit - sum(chargingMatrix(:,j).*chargingDrawRate').*timeIntervals(j)/energyStorage.n_battery - sum(hotelMatrix(:,j).*hotelDrawRate')*timeIntervals(j)/energyStorage.n_battery, 0);

                % Calc new 'current' battery for start of next interval
                eStorageBattInit = eStorageBattInit + pwrOverflowRate.*timeIntervals(j) - energyStorage.hotelLoad.*timeIntervals(j)/energyStorage.n_battery/energyStorage.n_powerTrnsfr - sum(chargingMatrix(:,j).*chargingDrawRate')*timeIntervals(j)/energyStorage.n_battery - sum(hotelMatrix(:,j).*hotelDrawRate')*timeIntervals(j)/energyStorage.n_battery - powerDumps(j);  % est. wec battery at beginning of next interval (prev batt + power gen - power draw - power dumped)
            
                % Intermediate battery level flag
                if eStorageBattInit < 0
                    eStorageNoBatteryFlag = 1; 
                end

            end

            totalDraw = sum(sum(chargingMatrix.*(chargingDrawRate').*timeIntervals/energyStorage.n_battery) + sum(hotelMatrix.*(hotelDrawRate').*timeIntervals/energyStorage.n_battery));  % sum of auv draw for each interval - UNDERESTIMATES for any auv not completing their charge

            totPwrDump = sum(powerDumps); 
            
            eStorageFutureBattery = energyStorage.battery + totPwrOverflow - totPwrDump - energyStorage.hotelLoad*(t_return - simTime)/energyStorage.n_battery/energyStorage.n_powerTrnsfr - totalDraw;

            if eStorageFutureBattery < 0 
                eStorageNoBatteryFlag = 1; 

            end

        end  % battery estimator fn


    end  % instance methods

end