% AUV defines the properties and methods (functions) saved in AUV objects.
% These objects store information on the Unmanned Underwater Vehicle's
% properties (model ID, mass, battery capacity, endurance, charging 
% rate, etc...), and can compute values relating to their battery level
% and mission duration. 
%
% To create a AUV object named 'auv', use the following syntax. 
% auv = AUV(model, mass, batteryCapacity, missionTime,
%       missionBatteryUsed, hotelLoad, chargeRate, chargeMethod, 
%       batteryRechargeThreshold, n_battery, n_powerTransfer);
%
% Ama Hartman - 2026


classdef AUV < handle   

    properties (GetAccess = public, SetAccess = public)  %%%%%%%%%%%% change set access to private ? %%%%%%%%%%%%%
        % User-Defined
        model           % Model of auv (i.e. 'Iver3')
        mass            % [kg] Mass of auv 
        % mission         % Mission number currently using / running
        % missionSpecs    % (1) time to complete [hr], and (2) percentage of battery used (0.XX)
        missionTime     % [hr] time to complete mission
        missionBattPrcnt % [0.XX] Percentage of battery used during mission
        batteryCapacity % [Wh] Total energy storage onboard 
        chargeRate      % [W] Rate of battery charge 
        chargeMethod    % Wired (1) or Wireless (0) charging
        hotelLoad       % [W] Baseline power usage of AUV
        rechargeThreshold   % [0.XX] Percentage of battery to recharge at. Default is 20%
        n_powerTransfer     % [0.XX] Power transfer efficiency 
        n_battery           % [0.XX] Battery efficiency (losses during charge & discharge)
        
        % Dependent
        chargeTime      % [hr] Column-vector of time to go to 100% charge after respective mission. NOTE: Mission must use >= 20% of battery for this calc to be accurate. 
        chargeLoad      % [Wh] Column-vector of total load put on WEC to charge auv after a given mission, assuming AUV started with a full battery
    end

    % Simulation variables
    properties (Access = public)
        opState         % Vector containing current operational state and simulation timestamp of respective state change in hours ([opState, t]) with 1) Executing AUV mission, 2) AUV recharging, 3) AUV docked & fully charged
        opTimeComplete  % [hr] Time current opState will be complete
        battery         % Vector containing battery level at the end of the last operational state (battery(1)) and the current battery level (battery(2)) [Wh] ([B(i-1); B(i)]). battery(1) corresponds to timestamp saved as opState(2)
        batteryTime     % [hr] Time corresponding to the current battery level (battery(2)) - NOTE: Not currently used externally as of 8/14/25    
        wecBatteryDraw  % [W] Energy draw level (0 if not charging) associated with current battery and batteryTime
    end 


    %% Instance Methods (need object as an input)
    methods

        %% Constructor: Creates & returns an object
        function auv = AUV(model, mass, batteryCapacity, missionTime, missionBatteryUsed, hotelLoad, chargeRate, chargeMethod, batteryRechargeThreshold, n_battery, n_powerTransfer)
            % AUV object constructor generates an object with the given properties. Default values are used if no inputs are given.
            arguments
            model (1,1) string = "Default AUV"  % String containing model name
            mass (1,1) {mustBeNumeric, mustBePositive} = 100  % [kg] Mass of AUV
            batteryCapacity (1,1) {mustBeNumeric, mustBePositive} = 1000  % [Wh] Total energy storage on board
            missionTime (1,1) {mustBeNumeric, mustBePositive} = 10  % [hr] Time to complete a mission
            missionBatteryUsed (1,1) {mustBeNumeric, mustBeInRange(missionBatteryUsed, 0, 1, 'exclude-lower')} = 0.8  % [0.XX] percent of battery used during mission
            hotelLoad (1,1) {mustBeNumeric, mustBeNonnegative} = 90  % [W] Baseline power usage of AUV
            chargeRate (1,1) {mustBeNumeric, mustBePositive} = 160  % [Wh] Rate of battery charge
            chargeMethod (1,1) {mustBeNumeric, mustBeMember(chargeMethod, [1,0])} = 1  % Wired (1) or Wireless (0) charging
            batteryRechargeThreshold (1,1) {mustBeNumeric, mustBePositive} = 0.20  % [0.XX] Percentage of battery to initiate recharge
            n_battery (1,1) {mustBeNumeric,  mustBeInRange(n_battery, 0, 1, 'exclude-lower')} = 0.90  % [0.XX] Battery efficiency
            n_powerTransfer (1,1) {mustBeNumeric,  mustBeInRange(n_powerTransfer, 0, 1, 'exclude-lower')} = 0.90  % [0.XX] Power transfer efficiency
            end

            auv.model = model;
            auv.mass = mass; 
            auv.batteryCapacity = batteryCapacity;
            % auv.missionSpecs = [1, missionTime, missionBatteryUsed]; 
            auv.missionTime = missionTime;
            auv.missionBattPrcnt = missionBatteryUsed; 
            auv.hotelLoad = hotelLoad; 
            auv.chargeRate = chargeRate;
            auv.chargeMethod = chargeMethod;
            auv.rechargeThreshold = batteryRechargeThreshold;
            auv.n_battery = n_battery;
            auv.n_powerTransfer = n_powerTransfer;

            % auv.mission = 1;  % Current AUVs only have one mission option, so always will be running 'mission 1'

            % Dependent Constants
            rateLinPowerTransfer = auv.chargeRate;  % Rate of power transfer in linear region
            auv.chargeTime = ( 0.8-(1-auv.missionBattPrcnt) )*auv.batteryCapacity/rateLinPowerTransfer  + 0.4*auv.batteryCapacity./rateLinPowerTransfer;  % Time to charge up to 80% + time to charge from 80% to 100% 
            auv.chargeLoad = (auv.batteryCapacity*auv.missionBattPrcnt + auv.hotelLoad*auv.chargeTime) / auv.n_powerTransfer / auv.n_battery;  % Load on WEC battery to charge AUV battery

            % Initialize battery levels & operational state
            auv.opState = [3, 0];  % AUV initialized to be docked with a full battery at time t = 0
            auv.battery = [auv.batteryCapacity; auv.batteryCapacity];  % Current & previous battery states initialized as full
            auv.wecBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;  % Drawing enough from WEC to stay at full charge
            auv.opTimeComplete = NaN; 
            auv.batteryTime = 0;

        end  % constructor function


        %% Re-Initialize auv Object
        function reset(auv)
            % Re-initializes auv variable properties
            auv.opState = [3, 0];  % AUV initialized to be docked with a full battery at time t = 0
            auv.battery = [auv.batteryCapacity; auv.batteryCapacity];  % Current & previous battery states initialized as full
            auv.wecBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;  % Drawing enough from WEC to stay at full charge
            auv.opTimeComplete = NaN; 
            auv.batteryTime = [];
        end


        %% Battery Tracker
        function calcBatteryLvl(auv, simTime, presetDrawRate)
            % Calculates AUV battery level at a given time (simTime)
            % assuming no changes in operational state since the last
            % battery calculation
            %
            % INPUTS: 
            % simTime - [hr] Simulation time scalar 
            % 
            % OUTPUTS: 
            % auv.battery - [Wh] Vector containing battery level at the end 
            %   of the last operational state (battery(1)), and the current
            %   battery level (battery(2)). battery(1) corresponds to
            %   the timestamp saved as opState(2), and battery(2) 
            %   corresponds to the timestep saved as wec.batteryTime
            % auv.batteryTime - [hr] timestep of most recent battery level
            %   saved as battery(2)
            % presetDrawRate - [W] AUV draw rate from central battery
            %   during recharge. Includes AUV efficiencies and hotel load.
            %   Overrides the set AUV recharge rate, used to slow recharge.

            arguments
                auv                             % AUV class object
                simTime (1,1) {mustBeNumeric}   % Simulation time
                presetDrawRate {mustBeNumeric} = []  % Draw rate to override auv's rate. Includes hotel load and auv efficiencies (battery & power transfer)
            end

            % Check draw-rate input
            if ~isempty(presetDrawRate)  % use given draw rate
                if presetDrawRate < 0 || ~isscalar(presetDrawRate)
                    error('Invalid AUV draw rate.')
                end
            end
            
            stateTime = simTime - auv.opState(2); 

            % Evaluate battery level
            if ~isequal(simTime, auv.batteryTime)  % If already calculated the battery level this timestep, don't calculate it again...
                switch auv.opState(1)
                    case 1  % AUV executing mission
                        ratePowerUse = auv.missionBattPrcnt*auv.batteryCapacity / auv.missionTime;  % rate = Battery Wh used / time to use energy
                        auv.battery(2) = auv.battery(1) - ratePowerUse*stateTime;
                    
                    case 2  % AUV docked & recharging
                        if ~isempty(presetDrawRate)
                            % Calculate recharge rate from draw override
                            newLinRecharge = presetDrawRate*auv.n_powerTransfer*auv.n_battery - auv.hotelLoad;
                            if newLinRecharge < 0
                                error('Nonphysical recharge rate calculated for the %s auv.', auv.model);
                            end

                            % Calculate battery level using given draw rate & recalculate opTimeComplete
                            dt = simTime - auv.batteryTime;  % timestep between this and the previous battery calculation (executed each timestep, so this is equal to the simulation timestep)
                            auv.battery(2) = auv.battery(2) + newLinRecharge*dt;
    
                            % Calculate new draw rate
                            remainingRecharge = auv.batteryCapacity - auv.battery(2);
    
                            if remainingRecharge < 0  % pulled more energy from central battery than needed for AUV recharge, oops
                                warning("%s Wh of energy lost due to an attempted AUV overcharge. Central battery c-rate is too low for regular operation of this system.", num2str(abs(remainingRecharge)))
                                auv.battery(2) = auv.batteryCapacity;
    
                            elseif remainingRecharge == 0  % pulled just enough energy to become fully-charged yay
    
                            elseif remainingRecharge >= 0.2*auv.batteryCapacity  % still in the linear recharge region
                                t_80 = ((0.8*auv.batteryCapacity) - auv.battery(2)) / auv.chargeRate; 
                                t_full = t_80 + (0.4*auv.batteryCapacity / auv.chargeRate);
        
                                auv.opTimeComplete = simTime + t_full; 
    
                            else  % remainingRecharge < 0.2*auv.batteryCapacity so in nonlinear charging slowdown region
                                t_current = (0.4*auv.batteryCapacity / auv.chargeRate)*(1 - sqrt( 1+ (4*auv.batteryCapacity - 5*auv.battery(2))/auv.batteryCapacity ) );  
                                t_full = (0.4*auv.batteryCapacity / auv.chargeRate);  % not adding t_80 because built into battery charge equation (> 80%) t_80 coincides with 0
    
                                auv.opTimeComplete = simTime + (t_full - t_current);
                            end
    
                        else
                            % Calculate battery level using auv charge-curve
                            tempBattery = NaN;  
                            rateLinPowerTransfer = auv.chargeRate;  % rate of charge in linear region (up to 80%) (no efficiencies applied yet)
                            
                            % Linear region
                            if auv.battery(2) < (0.8*auv.batteryCapacity) 
                                tempBattery = auv.battery(2) + rateLinPowerTransfer*(simTime - auv.batteryTime);
                            end
                            
                            % Charge > 80% region
                            if (0.8*auv.batteryCapacity) <= auv.battery(2) || (0.8*auv.batteryCapacity) < tempBattery    % last battery calc'd is > 80% or temp battery just calc'd is > 80%
        
                                if (0.8*auv.batteryCapacity) < tempBattery
                                    t_80 = (0.8*auv.batteryCapacity - auv.battery(2))/rateLinPowerTransfer + auv.batteryTime;  % 80% charge in simulation time
                                    t_newBatteryCalc = simTime - t_80;  % Time for calculation in >80% charge curve (relative to 80% Batt at t = 0)
        
                                else
                                    % Calc. time corresponding to the last-calculated battery level (t_current)
                                    t_current = (0.4*auv.batteryCapacity / rateLinPowerTransfer)*(1 - sqrt( 1+ (4*auv.batteryCapacity - 5*auv.battery(2))/auv.batteryCapacity ) );  
                                    dt = simTime - auv.batteryTime; 
                                    if dt < 0 
                                        error('dt outside of expected range during auv battery calculation')
                                    end
            
                                    t_newBatteryCalc = t_current + dt;
                                    
                                end
                                t_100 = 0.4*auv.batteryCapacity/rateLinPowerTransfer;
        
                                if t_newBatteryCalc >= t_100
                                    auv.battery(2) = auv.batteryCapacity;
        
                                else
                                    auv.battery(2) = rateLinPowerTransfer*( t_newBatteryCalc - (rateLinPowerTransfer*t_newBatteryCalc^2 / (0.8*auv.batteryCapacity)) ) + 0.8*auv.batteryCapacity;
    
                                end
        
                            else  % Battery calculation stays in linear region
                                auv.battery(2) = tempBattery; 
        
                            end
                        end
    
                    case 3  % AUV battery full or central battery in low power mode (recharging itself while supporting auv hotel needs)
    
                    otherwise
                        error("Unknown AUV operational state.")
                
                end  % operational cases
            end

            % Calculate Draw Rate
            switch auv.opState(1)
                case 1  % AUV executing mission
                    newBatteryDraw = 0; 
                
                case 2  % AUV docked & recharging
                    if ~isempty(presetDrawRate)
                        if remainingRecharge < 0  % pulled more energy from central battery than needed for AUV recharge, oops
                            newBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;

                        elseif remainingRecharge == 0  % pulled just enough energy to become fully-charged yay
                            newBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;

                        elseif remainingRecharge >= 0.2*auv.batteryCapacity  % still in the linear recharge region
                            newBatteryDraw = (auv.chargeRate + auv.hotelLoad) / auv.n_powerTransfer /auv.n_battery * ones(length(stateTime), 1);  % [W]

                        else  % remainingRecharge < 0.2*auv.batteryCapacity so in nonlinear charging slowdown region
                            t_newBatteryCalc = t_current + dt;
                            newBatteryDraw = (auv.chargeRate.*( 1 - auv.chargeRate*t_newBatteryCalc/(0.4*auv.batteryCapacity) ) + auv.hotelLoad) / auv.n_powerTransfer /auv.n_battery;  % [W] (d/dt of above + hotel load)

                        end

                    else  
                        if isequal(simTime, auv.batteryTime)
                            tempBattery = NaN;  
                            rateLinPowerTransfer = auv.chargeRate;  % rate of charge in linear region (up to 80%) (no efficiencies applied yet)
                            
                            % Linear region
                            if auv.battery(2) < (0.8*auv.batteryCapacity) 
                                tempBattery = auv.battery(2) + rateLinPowerTransfer*(simTime - auv.batteryTime);
                            end

                        end
                        
                        % Charge > 80% region
                        if (0.8*auv.batteryCapacity) <= auv.battery(2) || (0.8*auv.batteryCapacity) < tempBattery    % last battery calc'd is > 80% or temp battery just calc'd is > 80%
                            if (0.8*auv.batteryCapacity) < tempBattery
                                t_80 = (0.8*auv.batteryCapacity - auv.battery(2))/rateLinPowerTransfer + auv.batteryTime;  % 80% charge in simulation time
                                t_newBatteryCalc = simTime - t_80;  % Time for calculation in >80% charge curve (relative to 80% Batt at t = 0)
    
                            else
                                % Calc. time corresponding to the last-calculated battery level (t_current)
                                t_current = (0.4*auv.batteryCapacity / rateLinPowerTransfer)*(1 - sqrt( 1+ (4*auv.batteryCapacity - 5*auv.battery(2))/auv.batteryCapacity ) );  
                                dt = simTime - auv.batteryTime; 
                                if dt < 0 
                                    error('dt outside of expected range during auv battery calculation')
                                end
        
                                t_newBatteryCalc = t_current + dt;
                                
                            end
                            t_100 = 0.4*auv.batteryCapacity/rateLinPowerTransfer;
                            
                            if t_newBatteryCalc >= t_100
                                newBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;
    
                            else
                                newBatteryDraw = (rateLinPowerTransfer.*( 1 - rateLinPowerTransfer*t_newBatteryCalc/(0.4*auv.batteryCapacity) ) + auv.hotelLoad) / auv.n_powerTransfer /auv.n_battery;  % [W] (d/dt of above + hotel load)

                            end
    
                        else  % Battery calculation stays in linear region
                            newBatteryDraw = (rateLinPowerTransfer + auv.hotelLoad) / auv.n_powerTransfer /auv.n_battery * ones(length(stateTime), 1);  % [W]
    
                        end
                    end

                case 3  % AUV battery full or central battery in low power mode (recharging itself while supporting auv hotel needs)
                    newBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;  % Drawing enough from central battery to supply hotel load 
            
            end  % draw rate calc
            
            % Update AUV object properties
            auv.batteryTime = simTime; 

            if abs(newBatteryDraw - auv.wecBatteryDraw) > 0.01
                auv.wecBatteryDraw = newBatteryDraw; 
            end
            
        end  % fn calcBatteryLvl


        %% Change Operational State
        function changeOpState(auv, newOpState, simTime)
            % Manages and tracks values given a change in operational state
            % Changes logged operational state of auv, updates saved 
            % battery levels, resets stateTime, and updates batteryDraw. 
            % 
            % INPUTS: 
            % newOpState - New operational state - 1) Executing AUV 
            %   mission, 2) AUV recharging, 3) AUV docked & fully charged, 
            %   4) Final 10% AUV recharge
            % simTime - [hr] Simulation time at which operational state  
            %   change occurs (scalar)
            
            % Update to new operational state
            auv.opState = [newOpState, simTime];
            
            % Set battery(1) to battery at end of last operational state. NOTE: auv.battery(1) corresponds to auv.opState(2) timestamp
            auv.battery(1) = auv.battery(2); 

            % Calculate time operation will be complete
            rateLinPowerTransfer = auv.chargeRate;  % rate of charge in linear region (up to 80%)

            switch auv.opState(1)
                case 1  % on mission
                    auv.opTimeComplete = simTime + auv.missionTime; 
                
                case 2  % recharging
                    if auv.battery(2) < 0.8*auv.batteryCapacity
                        t_80 = ((0.8*auv.batteryCapacity) - auv.battery(2)) / rateLinPowerTransfer; 
                        t_full = t_80 + (0.4*auv.batteryCapacity / rateLinPowerTransfer);

                        auv.opTimeComplete = simTime + t_full; 
                    
                    else  % auv.battery(2) >= 0.8*auv.maxBattery
                        t_full = (0.4*auv.batteryCapacity / rateLinPowerTransfer);  % not adding t_80 because built into battery charge equation (> 80%) t_80 coincides with 0
                        t_current = (0.4*auv.batteryCapacity / rateLinPowerTransfer)*(1 - sqrt( 1+ (4*auv.batteryCapacity - 5*auv.battery(2))/auv.batteryCapacity ) );
                        
                        if t_full < t_current
                            error('Problem with opTimeComplete for charging case')
                        end

                        auv.opTimeComplete = simTime + (t_full - t_current);

                    end
                
                case 3  % AUV docked but not charging
                    auv.opTimeComplete = NaN;  % Dependant on WEC

            end
            auv.calcBatteryLvl(simTime);  % To set wec battery draws

        end  % fn changeOpState

    end  % Instance Methods

end  % class def