% Ama Hartman

%% Define Class
classdef AUV < handle   
    % AUV defines the properties and methods (functions) saved in AUV objects.
    % These objects store information on the Unmanned Underwater Vehicle's
    % properties (model ID, mass, energy storage, endurance, charging rate,
    % etc...) and can compute values relating to their battery level and
    % mission duration. 
    %
    % To create a AUV object named 'auv', use the following syntax. 
    % auv = AUV(model);
    % 
    % Currently supported models 'A' - 'U' are generic examples

    properties
        model           % Model of auv (i.e. 'Iver3')
        mass            % [kg] Mass of auv 
        mission         % Mission number currently using / running
        missionSpecs    % Column-vector matrix containing mission number (1,2..), time to complete [hr], and percentage of battery used (0.XX)
        maxBattery      % [Wh] Total energy storage onboard 
        chargeRate      % [W] Rate of battery charge 
        chargeTime      % [hr] Column-vector of time to go to 100% charge after respective mission. NOTE: Mission must use >= 20% of battery for this calc to be accurate. 
        chargeMethod    % Wired (1) or Wireless (2) charging
        chargeLoad      % [Wh] Column-vector of total load put on WEC to charge auv after a given mission, assuming AUV started with a full battery
        hotelLoad       % [W] Baseline power usage of AUV
        n_powerTransfer % [0.XX] Power transfer efficiency 
        n_battery       % [0.XX] Battery efficiency (losses during charge & discharge)
        opState         % Vector containing current operational state and simulation timestamp of respective state change in hours ([opState, t]) with 1) Executing AUV mission, 2) AUV recharging, 3) AUV docked & fully charged
        opTimeComplete  % [hr] Time current opState will be complete
        battery         % Vector containing battery level at the end of the last operational state (battery(1)) and the current battery level (battery(2)) [Wh] ([B(i-1); B(i)]). battery(1) corresponds to timestamp saved as opState(2)
        batteryTime     % [hr] Time corresponding to the current battery level (battery(2)) - NOTE: Not currently used externally as of 8/14/25    
        wecBatteryDraw  % [W] Time series of energy draw level (0 if not charging) with wecBatteryDraw(end) corresponding with current draw level
        rechargeThreshold % [0.XX] Percentage of battery to recharge at. Default is 20%

    end


    %% Instance Methods (need object as an input)
    methods

        %% Constructor: Creates & returns an object
        function auv = AUV(model)
            % Creates an object given AUV model

            auv.rechargeThreshold = 0.20;  
            
            auv.mission = 1;  % Current auv's only have one mission option, so always will be running 'mission 1

            
            switch model 
                case 'A'  % Source: Driscol '19
                    auv.model = model;  
                    auv.mass = 27;  % Iver3 model options range from 27-38.5 kg
                    auv.missionSpecs = [1, 8*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];  % Iver3 model options range from 8-14 hrs until recharge needed
                    auv.maxBattery = 800; 
                    auv.chargeRate = 160;
                    auv.chargeMethod = 1;  
                    auv.hotelLoad = 90;  

                case 'B'  % Source: Driscol '19
                    auv.model = model;
                    auv.mass = 38.5;  % range from 27-38.5 kg
                    auv.missionSpecs = [1, 14*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];  % range from 8-14 hrs
                    auv.maxBattery = 800; 
                    auv.chargeRate = 160;
                    auv.chargeMethod = 1;
                    auv.hotelLoad = 90;  

                case 'C'  % *all UUVs* Website source: https://hii.com/what-we-do/capabilities/unmanned-systems/remus-uuvs/
                    auv.model = model;
                    auv.mass = 38.6;
                    auv.missionSpecs = [1, 10*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 1500; 
                    auv.chargeRate = auv.maxBattery / 6; 
                    auv.chargeMethod = 1;
                    auv.hotelLoad = 90;  

                case 'D'  % 4
                    auv.model = model;
                    auv.mass = 58.5;
                    auv.missionSpecs = [1, 20*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 3000; 
                    auv.chargeRate = auv.maxBattery / 6;
                    auv.chargeMethod = 1;

                case 'E'  % 5
                    auv.model = model;
                    auv.mass = 70.3;
                    auv.missionSpecs = [1, 30*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 4500; 
                    auv.chargeRate = auv.maxBattery / 6;
                    auv.chargeMethod = 1;

                case 'F'  % 6
                    auv.model = model;
                    auv.mass = 210;
                    auv.missionSpecs = [1, 42*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 9600; 
                    auv.chargeRate = auv.maxBattery / 8;
                    auv.chargeMethod = 1;

                case 'G'  % 7
                    auv.model = model;
                    auv.mass = 279;
                    auv.missionSpecs = [1, 80*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 19300; 
                    auv.chargeRate = auv.maxBattery / 10;
                    auv.chargeMethod = 1;

                case 'H'
                    auv.model = model;
                    auv.mass = 347;
                    auv.missionSpecs = [1, 110*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 28900; 
                    auv.chargeRate = auv.maxBattery / 12;
                    auv.chargeMethod = 1;

                case 'I'  %9
                    auv.model = model;
                    auv.mass = 1630;
                    auv.missionSpecs = [1, 25*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 17550; 
                    auv.chargeRate = auv.maxBattery / 24;
                    auv.chargeMethod = 1;

                case 'J'
                    auv.model = model;
                    auv.mass = 70;
                    auv.missionSpecs = [1, 8*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 1900;
                    auv.chargeRate = auv.maxBattery / 6;
                    auv.chargeMethod = 1;

                case 'K'
                    auv.model = model;
                    auv.mass = 250;
                    auv.missionSpecs = [1, 36*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 7600;
                    auv.chargeRate = auv.maxBattery / 6;
                    auv.chargeMethod = 1;

                case 'L'  % 12
                    auv.model = model;
                    auv.mass = 750;
                    auv.missionSpecs = [1, 25*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 13500;
                    auv.chargeRate = auv.maxBattery / 6;
                    auv.chargeMethod = 1;  

                case 'M'
                    auv.model = model;
                    auv.mass = 72.6;
                    auv.missionSpecs = [1, 3.5*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 1500;
                    auv.chargeRate = auv.maxBattery / 6;
                    auv.chargeMethod = 1;

                case 'N'  % 14
                    auv.model = model;
                    auv.mass = 2200;
                    auv.missionSpecs = [1, 72*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 62500;
                    auv.chargeRate = auv.maxBattery / 8;
                    auv.chargeMethod = 1;

                case 'O'
                    auv.model = model;
                    auv.mass = 8000;
                    auv.missionSpecs = [1, 360*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 400000;
                    auv.chargeRate = auv.maxBattery / 8;
                    auv.chargeMethod = 1;

                case 'P'  % 16
                    auv.model = model;
                    auv.mass = 1000;
                    auv.missionSpecs = [1, 24*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 24000;
                    auv.chargeRate = auv.maxBattery / 8;
                    auv.chargeMethod = 1;

                case 'Q'
                    auv.model = model;
                    auv.mass = 1550;
                    auv.missionSpecs = [1, 74*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 48000;
                    auv.chargeRate = auv.maxBattery / 8;
                    auv.chargeMethod = 1;

                case 'R'  % 18
                    auv.model = model;
                    auv.mass = 25;
                    auv.missionSpecs = [1, 10*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 600;
                    auv.chargeRate = auv.maxBattery / 5;
                    auv.chargeMethod = 1;
                    auv.hotelLoad = 1;  % [W] from email correspondance with boxfish

                case 'S'
                    auv.model = model;
                    auv.mass = 28;
                    auv.missionSpecs = [1, 10*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 600;
                    auv.chargeRate = auv.maxBattery / 4;
                    auv.chargeMethod = 2;

                case 'T'  %20
                    auv.model = model; 
                    auv.mass = 650;
                    auv.missionSpecs = [1, 10.8*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 12000; 
                    auv.chargeRate = auv.maxBattery / 3.64;
                    auv.chargeMethod = 1;

                case 'U'
                    auv.model = model; 
                    auv.mass = 1300;
                    auv.missionSpecs = [1, 21.6*(1-auv.rechargeThreshold), (1-auv.rechargeThreshold)];
                    auv.maxBattery = 30000; 
                    auv.chargeRate = auv.maxBattery / 9.09;
                    auv.chargeMethod = 1;

                otherwise
                    error('AUV model %s not yet supported. Update AUV class file to support this model.', model);

            end

            % Hotel Load (extrapolate from Boxfish AUV)
            interpHotelLoad = 1;  % 1 - yes, interpolate
            
            if isempty(auv.hotelLoad)
                if interpHotelLoad == 1 
                    auv.hotelLoad = round(auv.mass / 25);  % Interpolates based on auv mass relative to hotel load of Boxfish AUV
                
                else
                    auv.hotelLoad = 90;  % [W] Generic, used in 'Wave-Powered AUV Recharging: A Feasibility Study' by B. P. Driscol, A. Gish, and R. G. Coe in 2019

                end
            end
            

            % Efficiencies
            auv.n_battery = 0.90;  % from 'A unified simulation framework for wave energy powered underwater vehicle docking and charging' by M. Chen Et. Al. in 2024

            if auv.chargeMethod == 1
                auv.n_powerTransfer = 0.9;  % https://www.researchgate.net/publication/264124522_Efficiency_Comparison_of_Wire_and_Wireless_Battery_Charging_Based_on_Connection_Probability_Analysis

            elseif auv.chargeMethod == 2
                auv.n_powerTransfer = 0.50;  % Baseline efficiency for wireless power transfer from 'Adaptive Wireless Power for Subsea Vehicles' by D. Manalang Et. Al. in 2022

            else
                error('Unsupported charge method. Please specify 1 for wired charging or 2 for wireless charging.');

            end


            % Calculate chargeTime (from rechargeThreshold to 100%)
            rateLinPowerTransfer = auv.chargeRate;  % Rate of power transfer in linear region
            auv.chargeTime = ( 0.8-(1-auv.missionSpecs(:,3)) )*auv.maxBattery/rateLinPowerTransfer  + 0.4*auv.maxBattery./rateLinPowerTransfer;  % Time to charge up to 80% + time to charge from 80% to 100% 

            % calculate load on WEC battery to charge AUV battery [Wh]
            auv.chargeLoad = (auv.maxBattery*auv.missionSpecs(:, 3) + auv.hotelLoad*auv.chargeTime) / auv.n_powerTransfer / auv.n_battery;  

            % Initialize battery levels
            auv.opState = [3, 0];  % AUV initialized to be docked with a full battery at time t = 0
            auv.battery = [auv.maxBattery; auv.maxBattery];  % Current & previous battery states initialized as full
            auv.wecBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;  % Drawing enough from WEC to stay at full charge

            % Initialize opTimeComplete
            auv.opTimeComplete = NaN; 

        end  % constructor function


        %% Battery Tracker
        function calcBatteryLvl(auv, simTime)
            % Calculates AUV battery level at a given time (simTime)
            % assuming no changes in operational state since the last
            % battery calculation
            %
            % INPUTS: 
            % simTime - [hr] Simulation time scalar 
            % 
            % OUTPUTS: 
            % wec.battery - [Wh] Vector containing battery level at the end 
            %   of the last operational state (battery(1)), and the current
            %   battery level (battery(2)). battery(1) corresponds to
            %   the timestamp saved as opState(2), and battery(2) 
            %   corresponds to the timestep saved as wec.batteryTime
            % wec.batteryTime - [hr] timestep of most recent battery level
            %   saved as battery(2)


            stateTime = simTime - auv.opState(2); 

            switch auv.opState(1)
                case 1  % AUV executing mission
                    ratePowerUse = auv.missionSpecs(auv.mission, 3)*auv.maxBattery / auv.missionSpecs(auv.mission, 2);  % rate = Battery Wh used / time to use energy
                    auv.battery(2) = auv.battery(1) - ratePowerUse*stateTime;
                    auv.batteryTime = simTime;

                    auv.wecBatteryDraw = 0; 
                
                case 2  % AUV docked & recharging
                    tempBattery = NaN;  
                    rateLinPowerTransfer = auv.chargeRate;  % rate of charge in linear region (up to 80%) (no efficiencies applied yet)
                    
                    % Linear region
                    if auv.battery(2) < (0.8*auv.maxBattery) 
                        tempBattery = auv.battery(2) + rateLinPowerTransfer*(simTime - auv.batteryTime);
                    end
                    
                    % Charge > 80% region
                    if (0.8*auv.maxBattery) <= auv.battery(2) || (0.8*auv.maxBattery) < tempBattery    % last battery calc'd is > 80% or temp battery just calc'd is > 80%

                        if (0.8*auv.maxBattery) < tempBattery
                            t_80 = (0.8*auv.maxBattery - auv.battery(2))/rateLinPowerTransfer + auv.batteryTime;  % 80% charge in simulation time
                            t_newBatteryCalc = simTime - t_80;  % Time for calculation in >80% charge curve (relative to 80% Batt at t = 0)

                        else
                            % Calc. time corresponding to the last-calculated battery level (t_current)
                            t_current = (0.4*auv.maxBattery / rateLinPowerTransfer)*(1 - sqrt( 1+ (4*auv.maxBattery - 5*auv.battery(2))/auv.maxBattery ) );  
                            dt = simTime - auv.batteryTime; 
                            if dt < 0 
                                error('dt outside of expected range during auv battery calculation')
                            end
    
                            t_newBatteryCalc = t_current + dt;
                            
                        end
                            t_100 = 0.4*auv.maxBattery/rateLinPowerTransfer;
    
                            if t_newBatteryCalc >= t_100
                                auv.battery(2) = auv.maxBattery;
                                auv.wecBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;
    
                            else
                                auv.battery(2) = rateLinPowerTransfer*( t_newBatteryCalc - (rateLinPowerTransfer*t_newBatteryCalc^2 / (0.8*auv.maxBattery)) ) + 0.8*auv.maxBattery;
                                auv.wecBatteryDraw = (rateLinPowerTransfer.*( 1 - rateLinPowerTransfer*t_newBatteryCalc/(0.4*auv.maxBattery) ) + auv.hotelLoad) / auv.n_powerTransfer /auv.n_battery;  % [W] (d/dt of above + hotel load)

                            end

                    else  % Battery calculation stays in linear region
                        if isnan(tempBattery)
                            error('AUV battery calculation error')
                        end
                        
                        auv.battery(2) = tempBattery; 
                        auv.wecBatteryDraw = (rateLinPowerTransfer + auv.hotelLoad) / auv.n_powerTransfer /auv.n_battery * ones(length(stateTime), 1);  % [W]

                    end

                    auv.batteryTime = simTime;

                case 3  % AUV battery full or WEC in low power mode, WEC recharging itself while supporting auv hotel needs
                    auv.battery(2) = auv.battery(2);  
                    auv.batteryTime = simTime; 
                    auv.wecBatteryDraw = auv.hotelLoad / auv.n_powerTransfer / auv.n_battery;  % Drawing enough from WEC to stay at full charge 
            
            end  % operational cases
            
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
            %   change occurs 

            
            % Update to new operational state
            auv.opState = [newOpState, simTime];
            
            % Set battery(1) to battery at end of last operational state. NOTE: auv.battery(1) corresponds to auv.opState(2) timestamp
            auv.battery(1) = auv.battery(2); 

            % Calculate time operation will be complete
            rateLinPowerTransfer = auv.chargeRate;  % rate of charge in linear region (up to 80%)

            switch auv.opState(1)
                case 1  % on mission
                    auv.opTimeComplete = simTime + auv.missionSpecs(auv.mission, 2); 
                
                case 2  % recharging
                    if auv.battery(2) < 0.8*auv.maxBattery
                        t_80 = ((0.8*auv.maxBattery) - auv.battery(2)) / rateLinPowerTransfer; 
                        t_full = t_80 + (0.4*auv.maxBattery / rateLinPowerTransfer);

                        auv.opTimeComplete = simTime + t_full; 
                    
                    else  % auv.battery(2) >= 0.8*auv.maxBattery
                        t_full = (0.4*auv.maxBattery / rateLinPowerTransfer);  % not adding t_80 because built into battery charge equation (> 80%) t_80 coincides with 0
                        t_current = (0.4*auv.maxBattery / rateLinPowerTransfer)*(1 - sqrt( 1+ (4*auv.maxBattery - 5*auv.battery(2))/auv.maxBattery ) );
                        
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