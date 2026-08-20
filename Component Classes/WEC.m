% WEC defines the properties and methods (functions) saved in WEC
% objects. These objects store information on the Wave Energy
% Converter's properties (model, etc...) and can compute values
% relating to their energy generation and battery level.
%
% To create a WEC object named 'wec', use the following syntax.
% wec = WEC(model, characteristicDimension, batteryCapacity, hotelLoad, n_hydro, n_gen, n_battery);
%
% Ama Hartman


classdef WEC < handle 

    % Constants
    properties (GetAccess = public, SetAccess = private)
        model           % Model of WEC
        charDim         % [m] Characteristic dimension (sqrt(4*A/pi) w/ A = device's max horizontal cross-sectional area)
        batteryCapacity % [Wh] Total energy storage onboard
        hotelLoad       % [W] Baseline power usage of WEC
        n_hydro         % [0.XX] Hydrodynamic efficiency (used in calcPowerGen if model input is wave resource information)
        n_gen           % [0.XX] Generator efficiency (same as above)
        n_battery       % [0.XX] Battery efficiency (losses during charge & discharge)
    end

    % Variables
    properties (Access = public)
        battery         % [Wh] Value of current battery level
        powerGenMeans   % [W] Rate of energy generation given wave profile at simulation timesteps (mean values of power at data timesteps)
        meanPowerGen    % [W] Mean value of power generation 
        lowPowerGen     % [W] Mean of lowest generation window (window sized to accomodate AUV charging) OR 0.75% of mean for constant resource inputs
    end
    
    properties (Dependent)
        lowBatteryLvl   % [Wh] Battery level that triggers 'low power mode' (stops charging AUV(s), only supports hotel loads)
    end


    %% Instance Methods (need object as an input)
    methods

        %% Constructor: Creates & returns an object
        function wec = WEC(model, characteristicDimension, batteryCapacity, hotelLoad, n_hydro, n_gen, n_battery)
            % WEC object constructor generates an object with the given
            % properties. Default values are used if no inputs are given.
            arguments
                model (1,1) string = 'Default WEC'  % String containing model name
                characteristicDimension (1,1) {mustBeNumeric, mustBePositive} = 1.5  % [m] Characteristic dimension of WEC
                batteryCapacity (1,1) {mustBeNumeric, mustBePositive} = 500  % [Wh] Total energy storage onboard
                hotelLoad (1,1) {mustBeNumeric, mustBeNonnegative} = 50  % [W] Baseline power usage
                n_hydro (1,1) {mustBeNumeric, mustBeInRange(n_hydro, 0, 1, 'exclude-lower')} = []  % [0.XX] Hydrodynamic efficiency
                n_gen (1,1) {mustBeNumeric, mustBeInRange(n_gen, 0, 1, 'exclude-lower')} = 0.8  % [0.XX] Generator efficiency
                n_battery   (1,1) {mustBeNumeric, mustBeInRange(n_battery, 0, 1, 'exclude-lower')} = 0.9  % [0.XX] Battery efficiency                
            end

            % Assign user-input and/or default values
            wec.model = model;
            wec.charDim = characteristicDimension;
            wec.batteryCapacity = batteryCapacity;
            wec.hotelLoad = hotelLoad; 
            wec.n_hydro = n_hydro;
            wec.n_gen = n_gen;
            wec.n_battery = n_battery;

            if isempty(wec.n_hydro)
                wec.n_hydro = (1.3*wec.charDim + 5.6)/100;  % Babarit, A. 2015 from Driscol '19
            end

            % Initialize WEC battery
            wec.battery = wec.batteryCapacity;  % Initialize at 100% charge

        end % constructor fn


        %% Dependent property calcs..
        function lowBatteryLvl = get.lowBatteryLvl(wec)
            % Sets the low battery level to 5% of the maximum capacity
            lowBatteryLvl = 0.05*wec.batteryCapacity;  
        end

        
        %% Re-Initialize wec object
        function reset(wec)
            wec.battery = wec.batteryCapacity;
        end

        
        % %% Modify hotel load - retired Aug 2026
        % function modifyHotelLoad(wec, newHotelLoad)
        %     % re-writes the hotel load of the wec object as the given value
        %     wec.hotelLoad = newHotelLoad;
        % end

        
        % %% Calculate 'low power' generation threshold - retired Aug 2026
        % function calcLowPower(wec, resourceDataType, dt, auv)
        %     switch resourceDataType
        %         case{1, 2, 4, 6}
        %             chargeWindow = ceil([auv.chargeTime]./dt);
        %             if numel(chargeWindow) > 1
        %                 chargeWindow = mean(chargeWindow);
        %             end
        %             tmpPwrGen = movmean(wec.powerGenMeans, chargeWindow);
        %             wec.lowPowerGen = min(tmpPwrGen(floor(0.5*chargeWindow) : end - floor(0.5*chargeWindow)) );  % Exclude truncated means
        % 
        %         otherwise
        %             wec.lowPowerGen = 0.75*wec.meanPowerGen;
        %     end
        % 
        % end  % low power threshold fn
        

        %% Energy Generation
        function reshapePowerGen(wec, powerData, timeData, simHrs, dt)
            % Reshapes power generation as a function of time given power 
            % timeseries and sea state. Reshapes power generation to fit 
            % simulation timesteps.
            %
            % Inputs: 
            % - powerData: [W] nx1 vector of power as a function of time
            % - timeData: [hr] nx1 vector of time corresponding to powerData
            % - simHrs: [hr] Simulation hours
            % - dt: [hr] Simulation timestep

            % Time series calculations
            powerData_dt = mean(diff(timeData));  

            % Resample data to fit simulation timestep (dt)
            if powerData_dt ~= dt
                timeDataResampled = (dt:dt:timeData(end))';
                powerDataResampled = interp1(timeData, powerData, timeDataResampled, 'linear', powerData(1)) * wec.n_gen;  % Fills NaN values at beginning (between t = 0 and t = timeData(1)) with powerData(1)
            else
                timeDataResampled = timeData;
                powerDataResampled = powerData * wec.n_gen;
            end

            % Resize time series to fit simulation length (simHrs)
            if timeDataResampled(end) < simHrs
                % repeat data 
                timeDataResized = (dt:dt:simHrs)';  % new timestep matches simulation (dt)
                powerDataResized = powerDataResampled(mod(0:length(timeDataResized)-1, length(powerDataResampled)) +1);

            elseif timeDataResampled(end) > simHrs 
                % Need to trim data to simulation length
                trimIndx = find(timeDataResampled <= simHrs, 1, 'last');
                powerDataResized = powerDataResampled(1:trimIndx);

            else
                % timeData == simHrs --> no resizing needed
                powerDataResized = powerDataResampled;
            end

            wec.powerGenMeans = powerDataResized; 
            wec.meanPowerGen = mean(powerDataResized);

        end  % power gen fn


        %% Energy Generation from Wave Specs
        function calcPowerGen(wec, waveData, dataWindowSelection, simTime, plotPwrGen, windowOverrideIndx)
            % Calculates power generation as a function of simulation time
            % gven wave profile characteristics
            % 
            % INPUTS: 
            % - waveData: Struct including 'sigWaveHeight',
            %   'waveEnergyPeriod', 'peakPeriod', and 'dataTime' vectors
            % - dataWindowSelection: either 'minPwr', 'maxPwr', or
            %   'meanPwr' indicating whether the user wants to pull the
            %   min/max/mean power generation window.
            % - simTime: [hr] Time series vector used for simulation
            % - WindowOverrideIndx: Index of wave data time series to use
            %   for the simulation starting point, 0 to use window selection
            %   based on min/max/mean power generation window. Set to 403 and
            %   use Oregon (US West Coast) dataset to replicate JOE paper

            % Constants
            rho = 1025;  % [kg/m^3] Seawater Density
            g = 9.81;  % [m/s^2] Gravitational acceleration
            
            % Power Calculations
            wavePwrDensity = rho * g^2 * waveData.sigWaveHeight.^2 .* waveData.waveEnergyPeriod / (64*pi);  % [W/m]
            pwrGen = wavePwrDensity * wec.charDim * wec.n_hydro * wec.n_gen; 
            budalLimit = rho * g * pi^2 * (wec.charDim/2)^3 * waveData.sigWaveHeight ./ (6 * waveData.peakPeriod); 

            wecPwrGen = min(pwrGen, budalLimit);

            if ~isscalar(wecPwrGen)  

                dataDt = waveData.dataTime(2) - waveData.dataTime(1); 
    
                % If data shorter than simulation, need to repeat.. 
                if waveData.dataTime(end) < simTime(end)  % Data window < simulation window & need to do some repeats. This still needs to be tested.
                    dataTimeSeries = dataDt:dataDt:simTime(end);
                    wecPowerGenWindowed = wecPwrGen(mod(0:length(dataTimeSeries)-1, length(wecPwrGen)) +1);  % repeating power for length of simulation

                elseif windowOverrideIndx == 0
                    % Pull time period for simulation
                    window = ceil(simTime(end)/dataDt);
                    tmpMeans = movmean(wecPwrGen, window);
                    tmpMeans(1:floor(0.5*window)) = NaN;  % Cut out truncated windows
                    tmpMeans(end-floor(0.5*window):end) = NaN;
                    dataTimeWindowed = waveData.dataTime(1:window+1);
        
                    switch dataWindowSelection
                        case 'minPwr'
                            [~, indx] = min(tmpMeans); 
                        case 'maxPwr'
                            [~, indx] = max(tmpMeans); 
                        case 'meanPwr'
                            [~, indx] = min(abs(mean(tmpMeans,'omitmissing') - tmpMeans)); 
                    end
        
                    if rem(window, 2) == 0  % (if window is divisible by 2)
                        wecPowerGenWindowed = wecPwrGen(indx - (0.5*window)+1 : indx + (0.5*window) +1 );
    
                    else
                        wecPowerGenWindowed = wecPwrGen(indx - (floor(0.5*window)) : indx + floor(0.5*window)+1);
    
                    end

                else  % windowOverrideIndx ~= 0  % Pull time period for simulation as defined by given override index
                    window = ceil(simTime(end)/dataDt);
                    dataTimeWindowed = waveData.dataTime(1:window+1);
                    wecPowerGenWindowed = wecPwrGen(windowOverrideIndx : windowOverrideIndx + window);

                end
    
                simDt = simTime(2) - simTime(1);

                if dataDt ~= simDt
                    % resample to match simulation timestep
                    wec.powerGenMeans = interp1((dataTimeWindowed-dataDt), wecPowerGenWindowed, simTime, 'linear', wecPowerGenWindowed(1));
                else
                    % no resampling needed
                    wec.powerGenMeans = wecPowerGenWindowed;
                end

                % Save mean value
                wec.meanPowerGen = mean(wec.powerGenMeans); 

                if plotPwrGen == 1
                    %% Plot Power Generation during simulation window 
                    set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

                    if windowOverrideIndx == 0
                        tempDatetimes = interp1(dataTimeWindowed, dataDateTimes(indx-(0.5*window)+1 : indx + (0.5*window)+1), simTime);  
                    else
                        tempDatetimes = interp1(dataTimeWindowed, dataDateTimes(windowOverrideIndx : windowOverrideIndx + window), simTime);  
                    end
                    figure; plot(tempDatetimes, wec.powerGenMeans)
                    grid on; 
                    ylabel('Power Generated [W]'); 
                    xlabel('Simulation Window'); 
                    sgtitle(''); 
                    hold on; plot(tempDatetimes, wec.meanPowerGen * ones(size(tempDatetimes)));
                end

            else 
                wec.meanPowerGen = wecPwrGen;
                wec.powerGenMeans = wecPwrGen * ones(size(simTime));
                
            end

        end  % calc power gen fn

    end  % Instance Methods
    
end  % class def
