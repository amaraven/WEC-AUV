% Ama Hartman 

classdef WEC < handle 
    % WEC defines the properties and methods (functions) saved in WEC 
    % objects. These objects store information on the Wave Energy 
    % Converter's properties (model, etc...) and can compute values 
    % relating to their energy generation and battery level. 
    % 
    % To create a WEC object named 'wec', use the following syntax. 
    % wec = WEC(model);
    % 
    % Currently supported models include the following: 
    % 'generic' - A generic placeholder! How exciting.

    properties
        model           % Model of WEC
        charDim         % [m] Characteristic dimension (sqrt(4*A/pi) w/ A = device's max horizontal cross-sectional area)
        maxBattery      % [Wh] Total energy storage onboard
        powerGen        % [W] Rate of energy generation (given wave profile)
        powerGenTime    % [h] Time series corresponding to energy generation
        hotelLoad       % [W] Baseline power usage of WEC
        n_hydro         % [0.XX] Hydrodynamic efficiency (used in calcPowerGen if model input is wave resource information)
        n_gen           % [0.XX] Generator efficiency (same as above)
        n_battery       % [0.XX] Battery efficiency (losses during charge & discharge)
        battery         % [Wh] Vector of battery levels throughout simulation
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
        function wec = WEC(model)
            % Creates an object given the WEC model 

            switch model
                case 'generic'
                    wec.model = model;
                    wec.charDim = 1.5;  
                    wec.hotelLoad = 50;  % generic, colleague's estimate (thanks to Corey Crisp)
                    wec.n_hydro = (1.3*wec.charDim + 5.6)/100;  % Babarit, A. 2015 from Driscol '19
                    wec.n_gen = 0.8;  % Driscol '19,  Chen '24
                    wec.n_battery = 0.90;  % Chen '24
                    wec.maxBattery = 500;  % [Wh] generic small battery to support a few hours of hotel load
                otherwise
                    error('WEC model not yet supported. Update WEC class file to support this model.')
            end

        end % constructor fn


        %% Dependent property calcs..
        function lowBatteryLvl = get.lowBatteryLvl(wec)

            % Sets the low battery level to 5% of the maximum capacity
            lowBatteryLvl = 0.05*wec.maxBattery;  
        
        end
        

        %% Energy Generation
        function reshapePowerGen(wec, powerData, timeData, simHrs, dt)
            % Reshapes power generation as a function of time given power 
            % timeseries and sea state. Reshapes power generation to fit 
            % simulation timesteps.
            %
            % Inputs: 
            % - powerData: [W] nx1 vector of power as a function of time
            % - timeData: [s] nx1 vector of time corresponding to powerData
            % - simHrs: [hr] Simulation hours
            % - dt: [hr] Simulation timestep

            % Time series calculations
            powerData_dt = mean(diff(timeData));  
            simSeconds = powerData_dt:powerData_dt:simHrs*60*60;  % Time series for length of simulation with timestep defined by wave resource data
            wec.powerGenTime = simSeconds'/60/60;  % Time series corresponding to powerGen in hours

            % Create vector of repeating power generation for the length of simulation.
            wec.powerGen = powerData(mod(0:length(simSeconds)-1, length(powerData)) +1) * wec.n_gen;  

            % Reshape power generation according to simulation timesteps
            powerGenReshaped = reshape(wec.powerGen, dt/(wec.powerGenTime(2)-wec.powerGenTime(1)), []);
            wec.powerGenMeans = (mean(powerGenReshaped, 1))';  %  [W]   
            wec.meanPowerGen = mean(wec.powerGenMeans);  % [W]

        end  % power gen fn


        %% Energy Generation from Wave Specs
        function calcPowerGen(wec, waveData, dataWindowSelection, simTime, plotPwrGen, windowOverrideIndx)
            % Calculates power generation as a function of simulation time
            % gven wave profile characteristics
            % 
            % INPUTS: 
            % - waveData: Table containing at least 'SignificantWaveHeight', 'EnergyPeriod' and 'PeakPeriod'
            % - dataWindowSelection: either 'minPwr', 'maxPwr', or
            %   'meanPwr' indicating whether the user wants to pull the
            %   min/max/mean power generation window.
            % - simTime: [hr] Time series vector used for simulation
            % - WindowOverrideIndx: Index of wave data time series to use
            % for the simulation starting point, 0 to use window selection
            % based on min/max/mean power generation window

            % Constants
            rho = 1025;  % [kg/m^3] Seawater Density
            g = 9.81;  % [m/s^2] Gravitational acceleration
            
            % Power Calculations
            wavePwrDensity = rho * g^2 * waveData.SignificantWaveHeight.^2 .* waveData.EnergyPeriod / (64*pi);  % [W/m]
            pwrGen = wavePwrDensity * wec.charDim * wec.n_hydro * wec.n_gen; 
            budalLimit = rho * g * pi^2 * (wec.charDim/2)^3 * waveData.SignificantWaveHeight ./ (6 * waveData.PeakPeriod); 

            wecPwrGen = min(pwrGen, budalLimit);

            if all(ismember({'Year', 'Month', 'Day', 'Hour', 'Minute', 'EnergyPeriod', 'PeakPeriod', 'SignificantWaveHeight'}, waveData.Properties.VariableNames))
                
                % Data time vector (entire year)
                dataDateTimes = datetime(waveData.Year, waveData.Month, waveData.Day, waveData.Hour, waveData.Minute, zeros(size(waveData.Year)));
                dataTime = hours(dataDateTimes-dataDateTimes(1)+ (dataDateTimes(2)-dataDateTimes(1)));

                dataDt = dataTime(2) - dataTime(1); 
    
                % If data shorter than simulation, need to repeat.. 
                if dataTime(end) < simTime(end)  % Data window < simulation window & need to do some repeats. This still needs to be tested.
                    dataTimeSeries = dataDt:dataDt:simTime(end);
                    wecPowerGenWindowed = wecPwrGen(mod(0:length(dataTimeSeries)-1, length(wecPwrGen)) +1);  % repeating power for length of simulation

                elseif windowOverrideIndx == 0
                    % Pull time period for simulation
                    window = ceil(simTime(end)/dataDt);
                    tmpMeans = movmean(wecPwrGen, window);
                    tmpMeans(1:floor(0.5*window)) = NaN;  % Cut out truncated windows
                    tmpMeans(end-floor(0.5*window):end) = NaN;
                    dataTimeWindowed = dataTime(1:window+1);
        
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
                    dataTimeWindowed = dataTime(1:window+1);
                    wecPowerGenWindowed = wecPwrGen(windowOverrideIndx : windowOverrideIndx + window);

                end
    
                if hours(dataDateTimes(2)-dataDateTimes(1)) > 1
                    % Interpolate to up-sample data to simTime 
                    wec.powerGenMeans = interp1((dataTimeWindowed- dataDt), wecPowerGenWindowed, simTime);

                elseif hours(dataDateTimes(2)-dataDateTimes(1)) < 1
                    % Down-sample data to simTime
                    wec.powerGenMeans = reshape(wecPowerGenWindowed, (simTime(2)-simTime(1))/(dataDt), []);
                end

                % Save final power gen. vector
                wec.meanPowerGen = mean(wec.powerGenMeans); 

                if plotPwrGen == 1
                    %% Plot Power Generation during simulation window 
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

            elseif isequal({'SignificantWaveHeight', 'EnergyPeriod', 'PeakPeriod'}, waveData.Properties.VariableNames)  % If input is a table with single values...
                wec.meanPowerGen = wecPwrGen;

            else 
                error('Problem calculating power generation. Unrecognized wave resource data format.')
            end

        end  % calc power gen fn
        

        %% Calculate 'low power' generation threshold
        function calcLowPower(wec, resourceDataType, dt, auv)
            switch resourceDataType
                case{1, 2, 4, 6}
                    chargeWindow = ceil(auv.chargeTime(auv.mission)/dt);
                    tmpPwrGen = movmean(wec.powerGenMeans, chargeWindow);
                    wec.lowPowerGen = min(tmpPwrGen(floor(0.5*chargeWindow) : end - floor(0.5*chargeWindow)) );  % Exclude truncated means
                    
                otherwise
                    wec.lowPowerGen = 0.75*wec.meanPowerGen;
            end

        end  % low power threshold fn

    end  % Instance Methods
end  % class def
