% Ama Hartman

classdef ModelOutput < handle
    % Model Output defines the properties and methods saved in ModelOutput 
    % objects. These objects store information on the model inputs &
    % outputs, and can generate plots from the saved data.
    % 
    % To create a ModelOutput object named modOut, use the following syntax
    % modOut = ModelOutput(simTime, depVar, resourceDataType)
    
    properties
        % Model Input Properties
        simTime             % Time series vector used for simulation    
        depVar              % Dependent variable chosen for simulation
        seaState            % sea state
        meanPowerGen        % [W] Mean value
        powerGenMeans       % [W] Time series of power gen throughout simulation
        resourceDataType    % String describing input type
        dataIn              % Input data specifics
        userDefinedBattery  % 1 - User input max battery / 0 - model outputs battery capacity
        centralBatteryCapacity  % Vector containing total energy storage of central battery (1) [Wh]
        incorpStagger       % 1 - Incorporates even stagger between AUV deployments. 0 - No stagger introduced into the mission scheduling logic
        maxFleetSize        % 0 - No maximum fleet size imposed. <nonzero> - maximum AUV fleet size

        % Model Output Properties
        fleetSize         % (1xm) Array with number of AUVs in fleet for each test case
        auvMissionLength    % [hr] Could be 1x1 or 1xm depending on if mission length changes between test cases
        ratePwrUsed         % [W] Rate AUV uses power (power used during mission + recharge [Wh] / mission + recharge time [h])
        energyStorageBatteryLvl  % [Wh] Battery level of central energy storage as a function of time with rows corresponding to time steps and columns corresponding to test cases
        wecBatteryLvl       % [Wh] Battery level of WEC as a function of time with rows corresponding to time steps and columns corresponding to test cases
        auvBatteryLvl       % [Wh] Battery level of AUV(s) as a function of time. Cell array with columns corresponding to test case. Each array item is a nxm matrix with n = timesteps and m = auv number
        auvSchedule         % Operational state schedule of AUV(s) as a function of time with 1) Executing AUV mission, 2) AUV recharging, 3) AUV docked & fully charged. Organized in the same manner as 'auvBatteryLvl'. 
        auvTimeOnMission    % [h] Time each AUV spends 'on-mission' during the simulation. Cell array with columns corresponding to test case. Each array item is a 1xm matrix with m corresponding to the number of AUVs in fleet.
        auvTimeOnMissionCorrected   % [h] Time each AUV spends 'on-mission' after all AUVs are deployed in the case of a staggered deployment. Empty if deployment stagger not incorporated. Used for performance calculations using an adjusted domain to exclude time when AUVs are artificially held at dock
        auvModels           % Model(s) of auv used for simulation
        auvFleet            % Cell array containing auv fleet for each test case
    end

    properties (Hidden = true)
        fleetNumber         % Legacy variable name for fleetSize, needed to read old data
    end
    
    
    %% Instance Methods (need object as an input)
    methods

        %% Constructor
        function modOut = ModelOutput(simTime, depVar, resourceDataType)
            % MODELOUT Constructs an instance of this class

            modOut.simTime = simTime;
            modOut.depVar = depVar;
            modOut.resourceDataType = resourceDataType;

            switch depVar
                case 'AUV Model'
                    modOut.auvModels = [{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}, {'G'}, {'H'}, {'I'}, {'J'}, {'K'}, {'L'}, {'M'}, {'N'}, {'O'}, {'P'}, {'Q'}, {'R'}, {'S'}, {'T'}, {'U'}];
                
                case 'WEC Power Gen / Wave Resource'
                case 'Battery Specs'
            end

        end  % constructor fn

        %% Plot Battery Tracker
        function plotBatteryTracker(modOut)
            %plotBatterytracker generates a plot of the central and AUV
            %batteries as a function of time during a simulation.
            set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

            % Generate figure, set title
            for figNum = 1:length(modOut.auvBatteryLvl)
                figure; 
                switch modOut.depVar
                    case 'AUV Model'
                        sgtitle(modOut.auvModels{figNum});
                    case 'WEC Power Gen / Wave Resource'
                        switch modOut.resourceDataType
                            case 1
                                sgtitle(strcat('SeaState_', num2str(modOut.seaState(figNum))));
                            case 3
                                sgtitle(strcat(num2str(modOut.meanPowerGen(figNum)), ' [W] Mean Power Generated'))
                            case 4
                                sgtitle(strcat('Data File: ', extractBefore(modOut.dataIn{figNum}, "_")))
                            case 5
                                sgtitle(strcat('Hs: ', num2str(modOut.dataIn.SignificantWaveHeight(figNum)), ' Tw: ', num2str(modOut.dataIn.EnergyPeriod(figNum))))
                        end
                end

                % Plot battery data
                xlabel('Simulation Time [h]');  grid on;
                yyaxis left; ylabel('AUV Battery Level(s) [Wh]'); hold on; xlim([0,170]);
                for auvNum = 1:modOut.fleetSize(figNum)
                    plot(modOut.simTime, modOut.auvBatteryLvl{figNum}(:,auvNum)); 
                end
                map = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]/255;
                set(gca, 'ColorOrder', map, 'NextPlot', 'replacechildren');
                yyaxis right; plot(modOut.simTime, modOut.energyStorageBatteryLvl(:, figNum),'k','LineWidth',1); ylabel('Central Battery Level [Wh]'); 
                ax = gca; % Get current axes
                ax.YColor = 'k'; % Set the color of the right y-axis label and ticks to black % plot(modOut.simTime, modOut.wecBatteryLvl(:, figNum));
                legend([string(num2cell(1:modOut.fleetSize(figNum))), 'Central Battery'], 'Location','northeast');
            end
        
        end  % battery track plot fn
   
        
        %% Plot AUV Time on Mission
        function plotTimeOnMission(modOut)
            %plotTimeOnMission generates plots comparing the time spent 'on
            %mission' between test cases. 
            set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    
            if isempty(modOut.fleetSize)
                modOut.fleetSize = modOut.fleetNumber; 
            end
            
            for i = 1:length(modOut.fleetSize)  % For each system in the simulation batch  
                    t_m = modOut.auvFleet{1,i}{1,1}.missionSpecs(2); 
                    t_r = modOut.auvFleet{1,i}{1,1}.chargeTime;
                    t_c(i) = t_m+t_r;
            end
            
            % Calc aggregate time spent 'on-mission' and determine domain
            if modOut.incorpStagger == 1
                totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMissionCorrected);
                domainTime = (modOut.simTime(end) - max((modOut.fleetSize-1).*(t_c)/modOut.fleetSize, 0) );
            
            else  % no stagger incorporated...
                totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMission);
                domainTime = modOut.simTime(end);

            end


            % Plot time on-mission vs Normalized Power w/ Fleet Size color bar 
            figure; scatter(modOut.ratePwrUsed./modOut.meanPowerGen, 100*totalTimeOnMission./domainTime, [], modOut.fleetSize);
            grid on; xlabel('AUV Power Use Normalized to Mean Power Gen.','Interpreter','latex'); ylabel('Relative Working Hours','Interpreter','latex'); 
            switch modOut.depVar
                case 'AUV Model'
                    switch modOut.resourceDataType
                        case 1
                            title(strcat('Varying AUV Model with Sea State #', modOut.seaState),'Interpreter','latex')
                        case 3
                            title(strcat('Varying AUV Model with constant mean power gen.'),'Interpreter','latex')
                        case {4, 5}
                            title(strcat('Varying AUV Model with constant wave resource specifications'),'Interpreter','latex')
                    end
                case 'WEC Power Gen / Wave Resource'
                    title(strcat('Varying Power Generation', ' with the "', modOut.auvModels, '" AUV'),'Interpreter','latex')
            end
            colorbar; map = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]; colormap(map/255)
            c = colorbar; % Create the colorbar
            c.Label.String = 'Fleet Size'; c.Label.Interpreter = 'lat'; % Assign the colorbar label and use LaTeX as interpreter


            % Plot NORMALIZED time on-mission vs Normalized Power w/ Fleet Size color bar ----------------------
            maxAggregateTime = zeros(size(modOut.fleetSize));  % preallocate
            for i = 1:length(modOut.fleetSize)  % For each system in the simulation batch  
                    t_m = modOut.auvFleet{1,i}{1,1}.missionSpecs(2); 

                    if modOut.incorpStagger == 1
                        t_auvDomain = modOut.simTime(end) - ((t_c(i))/modOut.fleetSize(i))*[0:1:(modOut.fleetSize(i)-1)];  % Time from initial deployment of AUV to the end of the simulation
        
                        % For each auv, how many mission+recharge cycles can they complete from the time they are first deployed to the simulation-end
                        numCycles = t_auvDomain./(t_c(i));
                        max_time_onMission_prelim = floor(numCycles).*t_m + t_m;  % if remainder is > mission time
                        cases = (numCycles - floor(numCycles)) < (t_m./(t_c(i)));
                        max_time_onMission_prelim(cases) = floor(numCycles(cases)).*t_m+ ( (numCycles(cases) - floor(numCycles(cases))) / (t_m/(t_c(i))) )*t_m;  % If remainder is < mission time
        
                        % For each auv, track time spent on-mission before the start of the domain (once all AUVs have been deployed)
                        t_auvTimeSubtract = ((t_c(i))/modOut.fleetSize(i)) * flip([0:1:(modOut.fleetSize(i)-1)]);  % max_time_onMission_prelim included C*(N-1) extra hours for the first auv, C*(N-1)-1 extra hours for the second...etc with C = time between deployments, and N = fleet size
                        t_auvTimeSubtract(t_auvTimeSubtract >= t_m) = t_m; 
        
                        % Final maximum time AUVs could spend on-mission
                        max_time_onMission = max_time_onMission_prelim - t_auvTimeSubtract; 
                        maxAggregateTime(i) = sum(max_time_onMission);  % Sum of max time each AUV in fleet can be on-mission given the initial deployment time and evaluation domain
                    else  % no stagger
                        numCycles = domainTime./(t_c(i));
                        if (numCycles - floor(numCycles)) < (t_m/(t_c(i)))
                            max_time_onMission = floor(numCycles).*t_m+ ( (numCycles - floor(numCycles)) / (t_m/(t_c(i))) )*t_m;
                        else
                            max_time_onMission = floor(numCycles).*t_m + t_m;  % if remainder is > mission time
                        end
                        
                        maxAggregateTime(i) = max_time_onMission*modOut.fleetSize(i);

                    end
            end
                          
            figure; scatter(modOut.ratePwrUsed./modOut.meanPowerGen, 100*totalTimeOnMission./(maxAggregateTime), [], modOut.fleetSize); ylabel('Deployment Efficiency Percentage') % Deployment Efficiency (realized / nominal time on-mission
            % figure; scatter(modOut.ratePwrUsed./modOut.meanPowerGen, 100*totalTimeOnMission./(domainTime.*modOut.fleetSize), [], modOut.fleetSize);  ylabel('Normalized \% On','Interpreter','latex');  % time on-mission / total time
            grid on; xlabel('AUV Power Use Normalized to Mean Power Gen.','Interpreter','latex');
            switch modOut.depVar
                case 'AUV Model'
                    switch modOut.resourceDataType
                        case 1
                            title(strcat('Varying AUV Model with Sea State #', modOut.seaState),'Interpreter','latex')
                        case 3
                            title(strcat('Varying AUV Model with constant mean power gen.'),'Interpreter','latex')
                        case {4, 5}
                            title(strcat('Varying AUV Model with constant wave resource specifications'),'Interpreter','latex')
                    end
                case 'WEC Power Gen / Wave Resource'
                    title(strcat('Varying Power Generation', ' with the "', modOut.auvModels, '" AUV'),'Interpreter','latex')
            end
            colorbar; map = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]; colormap(map/255)
            c = colorbar; % Create the colorbar
            % Set the colorbar label
            c.Label.String = 'Fleet Size'; c.Label.Interpreter = 'lat'; % Assign the label and use LaTeX as interpreter


            % Plot time spent on mission against AUV model labels
            figure; 
            switch modOut.depVar
                case 'AUV Model'
                    categories = categorical(modOut.auvModels); 
                     xlabel('AUV Model','Interpreter','latex');
                case 'WEC Power Gen / Wave Resource'
                    switch modOut.resourceDataType
                        case 1
                            categories = categorical(modOut.seaState);
                            xlabel('Sea State','Interpreter','latex')
                        case 3
                            categories = categorical(modOut.meanPowerGen);
                            xlabel('Mean Power Generated [W]','Interpreter','latex')
                        case 4
                            categories = categorical(extractBefore(modOut.dataIn, "_"));
                            xlabel('Data Files','Interpreter','latex')
                        case 5
                            categories = categorical(1:1:length(modOut.dataIn.EnergyPeriod));
                            xlabel('Wave Resource Case','Interpreter','latex')
                    end
            end
            % scatter(categories, 100*totalTimeOnMission./(maxAggregateTime)); ylabel('Deployment Efficiency Percentage')
            scatter(categories,  100*totalTimeOnMission./domainTime); ylabel('Aggregate Working Hours Percentage')
            grid on;
            switch modOut.depVar
                case 'AUV Model'
                    switch modOut.resourceDataType
                        case 1
                            title(strcat('Varying AUV Model with Sea State #', modOut.seaState))
                        case 3
                            title(strcat('Varying AUV Model with constant mean power gen.'))
                        case {4, 5}
                            title(strcat('Varying AUV Model with constant wave resource specifications'))
                    end
                case 'WEC Power Gen / Wave Resource'
                    title(strcat('Varying Power Generation', ' with the "', modOut.auvModels, '" AUV'))
            end
                    
            
            % % Plot time on-mission vs normalized power, no colorbar
            %{
            figure; scatter(modOut.ratePwrUsed./modOut.meanPowerGen, 100*totalTimeOnMission./domainTime);
            grid on; xlabel('AUV Power Use Rate Normalized to Mean Power Gen.','Interpreter','latex'); ylabel('Percent of Time Spent On Mission','Interpreter','latex'); 
            % Set title
            switch modOut.depVar
                case 'AUV Model'
                    switch modOut.resourceDataType
                        case 1
                            title(strcat('Varying AUV Model with Sea State #', modOut.seaState),'Interpreter','latex')
                        case 3
                            title(strcat('Varying AUV Model with constant mean power gen.'),'Interpreter','latex')
                        case {4, 5}
                            title(strcat('Varying AUV Model with constant wave resource specifications'),'Interpreter','latex')
                    end
                case 'WEC Power Gen / Wave Resource'
                    title(strcat('Varying Power Generation', ' with the "', modOut.auvModels, '" AUV'),'Interpreter','latex')
            end
            %}


        end  % plot time on mission fn


        %% Plot Power Generation
        function plotPowerGen(modOut)
            figure; 
            plot(modOut.simTime, modOut.powerGenMeans)
            grid on; 
            ylabel('Power Generation [W]'); 
            xlabel('Simulation Time [h]'); 
            sgtitle(''); 
            hold on; plot([0,modOut.simTime(end)], mean(modOut.powerGenMeans) * [1,1]);

        end  % plot pwr gen fn

    end  % methods
end  % class