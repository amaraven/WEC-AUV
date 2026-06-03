% Ama Hartman

classdef ModelOutput < handle
    % Model Output defines the properties and methods saved in ModelOutput 
    % objects. These objects store model output data, and can generate 
    % plots.
    % 
    % To create a ModelOutput object named modOut, use the following syntax
    % modOut = ModelOutput(simTime)
    
    properties
        simTime (:,1) {mustBeNumeric}  % Time series vector used for simulation. Default is 1 week with 30 s timesteps   
        meanPowerGen {mustBeNumeric}  % [W] Mean value of power gen. Default is generic 1.3 kW
        powerGenMeans {mustBeNumeric}  % [W] Time series of power gen throughout simulation
        dataIn              % Input data specifics
        centralBatteryCapacity (1,:) {mustBeNumeric}  % Vector containing total energy storage of central battery for each simulation [Wh]
        fleetSize (1,:) {mustBeNumeric}  % (1xm) Array with number of AUVs in fleet for each test case
        auvMissionLength (1,:) {mustBeNumeric}  % [hr] Could be 1x1 or 1xm depending on if mission length changes between test cases
        ratePwrUsed (1,:) {mustBeNumeric}  % [W] Rate AUV uses power (power used during mission + recharge [Wh] / mission + recharge time [h])
        energyStorageBatteryLvl {mustBeNumeric} % [Wh] Battery level of central energy storage as a function of time with rows corresponding to time steps and columns corresponding to test cases
        wecBatteryLvl {mustBeNumeric}  % [Wh] Battery level of WEC as a function of time with rows corresponding to time steps and columns corresponding to test cases
        auvBatteryLvl (1,:) cell  % [Wh] Battery level of AUV(s) as a function of time. Cell array with columns corresponding to test case. Each array item is a nxm matrix with n = timesteps and m = auv number
        auvSchedule (1,:) cell  % Operational state schedule of AUV(s) as a function of time with 1) Executing AUV mission, 2) AUV recharging, 3) AUV docked & fully charged. Organized in the same manner as 'auvBatteryLvl'. 
        auvTimeOnMission (1,:) cell  % [h] Time each AUV spends 'on-mission' during the simulation. Cell array with columns corresponding to test case. Each array item is a 1xm matrix with m corresponding to the number of AUVs in fleet.
        auvTimeOnMissionCorrected (1,:) cell  % [h] Time each AUV spends 'on-mission' after all AUVs are deployed in the case of a staggered deployment. Empty if deployment stagger not incorporated. Used for performance calculations using an adjusted domain to exclude time when AUVs are artificially held at dock
        auvFleet (1,:) cell  % Cell array containing auv fleet for each test case
    end
    
    
    %% Instance Methods (need object as an input)
    methods

        %% Constructor
        function modOut = ModelOutput(simTime)
            arguments
                simTime (:,1) {mustBeNumeric} = (30:30:(7*24*60*60))'/60/60
            end

            modOut.simTime = simTime;

        end  % constructor fn


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

        
        %% Plot Central Battery Tracker
        function plotCentralBatteryTracker(modOut, modIn)
            % Generates a plot of the central battery level as a function
            % of time 
            set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
            for figNum = 1:size(modOut.energyStorageBatteryLvl, 2)
            
               % Generate figure, set title
               figure; grid on; hold on;
               switch modIn.depVar
                   case "AUV Model"
                       sgtitle(strcat("AUV:", " ", modIn.auvModels{figNum}));
                   
                   case "WEC Power Gen. / Wave Resource"
                       sgtitle(strcat("Simulation", " ", num2str(figNum)));
               end

               % Plot battery data
               plot(modOut.simTime, modOut.energyStorageBatteryLvl(:, figNum)); 
               xlabel('Simulation Time [h]'); ylabel('Central Battery Level [Wh]'); 
               hold off;
            end

        end


        %% Plot AUV Battery Tracker(s)
        function plotAUVbatteryTracker(modOut, modIn)
            % Generates a plot of the AUV battery level(s) as a function of
            % time
            set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
            
            for figNum = 1:size(modOut.auvBatteryLvl, 2)
                % Generate figure, set title
                figure; hold on; grid on;
                switch modIn.depVar
                    case "AUV Model"
                        sgtitle(strcat("AUV:", " ", modIn.auvModels{figNum}));
    
                    case "WEC Power Gen. / Wave Resource"
                        switch modIn.resourceDataType
                            case 1
                                sgtitle(strcat("Sea State:", " ", num2str(modIn.resourceDataVars.seaState(figNum))));
    
                            otherwise
                                sgtitle(strcat("Simulation", " ", num2str(figNum)));
                        end
                end
    
                % Plot battery data
                for auvNum = 1:modOut.fleetSize(figNum)
                    plot(modOut.simTime, modOut.auvBatteryLvl{figNum}(:, auvNum));
                end
                xlabel("Simulation Time [h]"); ylabel("AUV Battery Level(s) [Wh]");
                customColors = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]/255;
                set(gca, "ColorOrder", customColors, "NextPlot", "replacechildren");
                legend(string(num2cell(1:modOut.fleetSize(figNum))), "Location", "northeast");
            end

        end


        %% Plot Power Gen, Central Battery, AUV Battery (on same panel)
        function plotPowerTrackerPanel(modOut, modIn)
            % Generates a panel including plots of the power generation,
            % central battery level, and AUV battery level(s) throughout a
            % simulation.
            set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

            for figNum = 1:size(modOut.auvBatteryLvl, 2)
                figure('Position', [100, 100, 560, 560]);
                t = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding' ,'compact');
                t.OuterPosition = [0, 0, 1, 1];  % fills the figure

                % Title
                switch modIn.depVar
                    case "AUV Model"
                        sgtitle(strcat("AUV:", " ", modIn.auvModels{figNum}),'Interpreter','latex');

                    case "WEC Power Gen. / Wave Resource"
                        sgtitle(strcat("Simulation", " ", num2str(figNum)),'Interpreter','latex');
                end

                % Power Generation
                ax1 = nexttile;  %subplot(3,1,1)
                grid on; hold on;
                switch modIn.depVar
                    case 'AUV Model'       
                        plot(modOut.simTime, modOut.powerGenMeans)
                        plot([0, modOut.simTime(end)], mean(modOut.powerGenMeans) * [1,1])

                    otherwise
                        plot(modOut.simTime, modOut.powerGenMeans(:, figNum))
                        plot([0, modOut.simTime(end)], mean(modOut.powerGenMeans(:, figNum)) * [1,1])
                end
                % xlabel('Simulation Time [h]');  
                ylabel('Power Generation [W]');
                legend("", "Mean")

                % Central Battery Level
                ax2 = nexttile;  %subplot(3,1,2)
                grid on; hold on;
                plot(modOut.simTime, modOut.energyStorageBatteryLvl(:, figNum));
                % xlabel('Simulation Time [h]'); 
                ylabel('Central Battery Level [Wh]');

                % AUV Battery Level(s)
                ax3 = nexttile;  %subplot(3,1,3)
                hold on; grid on;
                for auvNum = 1:modOut.fleetSize(figNum)
                    plot(modOut.simTime, modOut.auvBatteryLvl{figNum}(:, auvNum))
                end
                xlabel("Simulation Time [h]"); ylabel("AUV Battery Level(s) [Wh]");
                customColors = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]/255;
                set(gca, "ColorOrder", customColors, "NextPlot", "replacechildren");
                legend(string(num2cell(1:modOut.fleetSize(figNum))), "Location", "northeast");

                linkaxes([ax1, ax2, ax3], 'x')
            end

        end


        %% Plot Working Hours
        function plotWorkingHours(modOut, modIn)
            % Generates a plot of the aggregate working hours of the AUV(s)
            % in each simulated system
            set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    
            % Generate figure, set title, set xlabel
            figure; hold on; grid on;
            switch modIn.depVar
                case "AUV Model"
                    xcategories = categorical(string(modIn.auvModels));
                    xlabel("AUV Model");
                    
                    switch modIn.resourceDataType
                        case 1
                            title(strcat("Varying AUV Model with Sea State #", " ", modIn.resourceDataVars.seaState))
                        case {2, 3}
                            title(strcat('Varying AUV Model with Constant Power Gen.'))
                        case {4, 5}
                            title(strcat('Varying AUV Model with Constant Wave Resource Specs'))
                    end

                case "WEC Power Gen. / Wave Resource"
                    title(strcat('Varying Power Generation', ' with the "', string(modIn.auvModels), '" AUV'))

                    switch modIn.resourceDataType
                        case 1
                            xcategories = categorical(string(modIn.resourceDataVars.seaState), string(modIn.resourceDataVars.seaState));
                            xlabel("Sea State");

                        case {2, 4, 6}
                            fileParts = split(string(modIn.resourceDataVars.dataFiles), "\");
                            fileNames = fileParts(:,:, end);
                            xcategories = categorical(strrep(fileNames, "_", "\_"));
                            xlabel("Data Files")

                        case 3
                            xcategories = categorical(modOut.meanPowerGen);
                            xlabel("Mean Power Generated [W]");

                        case 5
                            labels = 1:size(modOut.centralBatteryCapacity, 2);  % Turn into labels with Hs=XX, Te=XX for each simulation...
                            xcategories = categorical(labels);
                    end
            end

            % Calculate aggregate time spent 'on-mission' and domain
            t_c = zeros(1, length(modOut.fleetSize));
            for i = length(modOut.fleetSize):-1:1  % For each system in the simulation batch
                if ~isempty(modOut.auvFleet{1,i})
                    % Can delete later on...
                    if isempty(modOut.auvFleet{1,i}(1,1).chargeTime)
                        rateLinPowerTransfer = modOut.auvFleet{1,i}(1,1).chargeRate;  % Rate of power transfer in linear region
                        modOut.auvFleet{1,i}(1,1).chargeTime = ( 0.8-(1-modOut.auvFleet{1,i}(1,1).missionSpecs(:,3)) )*modOut.auvFleet{1,i}(1,1).batteryCapacity/rateLinPowerTransfer  + 0.4*modOut.auvFleet{1,i}(1,1).batteryCapacity./rateLinPowerTransfer;  % Time to charge up to 80% + time to charge from 80% to 100% 
                    end
    
                    if ~isempty(modOut.auvFleet{1,i})
                        t_m = modOut.auvFleet{1,i}(1,1).missionSpecs(2);  % mission time
                        t_r = modOut.auvFleet{1,i}(1,1).chargeTime;  % recharge time
                        t_c(i) = t_m+t_r;  % cycle time
                    else
                        warning("Issue with domain time calculation during post-processing. Value likely inaccurate.")
                    end
                end
            end

            if modIn.incorpStagger == 1
                totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMissionCorrected);
                domainTime = (modOut.simTime(end) - max((modOut.fleetSize-1).*(t_c)./modOut.fleetSize, 0) );  % Exclude time AUVs are held at dock during the start of the simulation

            else  % no stagger incorporated...
                totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMission);
                domainTime = modOut.simTime(end);

            end

            % Plot global statistics
            scatter(xcategories, 100*totalTimeOnMission./domainTime);
            ylabel("Aggregate Working Hours Percentage");
     
        end


        %% Plot Battery Tracker (AUVs and Central Battery on same plot)
        function plotBatteryTracker(modOut, modIn)
            % plotBatterytracker generates a plot of the central and AUV
            % batteries as a function of time during a simulation. All
            % batteries are overlayed on the same plot
            set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

            % Generate figure, set title
            for figNum = 1:length(modOut.auvBatteryLvl)
                figure; 
                switch modIn.depVar
                    case 'AUV Model'
                        sgtitle(modIn.auvModels{figNum});
                    case 'WEC Power Gen. / Wave Resource'
                        switch modIn.resourceDataType
                            case 1
                                sgtitle(strcat('SeaState_', num2str(modIn.resourceDataVars.seaState(figNum))));
                            otherwise
                               sgtitle(strcat("Simulation", " ", num2str(figNum)));
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
    end  % methods
end  % class