%% Figure: %ON vs normalized Power w/ WEC Diameter colorbar. Max fleet 1 - AUV systems plotted as lines
set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

dataFiles = {'outputData\usWestCoast_maxFleet1_wecDiam0.6_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam0.7_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam0.8_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam0.9_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.0_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.1_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.2_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.3_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.4_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.5_06Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.6_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.7_24Nov25.mat'};
             
normalizedPowerUse = zeros(length(dataFiles), 21);
deploymentEfficiency = zeros(length(dataFiles), 21);
normalizedOnMissionTime = zeros(length(dataFiles), 21);
charDim = zeros(length(dataFiles), 1);
for i = 1:length(dataFiles)
    clearvars -except dataFiles i normalizedPowerUse deploymentEfficiency normalizedOnMissionTime charDim

    % Load model output data
    load(dataFiles{i});

    % Calculate stuff
    totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMissionCorrected);
    domainTime = (modOut.simTime(end) - max(modOut.fleetNumber-1, 0) );

    maxAggregateTime = zeros(size(modOut.fleetNumber));  % preallocate
    for j = 1:length(modOut.fleetNumber)  % For each system in the simulation batch  

        if isempty(modOut.auvFleet{j})
            maxAggregateTime(j) = 1;

        else
            t_m = modOut.auvFleet{1,j}{1,1}.missionSpecs(2); 
            t_r = modOut.auvFleet{1,j}{1,1}.chargeTime;
            t_auvDomain = modOut.simTime(end) - ((t_m+t_r)/modOut.fleetNumber(j))*[0:1:(modOut.fleetNumber(j)-1)];  % Time from initial deployment of AUV to the end of the simulation

            % For each auv, how many mission+recharge cycles can they complete from the time they are first deployed to the simulation-end
            numCycles = t_auvDomain./(t_m + t_r);
            max_time_onMission_prelim = floor(numCycles).*t_m + t_m;  % if remainder is > mission time
            cases = (numCycles - floor(numCycles)) < (t_m./(t_m+t_r));
            max_time_onMission_prelim(cases) = floor(numCycles(cases)).*t_m+ ( (numCycles(cases) - floor(numCycles(cases))) / (t_m/(t_m+t_r)) )*t_m;  % If remainder is < mission time

            % For each auv, track time spent on-mission before the start of the domain (once all AUVs have been deployed)
            t_auvTimeSubtract = ((t_m+t_r)/modOut.fleetNumber(j)) * flip([0:1:(modOut.fleetNumber(j)-1)]);  % max_time_onMission_prelim included C*(N-1) extra hours for the first auv, C*(N-1)-1 extra hours for the second...etc with C = time between deployments, and N = fleet size
            t_auvTimeSubtract(t_auvTimeSubtract >= t_m) = t_m; 

            % Final maximum time AUVs could spend on-mission
            max_time_onMission = max_time_onMission_prelim - t_auvTimeSubtract; 
            maxAggregateTime(j) = sum(max_time_onMission);  % Sum of max time each AUV in fleet can be on-mission given the initial deployment time and evaluation domain

        end
    end

    % save datapoints
    normalizedPowerUse(i,:) = modOut.ratePwrUsed./modOut.meanPowerGen;         
    deploymentEfficiency(i,:) = 100*totalTimeOnMission./(maxAggregateTime);
    normalizedOnMissionTime(i,:) = 100*totalTimeOnMission./(domainTime);
    charDim(i) = wec.charDim;


end

for i = 1:(size(normalizedOnMissionTime,1)-1)
    normalizedOnMissionTime(i, normalizedOnMissionTime(i+1,:) == 0) = NaN;  % Set all preceeding 0s to NaN except for the last one 
    deploymentEfficiency(i, deploymentEfficiency(i+1, :) == 0) = NaN;

end

% Keep only ONE 'max deployment efficiency' datapoint (the one correlating to the lowest WEC size, as long as dataFiles is in order least -> most)
[~, maxIndx] = max(deploymentEfficiency, [], 1);
for jj = 1:length(maxIndx)
    deploymentEfficiency(maxIndx(jj)+1:end, jj) = NaN;
end


% Figure - select AUVs
figure; hold on;
for k = [3,4,9,13,21]%[3,4,9,12,13]%1:21 % [3,4,6,9,12,13]
% scatter(modOut.ratePwrUsed.*modOut.fleetNumber./modOut.meanPowerGen, 100*totalTimeOnMission./(maxAggregateTime), [], wec.charDim*ones(size(modOut.ratePwrUsed)));
    % plot(normalizedPowerUse(:,k), deploymentEfficiency(:,k), 'k-','LineWidth', 0.5); ylabel('Deployment Efficiency Percentage')
    % scatter(normalizedPowerUse(:,k), deploymentEfficiency(:,k), [], charDim, 'LineWidth',1.5); ylabel('Deployment Efficiency Percentage')
    plot(normalizedPowerUse(:,k), normalizedOnMissionTime(:,k), 'k-','LineWidth', 0.5); 
    scatter(normalizedPowerUse(:,k), normalizedOnMissionTime(:,k), [], charDim, 'LineWidth',1.5); ylabel('Working Hour Percentage','Interpreter','latex'); 
    % scatter(normalizedPowerUse(:,k), normalizedOnMissionTime(:,k), [], wec.charDim*ones(size(modOut.ratePwrUsed)));
    grid on; xlabel('AUV Power Use Normalized to Mean Power Gen.','Interpreter','latex'); 
    title('Varying AUV Model');
    colorbar; map = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]; colormap(map/255)
    c = colorbar; % Create the colorbar
    % Set the colorbar label
    c.Label.String = 'WEC Characteristic Dimension'; c.Label.Interpreter = 'lat'; % Assign the label and use LaTeX as interpreter
end
% Figure - All AUVs
figure; hold on;
for l = 1:21
    plot(normalizedPowerUse(deploymentEfficiency(:,l)~=0,l), deploymentEfficiency(deploymentEfficiency(:,l)~=0,l),'color',[0.25,0.25,0.25],'linewidth',0.5); 
    scatter(normalizedPowerUse(deploymentEfficiency(:,l)~=0,l), deploymentEfficiency(deploymentEfficiency(:,l)~=0,l), [], charDim(deploymentEfficiency(:,l)~=0), 'o','filled','LineWidth',0.25); ylabel('Deployment Efficiency Percentage')
    ylim([0 100])
    grid on; xlabel('Normalized AUV Power Use','Interpreter','latex'); 
    colorbar; map = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]; colormap(map/255)
    c = colorbar; % Create the colorbar
    % Set the colorbar label
    c.Label.String = 'WEC Characteristic Dimension [m]'; c.Label.Interpreter = 'lat'; % Assign the label and use LaTeX as interpreter
end

%% Figure: %ON vs normalized Power w/ WEC Diameter colorbar. Max fleet 1
%{
set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

dataFiles = {'outputData\usWestCoast_maxFleet1_wecDiam0.6_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam0.7_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam0.8_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam0.9_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.0_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.1_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.2_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.3_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.4_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.5_06Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.6_24Nov25.mat'
             'outputData\usWestCoast_maxFleet1_wecDiam1.7_24Nov25.mat'};
             
figure; hold on;

for i = 1:length(dataFiles)
    clearvars -except dataFiles i

    % Load model output data
    load(dataFiles{i});

    % Calculate stuff
    totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMissionCorrected);
    domainTime = (modOut.simTime(end) - max(modOut.fleetNumber-1, 0) );

    maxAggregateTime = zeros(size(modOut.fleetNumber));  % preallocate
    for j = 1:length(modOut.fleetNumber)  % For each system in the simulation batch  

        if isempty(modOut.auvFleet{j})
            maxAggregateTime(j) = 1;

        else
            t_m = modOut.auvFleet{1,j}{1,1}.missionSpecs(2); 
            t_r = modOut.auvFleet{1,j}{1,1}.chargeTime;
            t_auvDomain = modOut.simTime(end) - ((t_m+t_r)/modOut.fleetNumber(j))*[0:1:(modOut.fleetNumber(j)-1)];  % Time from initial deployment of AUV to the end of the simulation

            % For each auv, how many mission+recharge cycles can they complete from the time they are first deployed to the simulation-end
            numCycles = t_auvDomain./(t_m + t_r);
            max_time_onMission_prelim = floor(numCycles).*t_m + t_m;  % if remainder is > mission time
            cases = (numCycles - floor(numCycles)) < (t_m./(t_m+t_r));
            max_time_onMission_prelim(cases) = floor(numCycles(cases)).*t_m+ ( (numCycles(cases) - floor(numCycles(cases))) / (t_m/(t_m+t_r)) )*t_m;  % If remainder is < mission time

            % For each auv, track time spent on-mission before the start of the domain (once all AUVs have been deployed)
            t_auvTimeSubtract = ((t_m+t_r)/modOut.fleetNumber(j)) * flip([0:1:(modOut.fleetNumber(j)-1)]);  % max_time_onMission_prelim included C*(N-1) extra hours for the first auv, C*(N-1)-1 extra hours for the second...etc with C = time between deployments, and N = fleet size
            t_auvTimeSubtract(t_auvTimeSubtract >= t_m) = t_m; 

            % Final maximum time AUVs could spend on-mission
            max_time_onMission = max_time_onMission_prelim - t_auvTimeSubtract; 
            maxAggregateTime(j) = sum(max_time_onMission);  % Sum of max time each AUV in fleet can be on-mission given the initial deployment time and evaluation domain

        end
    end
                  
    scatter(modOut.ratePwrUsed./modOut.meanPowerGen, 100*totalTimeOnMission./(maxAggregateTime), [], wec.charDim*ones(size(modOut.ratePwrUsed))); ylabel('Deployment Efficiency Percentage')
    % scatter(modOut.ratePwrUsed.*modOut.fleetNumber./modOut.meanPowerGen, 100*totalTimeOnMission./(domainTime), [], wec.charDim*ones(size(modOut.ratePwrUsed))); ylabel('Normalized \% On'); 
    grid on; xlabel('AUV Power Use Normalized to Mean Power Gen.','Interpreter','latex'); ylabel('Normalized \% On','Interpreter','latex'); 
    title('Varying AUV Model');
    colorbar; map = [61 32 44; 120 63 87; 35 87 137; 80 145 145; 91 140 90; 162 154 58; 252 171 16; 219 110 48; 202 80 64; 185 49 79]; colormap(map/255)
    c = colorbar; % Create the colorbar
    % Set the colorbar label
    c.Label.String = 'WEC Characteristic Dimension'; c.Label.Interpreter = 'lat'; % Assign the label and use LaTeX as interpreter

end
%}


%% Figure: %ON vs Fleet Number
set(groot, 'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% Plot for Iver3-27, AUV #1
for auvNum = 3%1:21;   % 18: Boxfish AUV, max fleet 14  % 13: bluefin HAUV, N=
dataFiles = {'outputData\usWestCoast_maxFleet1_05Nov25.mat'
             'outputData\usWestCoast_maxFleet2_05Nov25.mat'
             'outputData\usWestCoast_maxFleet3_05Nov25.mat'
             'outputData\usWestCoast_maxFleet4_05Nov25.mat'
             'outputData\usWestCoast_maxFleet5_05Nov25.mat'
             'outputData\usWestCoast_maxFleet6_05Nov25.mat'
             'outputData\usWestCoast_maxFleet7_05Nov25.mat'
             'outputData\usWestCoast_maxFleet8_05Nov25.mat'
             'outputData\usWestCoast_maxFleet9_05Nov25.mat'
             'outputData\usWestCoast_maxFleet10_05Nov25.mat'
             'outputData\usWestCoast_maxFleet11_05Nov25.mat'};
             % 'outputData\usWestCoast_maxFleet12_31Oct25.mat'
             % 'outputData\usWestCoast_maxFleet13_31Oct25.mat'
             % 'outputData\usWestCoast_maxFleet14_31Oct25.mat'
             % 'outputData\usWestCoast_maxFleet15_31Oct25.mat'};

for i = 1:length(dataFiles)
    clearvars -except dataFiles auvNum i

    % Load model output data
    load(dataFiles{i});

    % Calculate stuff
    if modOut.incorpStagger == 1
        totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMissionCorrected);
        for k = 1:length(modOut.fleetNumber)
            staggerHours(k) = (modOut.auvFleet{1,k}{1,1}.missionSpecs(2) + modOut.auvFleet{1,k}{1,1}.chargeTime) / modOut.fleetNumber(k);
        end

        domainTime = modOut.simTime(end)-staggerHours.*(modOut.fleetNumber-1);


    else  % no stagger incorporated...
        totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMission);
        domainTime = modOut.simTime;
    end

    maxAggregateTime = zeros(size(modOut.fleetNumber));
    for j = 1:length(modOut.fleetNumber)  % For each system in the simulation batch  
        t_m = modOut.auvFleet{1,j}{1,1}.missionSpecs(2); 
        t_r = modOut.auvFleet{1,j}{1,1}.chargeTime;
        t_auvDomain = modOut.simTime(end) - ((t_m+t_r)/modOut.fleetNumber(j))*[0:1:(modOut.fleetNumber(j)-1)];  % Time from initial deployment of AUV to the end of the simulation

        % For each auv, how many mission+recharge cycles can they complete from the time they are first deployed to the simulation-end
        numCycles = t_auvDomain./(t_m + t_r);
        max_time_onMission_prelim = floor(numCycles).*t_m + t_m;  % if remainder is > mission time
        cases = (numCycles - floor(numCycles)) < (t_m./(t_m+t_r));
        max_time_onMission_prelim(cases) = floor(numCycles(cases)).*t_m+ ( (numCycles(cases) - floor(numCycles(cases))) / (t_m/(t_m+t_r)) )*t_m;  % If remainder is < mission time

        % For each auv, track time spent on-mission before the start of the domain (once all AUVs have been deployed)
        t_auvTimeSubtract = ((t_m+t_r)/modOut.fleetNumber(j)) * flip([0:1:(modOut.fleetNumber(j)-1)]);  % max_time_onMission_prelim included C*(N-1) extra hours for the first auv, C*(N-1)-1 extra hours for the second...etc with C = time between deployments, and N = fleet size
        t_auvTimeSubtract(t_auvTimeSubtract >= t_m) = t_m; 

        % Final maximum time AUVs could spend on-mission
        max_time_onMission = max_time_onMission_prelim - t_auvTimeSubtract; 
        maxAggregateTime(j) = sum(max_time_onMission);
    end

    if i == 1
        figure(auvNum+100); 
        yyaxis left
        xlabel('Fleet Size');
        hold on; grid on;

        yyaxis right
        ylabel('Working Hours Percentage')% Normalized to Fleet Size')
        hold on; grid on;
    end

    figure(auvNum+100)

    % Plot performance given fleet number
    yyaxis left
    % scatter(modOut.fleetNumber(auvNum),  100*totalTimeOnMission(auvNum)./maxAggregateTime(auvNum));  ylabel('Deployment Efficiency Percentage'); % 'deployment efficiency'
    scatter(modOut.fleetNumber(auvNum),  100*totalTimeOnMission(auvNum)./(domainTime(auvNum).*modOut.fleetNumber(auvNum)));  ylabel('Normalized Working Hours Percentage') 

    % Plot performance given fleet number
    yyaxis right 
    scatter(modOut.fleetNumber(auvNum),  100*totalTimeOnMission(auvNum)./domainTime(auvNum));%./modOut.fleetNumber(auvNum)); 

end
end



%%
% totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMission);
% 
% if exist modOut.incorpStagger
%     if modOut.incorpStagger == 1
%         totalTimeOnMission = totalTimeOnMission - cellfun(@sum, modOut.auvTimeOutsideDomain);
%     end
% end
% 
% for i = 1:length(modOut.fleetNumber)
%     t_m(i) = modOut.auvFleet{1,i}{1,1}.missionSpecs(2);  % t_m = mission time
%     t_r(i) = modOut.auvFleet{1,i}{1,1}.chargeTime;  % t_r = recharge time
% end
% num_cycles = modOut.simTime(end)./(t_m + t_r);  % # mission+recharge cycles fit in a simulation
% 
% max_time_onMission = zeros(size(modOut.fleetNumber));  % preallocate
% for j = 1:length(modOut.fleetNumber)
%     if (num_cycles(j) - floor(num_cycles(j))) >= (t_m(j)/(t_r(j)+t_m(j)))
%         max_time_onMission(j) = floor(num_cycles(j))*t_m(j) + t_m(j);
% 
%     elseif (num_cycles(j) - floor(num_cycles(j))) < (t_m(j)/(t_r(j)+t_m(j)))
%         max_time_onMission(j) = floor(num_cycles(j))*t_m(j) + t_m(j)* ((num_cycles(j) - floor(num_cycles(j))) / (t_m(j)/(t_r(j)+t_m(j))));
%     else
%         warning('oh no')
%     end
% end
% 
% max_time_onMission2 = floor(num_cycles).*t_m + t_m; 
% cases = (num_cycles - floor(num_cycles)) < (t_m./(t_m+t_r));
% max_time_onMission2(cases) = floor(num_cycles(cases)).*t_m(cases) + ( (num_cycles(cases) - floor(num_cycles(cases))) ./ (t_m(cases)./(t_m(cases)+t_r(cases))) ).*t_m(cases);



% for i = 1:length(modOut.fleetNumber)
%     staggerHours(i) = modOut.auvFleet{1,i}{1,1}.missionSpecs(2)+modOut.auvFleet{1,i}{1,1}.chargeTime/modOut.fleetNumber(i);
%     int = 1:1:modOut.fleetNumber(i);
%     fnSum(i) = sum(int(1:end-1))*staggerHours(i); 
% 
% end

% %% Plots from model outputs
% 
% dataFolder = 'outputData';
% dataFiles = {'Hawaii_v_auvModel_C1.25_15Aug25.mat';  %'Hawaii_v_auvModel_C1_15Aug25.mat'; too small to support most auvs
%              'Hawaii_v_auvModel_C1.5_15Aug25.mat';
%              'Hawaii_v_auvModel_C1.75_15Aug25.mat';
%              'Hawaii_v_auvModel_C2_15Aug25.mat';
%              'Hawaii_v_auvModel_C2.25_15Aug25.mat';
%              'Hawaii_v_auvModel_C2.5_15Aug25.mat'
%              'Alaska_v_auvModel_C1.5_15Aug25.mat';
%              'usWestCoast_v_auvModel_C1.5_15Aug25.mat'
%              'Alaska_v_auvModel_C2_15Aug25.mat';
%              'usWestCoast_v_auvModel_C2_15Aug25.mat'};
% 
% figure; hold on; grid on;
% xlabel('AUV Power Use Rate Normalized to Mean Power Gen.'); ylabel('Percent of Time Spent On Mission');
% title('Varying AUV Model with Hawaii Hindcast Wave Resource')
% for i = 1:length(dataFiles)
%     load(fullfile(dataFolder, dataFiles{i})); 
% 
%    % Polynomial fitting is the wrong tool. looks more like y = a/x; 
%     % missionTimeFit(:,i) = polyfit(modOut.ratePwrUsed./modOut.meanPowerGen, 100*totalTimeOnMission/modOut.simTime(end), 3);
%     % missionTimeFitX = linspace(0,2);
%     % missionTimeFitY = polyval(missionTimeFit(:,i), missionTimeFitX);
% 
%     % Plot % time on mission against normalized AUV power use
%     totalTimeOnMission = cellfun(@sum, modOut.auvTimeOnMission);
% 
%     scatter(modOut.ratePwrUsed./(modOut.meanPowerGen), 100*totalTimeOnMission/modOut.simTime(end));
%     % scatter(modOut.ratePwrUsed./(modOut.meanPowerGen-wec.hotelLoad), 100*totalTimeOnMission/modOut.simTime(end));
%     % scatter(modOut.ratePwrUsed./modOut.meanPowerGen, modOut.fleetNumber);
% 
% 
% end
% 

%% FLEET # CALCULATION
% Using 'wave fleet number' and a high-level estimate 'possible fleet
% number', then manually checking if they work. This scratch is to tease
% out some other variables that may be at play. 
%{

%     waveFleetNumber = floor((powerGen40p * wec.n_battery - wec.hotelLoad/wec.n_battery) / (auv.hotelLoad / wec.n_battery / auv.n_battery / auv.n_powerTransfer));  % PowerGen - wecHotelLoad - auvHotelLoad*fleetNumber > 0
%     cycleTime = auv.chargeTime(auv.mission) + auv.missionSpecs(auv.mission,2);
%     possibleFleet = floor((meanPowerGen*wec.n_battery - wec.hotelLoad/wec.n_battery)*cycleTime / (auv.chargeLoad/wec.n_battery) + 0.5);


wec = WEC('generic');
auvModel =[{'Iver3-27'}, {'Iver3-38.5'}, {'REMUS 100'}, {'REMUS 300-58.5'}, {'REMUS 300-70.3'}, {'REMUS 620-210'}, {'REMUS 620-279'}, {'REMUS 620-347'}, {'REMUS 6000'}, {'Bluefin-9'}, {'Bluefin-12'}, {'Bluefin-21'}, {'Bluefin-HAUV'}, {'Hugin Superior'}, {'Hugin Endurance'}, {'Hugin 3000'}, {'Hugin 4500'}, {'Boxfish AUV'}, {'Boxfish ARV-i'}];

% Calculate power generation (vector)
seaState = 3; 
wec.calcPowerGen(RM3, seaState, simHrs, dt);  % Output is wec.powerGenMeans: Timeseries of power generation [W] corresponding to simulation time.
meanPowerGen = mean(wec.powerGenMeans);
powerGen40p = prctile(wec.powerGenMeans, 40);  % Using 40th percentile for fleetNumber calculation because 'mean' often yields a slight overestimation for the number of AUV's. 40th %ile adds a safety factor.
lowPowerGen = prctile(wec.powerGenMeans, 25);  % p ends up being less than 0.1%...

chargeRateOvrMeanPwr(:,1) = 1:length(auvModel);
chargeLoadOvrMeanPwr(:,1) = 1:length(auvModel);
pwrAvailableFleet(:,1) = 1:length(auvModel);
maxHotelFleet(:,1) = 1:length(auvModel);

for iter = 1:length(auvModel)
    auv = AUV(auvModel{iter}); 

    t_wecRecharge = max((auv.chargeLoad(auv.mission)/wec.n_battery - (meanPowerGen*wec.n_battery - wec.hotelLoad/wec.n_battery)*auv.chargeTime(auv.mission) ) / (wec.n_battery*meanPowerGen - wec.hotelLoad/wec.n_battery), 0);  % If negative, generated more power than used during auv recharge, power dumped instead, and wec recovery time is 0
    
    % Fleet num based on number of times WEC can recharge (after a full AUV charge), within a single AUV mission  
    % fleetNumber = max(ceil(auv.missionSpecs(auv.mission,2) / (auv.chargeLoad/(wec.n_battery*meanPowerGen)))+1, 1); 

    % fleet num based on total power available
    cycleTime = auv.chargeTime(auv.mission) + auv.missionSpecs(auv.mission,2);  % Time for an auv to go on mission + recharge
    pwrAvailableFleet(iter,2) = (meanPowerGen*wec.n_battery - wec.hotelLoad/wec.n_battery)*cycleTime / (auv.chargeLoad/wec.n_battery) + 0.5;  % number of auv's that can recharge given the power available + a little, to account for docked auvs

    % Calculate MAX number of AUV's WEC could support in dock with given/calculated power generation
    maxHotelFleet(iter, 2) = (meanPowerGen * wec.n_battery - wec.hotelLoad/wec.n_battery) / (auv.hotelLoad / wec.n_battery / auv.n_battery / auv.n_powerTransfer);  % PowerGen - wecHotelLoad - auvHotelLoad*fleetNumber > 0    

    % Ratio between charge rate & mean power gen
    chargeRateOvrMeanPwr(iter,2) = auv.chargeRate / meanPowerGen;
    chargeLoadOvrMeanPwr(iter,2) = auv.chargeLoad / meanPowerGen;
    % numRecharge(iter,1) = max(ceil(auv.missionSpecs(auv.mission,2) /
    % (auv.chargeLoad/(wec.n_battery*meanPowerGen)))+1, 1);  % flawed.
    % chargeLoad is just the draw and doesn't take into account any power
    % generated. 

    numWECrecharge(iter, 1) = ceil(auv.missionSpecs(auv.mission, 2) / t_wecRecharge);
end

chargeRateOvrMeanPwr(:,3) = chargeRateOvrMeanPwr(:,2) > 0.6;
chargeRateOvrMeanPwr(:,4) = 2*(chargeRateOvrMeanPwr(:,2) > 3.0);  % fails on x3 instances
chargeLoadOvrMeanPwr(:,3) = chargeLoadOvrMeanPwr(:,2) > 5; 
chargeLoadOvrMeanPwr(:,4) = 2*(chargeLoadOvrMeanPwr(:,2) > 30);%.*chargeLoadOvrMeanPwr(:,3)  % Fails on x4 instances
%}

% if t_wecRecharge > auv.missionSpecs(auv.mission, 2)  % If an auv mission is shorter than the time it takes the WEC to recharge after an AUV-recharge...
%         fleetNumber(depVar) = 1;
%     else
% 
%         numWECrecharge = ceil(auv.missionSpecs(auv.mission,2) / t_wecRecharge);
%         cycleTime = auv.chargeTime(auv.mission) + auv.missionSpecs(auv.mission,2);
%         possibleFleet = floor((meanPowerGen*wec.n_battery - wec.hotelLoad/wec.n_battery)*cycleTime / (auv.chargeLoad/wec.n_battery) + 0.5);
% 
%        if possibleFleet <= maxHotelFleet  % If the number of times the wec can recharge within a single auv mission is less than the fleet number calculated from the wave resource..
%             % t_wecRechargeUpdated = auv.chargeLoad(auv.mission)/((wec.n_battery*meanPowerGen - wec.hotelLoad/wec.n_battery - auv.hotelLoad*(numWECrecharge/2)/wec.n_battery/auv.n_powerTransfer));  % hotelLoad*(fleetNumber(depVar)-1)
%             % fleetNumber(depVar) = ceil(auv.missionSpecs(auv.mission,2) / t_wecRechargeUpdated);
%             % disp('Lim. by wec recharge.')
%             fleetNumber(depVar) = possibleFleet;
%             disp('Lim. by possibleFleet')
%        else
%            fleetNumber(depVar) = maxHotelFleet;
%         % fleetNumber(depVar) = min( fleetNumber(depVar), ceil(auv.missionSpecs(auv.mission,2) / t_wecRecharge) ); %(auv.chargeLoad/(wec.n_battery*meanPowerGen))) );
%        end
% 
%     end