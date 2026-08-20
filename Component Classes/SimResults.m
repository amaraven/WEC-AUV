% SimResults defines the properties saved in a simResults object.
% simResults is a lightweight object which stores intermediate
% simulation values for the power modeling tool.
%
% Ama Hartman - 2026


classdef SimResults < handle

    properties 
        fleetSize = 0           % Size of AUV fleet
        centralBatteryCapacity = 0      % [Wh] Central battery capacity
        energyStorageBatteryLvl = []    % (nx1) vector of the energy storage battery level during the simulation
        wecBatteryLvl = []      % (nxm) matrix of WEC battery levels, with each column corresponding to a wec in the fleet
        auvBatteryLvl = []      % (nxm) matrix of AUV battery levels, with each column corresponding to battery levels of an auv in the fleet
        auvSchedule = []        % (nxm) matrix of AUV schedules, with each column corresponding to the timestep operational states of each AUV
        auvTimeOnMission = 0    % (1xm) [hr] Value of the amount of time each AUV spends on-mission during the simulation
        auvTimeOnMissionCorrected = 0   % (1xm) [hr] The amount of time each AUV spends on mission during the simulation minus any time on-mission before all AUVs are deployed
        auvFleet = []           % [1xm] array of AUV objects
        meanPowerGen = 0        % [W] Aggregate mean power generation of the WEC fleet
        powerGenMeans = []      % [W] Aggregate power generation of the WEC fleet for each timestep
        numWECs = 0             % Size of WEC fleet
        wecFleet = [];          % (1xp) Array of WEC objects
        aggWECHotelLoad = 0;    % [W] Aggregate WEC fleet hotel load
        lowPowerGen = 0         % [W] Aggregate low power generation of the WEC fleet
    end
end