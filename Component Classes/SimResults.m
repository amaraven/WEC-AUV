% Ama Hartman - 2026

classdef SimResults < handle
    % SimResults defines the properties saved in a simResults object.
    % simResults is a lightweight object which stores intermediate
    % simulation values.

    properties 
        fleetSize = 0
        centralBatteryCapacity = 0
        energyStorageBatteryLvl = []
        wecBatteryLvl = []      % nxm matrix of WEC battery levels, with each column corresponding to a wec in the fleet
        auvBatteryLvl = []      % nxm matrix of AUV battery levels, with each column corresponding to battery levels of an auv in the fleet
        auvSchedule = []
        auvTimeOnMission = 0
        auvTimeOnMissionCorrected = 0
        auvFleet = []
        meanPowerGen = 0        % Aggregate mean power generation of the WEC fleet
        powerGenMeans = []      % Aggregate power generation of the WEC fleet for each timestep
        numWECs = 0             % Number of WECs in the fleet
        wecFleet = [];          % Array of WEC objects
        aggWECHotelLoad = 0;    % Aggregate WEC fleet hotel load
        lowPowerGen = 0         % Aggregate low power generation of the WEC fleet
    end
end