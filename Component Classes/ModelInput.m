% Ama Hartman

classdef ModelInput < handle
    % ModelInput defines the properties and methods saved in ModelInput
    % objects. These objects store all inputs neededfor the power model simulation. 
    %
    % User-Defined properties: 
    % - simHrs: [hr] simulation duration. Default 1 wk
    % - depVar: variable to change between simulation iterations. Default
    %   'AUV Model'
    % - resourceDataType: Type of resource data. Power in W, time in h,
    %   distance in m, wave periods in s.
    %   1. Modeled power gen. of 1.5 m WEC in different sea states (1-10)
    %   2. Power generation time series and time vector (*.mat file) 
    %   3. Value of mean power generation
    %   4. Time series of wave specifications (significant wave
    %   height, wave energy and peak periods, and time) (Default)
    %   5. Values of mean wave specifications 
    %   6. Power matrix and wave spec. time series
    % - resourceDataVars: struct containing metadata about user-provided
    %   resource data
    % - incorpStagger: Enforce a stagger between AUV deployments. Default 1
    %   (on)
    % - maxFleetSize: Enforce a maximum AUV fleet size. Default 0 (off)
    % - userDefinedBattery: Size of user-defined central battery capacity. 
    %   Default 0 (off)
    % - dtSec: timestep interval in seconds. Default 30 s
    % - auvModels: cell array containing AUV model name(s)
    % - outputPlots: Struct containing output plots selected to generate

    % User-defined properties
    properties (GetAccess = public, SetAccess = private)
        simHrs (1,1) {mustBePositive} = (7*24)  % default: 1 wk
        depVar (1,1) string {mustBeMember(depVar, ["AUV Model", "WEC Power Gen. / Wave Resource"])} = "AUV Model"  % default: compare AUV models
        resourceDataType (1,1) {mustBeMember(resourceDataType, [1,2,3,4,5,6])} = 4  % default: 4 (time series of wave specs)
        resourceDataVars struct = struct()  % Struct containing resource data information and variables
        incorpStagger (1,1) {mustBeMember(incorpStagger, [1,0])} = 1  % default: 1 (yes incorporate stagger into AUV mission scheduling)
        maxFleetSize {mustBeNonnegative} = 0  % default: 0 (no max fleet size)
        userDefinedBattery {mustBeNonnegative} = 0  % default: 0 (no user-defined battery)
        dtSec (1,1) {mustBeNonnegative} = 30  % default: 30-second timestep intervals
        auvModels (1,:) cell = [{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}, {'G'}, {'H'}, {'I'}, {'J'}, {'K'}, {'L'}, {'M'}, {'N'}, {'O'}, {'P'}, {'Q'}, {'R'}, {'S'}, {'T'}, {'U'}] % Cell containing all auvModels to use in simulations. Models are defined in 'AUV' class script.
        outputPlots (1,1) struct = struct()  % Struct containing output plots selected
    end

    % Internal properties
    properties (Dependent = true)  % properties dependent on other props
        dt % [hr]
        simTime % [hr]
    end


    %% Instance Methods (need object as an input)
    methods
        %% Constructor
        function mi = ModelInput(simHrs, depVar, resourceDataType, incorpStagger, maxFleetSize, userDefinedBattery, dtSec, auvModels, resourceDataVariables, outputPlotSelection)
            arguments
                simHrs (1,1) {mustBePositive} = (7*24)  % default: 1 wk
                depVar (1,1) string {mustBeMember(depVar, ["AUV Model", "WEC Power Gen. / Wave Resource"])} = "AUV Model"  % default: compare AUV models
                resourceDataType (1,1) {mustBeMember(resourceDataType, [1,2,3,4,5,6])} = 4  % default: 4 (time series of wave specs)
                incorpStagger (1,1) {mustBeMember(incorpStagger, [1,0])} = 1  % default: 1 (yes incorporate stagger)
                maxFleetSize {mustBeNonnegative} = 0  % default: 0 (no max fleet size)
                userDefinedBattery {mustBeNonnegative} = 0  % default: 0 (no user-defined battery)
                dtSec (1,1) {mustBeNonnegative} = 30  % default: 30-second timestep intervals
                auvModels (1,:) cell = [{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}, {'G'}, {'H'}, {'I'}, {'J'}, {'K'}, {'L'}, {'M'}, {'N'}, {'O'}, {'P'}, {'Q'}, {'R'}, {'S'}, {'T'}, {'U'}]
                resourceDataVariables struct = []  % Struct containing all resource data variables
                outputPlotSelection struct = struct();
            end

            if isempty(maxFleetSize)
                maxFleetSize = 0;
            end

            if isempty(userDefinedBattery)
                userDefinedBattery = 0;
            end

            mi.simHrs = simHrs;
            mi.depVar = depVar;
            mi.resourceDataType = resourceDataType; 
            mi.incorpStagger = incorpStagger;
            mi.maxFleetSize = maxFleetSize;
            mi.userDefinedBattery = userDefinedBattery;
            mi.dtSec = dtSec; 
            mi.auvModels = auvModels;
            mi.resourceDataVars = resourceDataVariables;
            mi.outputPlots = outputPlotSelection;
        end


        %% Dependent property calcs
        function dt = get.dt(mi)
            dt = mi.dtSec /60 /60;  % [hr]
        end

        function simTime = get.simTime(mi)
            simSec = mi.simHrs *60 *60;  % [s]
            simTime = (mi.dtSec:mi.dtSec:simSec)' /60 /60;  % [hr]
        end


        %% Load Resource Data, convert to powerGen time series
        function calcPowerGen(mi, wec, modOut, depVarCount)
            % From user-input data, calculates power generation vector and
            % saves in modOut object.
            % Inputs: 
            % * wec: WEC class object
            % * modOut: ModelOutput class object
            % * depVarCount: (Optional) iteration of the dependent variable
            %   loop, exclude if calling from outside of the loop
            arguments
                mi ModelInput
                wec WEC
                modOut % ModelOutput
                depVarCount (1,1) {mustBePositive, mustBeInteger} = 1
            end

            % Calculate and save power gen
            switch mi.resourceDataType
                case 1 % Proteus struct containing power generated during different sea states
                    % load data
                    load('inputData/RM3_seaState_Power.mat', 'RM3');  % Time series output of proteus WEC model for a single point absorber with no dock attached.
                    
                    % extract seaState to use for power calc.
                    switch mi.depVar
                        case 'AUV Model'
                            seaState = str2double(mi.resourceDataVars.seaState);
                        case 'WEC Power Gen. / Wave Resource'
                            seaState = depVarCount;
                    end
                    
                    % calc
                    wec.reshapePowerGen(RM3(seaState).Power, RM3(seaState).Time/60/60, mi.simHrs, mi.dt);  % Time must be in hours
                    
                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;  % save power gen
                    modOut.dataIn.seaState(depVarCount) = seaState;  % save sea state

                case 2 % Time series of power gen
                    % load data
                    dataTime = load(mi.resourceDataVars.dataFiles{depVarCount}, mi.resourceDataVars.tVarName);
                    pGen_dataTime = load(mi.resourceDataVars.dataFiles{depVarCount}, mi.resourceDataVars.pwrVarName);

                    % reshape power gen. to fit simulation timestep
                    wec.reshapePowerGen(pGen_dataTime.(mi.resourceDataVars.pwrVarName), dataTime.(mi.resourceDataVars.tVarName), mi.simHrs, mi.dt);  % Time must  be in hr
                     
                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;

                case 3  % Value of mean power
                    % save
                    if isempty(modOut.meanPowerGen)
                        modOut.meanPowerGen = mi.resourceDataVars.meanPower;
                    end
    
                    % save for use in current sim. iteration
                    wec.meanPowerGen = modOut.meanPowerGen(depVarCount);
                    wec.lowPowerGen = 0.75 * wec.meanPowerGen; 
                    wec.powerGenMeans = ones(size(mi.simTime)) * wec.meanPowerGen; 

                case 4 % Time series of wave specs (Hs, Te, Tp)
                    % Determine filetype, load data, & build data struct
                    [~, ~, ext] = fileparts(mi.resourceDataVars.dataFiles{depVarCount});
                    if strcmp(ext, '.csv')
                        mi.resourceDataVars.dataTable = readtable(mi.resourceDataVars.dataFiles{depVarCount},'VariableNamingRule','modify');
                        
                        waveData.sigWaveHeight = mi.resourceDataVars.dataTable.SignificantWaveHeight; 
                        waveData.waveEnergyPeriod = mi.resourceDataVars.dataTable.EnergyPeriod;
                        waveData.peakPeriod = mi.resourceDataVars.dataTable.PeakPeriod;
                        dataDateTimes = datetime(mi.resourceDataVars.dataTable.Year, mi.resourceDataVars.dataTable.Month, mi.resourceDataVars.dataTable.Day, mi.resourceDataVars.dataTable.Hour, mi.resourceDataVars.dataTable.Minute, zeros(size(mi.resourceDataVars.dataTable.Year)));
                        waveData.dataTime = hours(dataDateTimes-dataDateTimes(1) + (dataDateTimes(2)-dataDateTimes(1)));

                    elseif strcmp(ext, '.mat')
                        rawData = load(mi.resourceDataVars.dataFiles{depVarCount});

                        waveData.sigWaveHeight = rawData.(mi.resourceDataVars.HsVarName);
                        waveData.waveEnergyPeriod = rawData.(mi.resourceDataVars.TeVarName);
                        waveData.peakPeriod = rawData.(mi.resourceDataVars.TpVarName); 
                        waveData.dataTime = rawData.(mi.resourceDataVars.tVarName);

                    else
                        error('Data filetype is not currently supported. Supported filetypes include *.csv and *.mat.')
                    end

                    % calculate
                    % wec.calcPowerGen(mi.resourceDataVars.dataTable,'meanPwr', mi.simTime, 0, 403); %% set windowOverrideIndx to 403 (and use Oregon dataset) to replicate paper results, otherwise set to 0  
                    wec.calcPowerGen(waveData,'meanPwr', mi.simTime, 0, 403); %% set windowOverrideIndx to 403 (and use Oregon dataset) to replicate paper results, otherwise set to 0  

                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;

                case 5 % Mean Wave specs (Hs, Te, Tp)
                    waveData.sigWaveHeight = mi.resourceDataVars.sigWaveHeight(depVarCount); 
                    waveData.waveEnergyPeriod = mi.resourceDataVars.energyPeriod(depVarCount);
                    waveData.peakPeriod = mi.resourceDataVars.peakPeriod(depVarCount);
                    % calculate
                    % wec.calcPowerGen(mi.resourceDataVars.waveSpecTable(depVarCount,:), [], mi.simTime, 0, 0);
                    wec.calcPowerGen(waveData, [], mi.simTime, 0, 0);
                    wec.lowPowerGen = 0.75*wec.meanPowerGen; 
                    wec.powerGenMeans = ones(size(mi.simTime)) * wec.meanPowerGen;

                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;

                case 6 % Power matrix & Hs, Te time series
                    % load data
                    pMat = load(mi.resourceDataVars.dataFiles{depVarCount}, mi.resourceDataVars.pMatVarName);  pMat = pMat.pMat;
                    dataTime = load(mi.resourceDataVars.dataFiles{depVarCount}, mi.resourceDataVars.tVarName);  dataTime = dataTime.dataTime;
                    Hs = load(mi.resourceDataVars.dataFiles{depVarCount}, mi.resourceDataVars.HsVarName);  Hs = Hs.Hs;
                    Te = load(mi.resourceDataVars.dataFiles{depVarCount}, mi.resourceDataVars.TeVarName);  Te = Te.Te;
                     
                    % Calculate power from grid interpolation
                    pMatTableFn = griddedInterpolant( repmat(pMat(2:end,1), 1, size(pMat, 2)-1), repmat(pMat(1,2:end), size(pMat, 1)-1, 1), pMat(2:end, 2:end),'linear','nearest');  % X1 grid, X2 grid, valueMatrix, interpolation method, extrapolation method
                    pGen_dataTime = pMatTableFn(Hs, Te); 

                     % reshape power gen
                    wec.reshapePowerGen(pGen_dataTime, dataTime, mi.simHrs, mi.dt);  %%%%%%%%%%% time must be in hr
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;
            end
        end  % calc power gen fn
    end  % methods
end