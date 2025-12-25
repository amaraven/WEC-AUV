% Ama Hartman

classdef ModelInput < handle
    % UserInputs defines the properties and methods saved in UserInputs
    % objects. These objects store user-defined or default inputs for the
    % 'powerModel.m' simulation. 
    %
    % User-Defined properties: 
    % * simHrs: [hr] simulation duration. Default 1 wk
    % * depVar: variable to change between simulation iterations. Default
    %   'AUV Model'
    % * resourceDataType: Type of resource data. Power in W, time in h,
    %   distance in m, wave periods in s.
    %   1. Modeled power gen. of 1.5 m WEC in different sea states (1-10)
    %   2. Power generation time series and time vector (*.mat file) 
    %   3. Value of mean power generation
    %   4. Time series of wave specifications (significant wave
    %   height, wave energy and peak periods, and time) (Default)
    %   5. Values of mean wave specifications 
    %   6. Power matrix and wave spec. time series
    % * incorpStagger: Enforce a stagger between AUV deployments. Default 1
    %   (on)
    % * maxFleetSize: Enforce a maximum AUV fleet size. Default 0 (off)
    % * userDefinedBattery: Size of user-defined central battery capacity. 
    %   Default 0 (off)
    % * dtSec: timestep interval in seconds. Default 30 s

    % User-defined properties
    properties
        simHrs (1,1) {mustBePositive} = (7*24)  % default: 1 wk
        depVar (1,1) string {mustBeMember(depVar, ["AUV Model", "WEC Power Gen / Wave Resource"])} = "AUV Model"  % default: compare AUV models
        resourceDataType (1,1) {mustBeMember(resourceDataType, [1,2,3,4,5,6])} = 4  % default: 4 (time series of wave specs)
        incorpStagger (1,1) {mustBeMember(incorpStagger, [1,0])} = 1  % default: 1 (yes incorporate stagger)
        maxFleetSize (1,1) {mustBeNonnegative} = 0  % default: 0 (no max fleet size)
        userDefinedBattery (1,1) {mustBeNonnegative} = 0  % default: 0 (no user-defined battery)
        dtSec (1,1) {mustBeNonnegative} = 30  % default: 30-second timestep intervals
        auvModels (1,:) cell = [{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}, {'G'}, {'H'}, {'I'}, {'J'}, {'K'}, {'L'}, {'M'}, {'N'}, {'O'}, {'P'}, {'Q'}, {'R'}, {'S'}, {'T'}, {'U'}] % Cell containing all auvModels to use in simulations. Models are defined in 'AUV' class script.
    end

    % Internal properties
    properties (Dependent = true)  % properties dependent on other props
        dt % [hr]
        simTime % [hr]
    end

    properties (Hidden = true, Access = private)  % Hidden from property lists and property accessed only by 'ModelInput' class members
        rdVars struct = struct()

    end


    %% Instance Methods (need object as an input)
    methods
        %% Constructor
        function mi = ModelInput(simHrs, depVar, resourceDataType, incorpStagger, maxFleetSize, userDefinedBattery, dtSec, auvModels)
            arguments
                simHrs (1,1) {mustBePositive} = (7*24)  % default: 1 wk
                depVar (1,1) string {mustBeMember(depVar, ["AUV Model", "WEC Power Gen / Wave Resource"])} = "AUV Model"  % default: compare AUV models
                resourceDataType (1,1) {mustBeMember(resourceDataType, [1,2,3,4,5,6])} = 4  % default: 4 (time series of wave specs)
                incorpStagger (1,1) {mustBeMember(incorpStagger, [1,0])} = 1  % default: 1 (yes incorporate stagger)
                maxFleetSize (1,1) {mustBeNonnegative} = 0  % default: 0 (no max fleet size)
                userDefinedBattery (1,1) {mustBeNonnegative} = 0  % default: 0 (no user-defined battery)
                dtSec (1,1) {mustBeNonnegative} = 30  % default: 30-second timestep intervals
                auvModels (1,:) cell = [{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}, {'G'}, {'H'}, {'I'}, {'J'}, {'K'}, {'L'}, {'M'}, {'N'}, {'O'}, {'P'}, {'Q'}, {'R'}, {'S'}, {'T'}, {'U'}]
            end

            mi.simHrs = simHrs;
            mi.depVar = depVar;
            mi.resourceDataType = resourceDataType; 
            mi.incorpStagger = incorpStagger;
            mi.maxFleetSize = maxFleetSize;
            mi.userDefinedBattery = userDefinedBattery;
            mi.dtSec = dtSec; 
            mi.auvModels = auvModels;
           
        end


        %% Dependent property calcs
        function dt = get.dt(mi)
            dt = mi.dtSec /60 /60;  %[hr]
        end

        function simTime = get.simTime(mi)
            simSec = mi.simHrs *60 *60;  % [s]
            simTime = (mi.dtSec:mi.dtSec:simSec)' /60 /60;  % [hr]
        end



        %% Load Resource Data, convert to powerGen time series
        function extractDataInfo(mi)
            % Given user-inputs, extract metainformation needed to load and 
            % process input data

            switch mi.resourceDataType
                case 1  % Proteus struct containing power generated during different sea states
                    if strcmp(mi.depVar, 'AUV Model')
                        mi.rdVars.seaState = listdlg('PromptString',{'Specify the sea state you would like to use for this simulation'}, 'SelectionMode','single', ...
                                'ListString',{'1','2','3','4','5','6','7','8','9','10'});
                    end

                case 2  % Time series of power gen
                    switch mi.depVar
                        case 'AUV Model'
                            mi.rdVars.dataFiles = uigetfile('*.mat', 'Select *.mat data files','MultiSelect','on');
                            fileVars = whos('-file', mi.rdVars.dataFiles); 

                        case 'WEC Power Gen / Wave Resource' 
                            mi.rdVars.dataFiles = uigetfile('*.mat', 'Select *.mat data files','MultiSelect','on');
                            if ischar(mi.rdVars.dataFiles)  % if only one is selected
                                mi.rdVars.dataFiles = {mi.rdVars.dataFiles};
                            end
            
                            fileVars = whos('-file', mi.rdVars.dataFiles{1});
                    end
                    tIndx = listdlg('PromptString',{'Select the time series vector'}, 'SelectionMode','single', ...
                            'ListString',{fileVars.name});
                    pwrIndx = listdlg('PromptString',{'Select the power gen. time series'}, 'SelectionMode','single', ...
                            'ListString',{fileVars.name});

                    mi.rdVars.tVarName = fileVars(tIndx).name;
                    mi.rdVars.pwrVarName = fileVars(pwrIndx).name;

                case 3  % Value of mean power
                    mi.rdVars.meanPower = input('Enter mean power generated as a single value or a bracketed, comma-seperated list: ');

                case 4 % Time series of wave specs (Hs, Te, Tp)
                    switch mi.depVar
                        case 'AUV Model'
                            mi.rdVars.dataFiles = uigetfile('*.csv', 'Select *.csv Data File', 'MultiSelect','off');
                            mi.rdVars.dataFiles = {mi.rdVars.dataFiles};
                            
                        case 'WEC Power Gen / Wave Resource'
                            mi.rdVars.dataFiles = uigetfile('*.csv', 'Select *.csv Data Files', 'MultiSelect','on');
                            if ischar(mi.rdVars.dataFiles)
                                mi.rdVars.dataFiles = {mi.rdVars.dataFiles};
                            end
                    end

                case 5 % Mean Wave specs (Hs, Te, Tp)
                    % examples:
                    % sigWaveHeight = [1.8774    1.8800    1.9000    2.2000    1.6000]';
                    % energyPeriod = [8.9367    8.9400    9.0000    9.0000    8.8000]';
                    % peakPeriod = [13.1428   13.1400   13.1400   13.1400   10.1000]';
                    prompt = {'Enter significant wave height(s) [m]:', 'Enter wave energy period(s) [s]:','Enter peak period(s) [s]:'};
                    fieldsize = [1, 45; 1, 45; 1, 45];
                    switch mi.depVar
                        case 'AUV Model'
                            defaultInputs = {'1.8774','8.9367','13.1428'};
                        case 'WEC Power Gen / Wave Resource'
                            defaultInputs = {'[1.8774, 1.8800, 1.9000, 2.2000, 1.6000]','[8.9367, 8.9400, 9.0000, 9.0000, 8.8000]','[13.1428, 13.1400, 13.1400, 13.1400, 10.1000]'};
                    end
                    tryAgain = 1;  % preallocate
                    while tryAgain == 1
                        waveSpecDlg = inputdlg(prompt, 'Wave Resource Specifications', fieldsize, defaultInputs);
                        sigWaveHeight = str2num(waveSpecDlg{1});
                        energyPeriod = str2num(waveSpecDlg{2});
                        peakPeriod = str2num(waveSpecDlg{3});
            
                        if isequal( size(sigWaveHeight), size(energyPeriod), size(peakPeriod) )
                            tryAgain = 0;
                        else
                            warning('Each specification must have the same number of values. Please re-enter the wave resource specifications.')
                        end
                    end
                    
                    mi.rdVars.waveSpecTable = table(sigWaveHeight', energyPeriod', peakPeriod','VariableNames',{'SignificantWaveHeight','EnergyPeriod','PeakPeriod'});

                case 6 % Power matrix & Hs, Te time series
                    switch mi.depVar
                        case 'AUV Model'
                            mi.rdVars.dataFiles = uigetfile('*.mat', 'Select *.mat data file for time series','MultiSelect','off');
                            fileVars = whos('-file', mi.rdVars.dataFiles);
                            mi.rdVars.dataFiles = {mi.rdVars.dataFiles};
                        case 'WEC Power Gen / Wave Resource' 
                            mi.rdVars.dataFiles = uigetfile('*.mat', 'Select *.mat data files','MultiSelect','on');
                            if ischar(mi.rdVars.dataFiles)  % if only one is selected
                                mi.rdVars.dataFiles = {mi.rdVars.dataFiles};
                            end                
                            fileVars = whos('-file', mi.rdVars.dataFiles{1});
                    end
                    

                    pMatIndx = listdlg('PromptString',{'Select the power matrix variable'}, 'SelectionMode','single', ...
                            'ListString',{fileVars.name});
                    tIndx = listdlg('PromptString',{'Select the time series vector'}, 'SelectionMode','single', ...
                            'ListString',{fileVars.name});
                    HsIndx = listdlg('PromptString',{'Select the significant wave height time series'}, 'SelectionMode','single', ...
                            'ListString',{fileVars.name});
                    TeIndx = listdlg('PromptString',{'Select the wave energy period time series'}, 'SelectionMode','single', ...
                            'ListString',{fileVars.name});

                    mi.rdVars.pMatVarName = fileVars(pMatIndx).name;
                    mi.rdVars.tVarName = fileVars(tIndx).name;
                    mi.rdVars.HsVarName = fileVars(HsIndx).name;
                    mi.rdVars.TeVarName = fileVars(TeIndx).name;
            end
        end

        
        function calcPowerGen(mi, wec, modOut, depVarCount)
            % From user-input data, calculates power generation vector. 
            % Inputs: 
            % * wec: WEC class object
            % * modOut: ModelOutput class object
            % * depVarCount: (Optional) iteration of the dependent variable
            %   loop, exclude if calling from outside of the loop
            arguments
                mi ModelInput
                wec WEC
                modOut ModelOutput
                depVarCount (1,1) {mustBePositive, mustBeInteger} = 1
            end

            % If metadata has not yet been loaded..
            if isempty(mi.rdVars)
                    warning('Must load data first. Loading data...')
                    mi.extractDataInfo;
            end

            % Calculate and save power gen
            switch mi.resourceDataType
                case 1 % Proteus struct containing power generated during different sea states
                    % load data
                    load('inputData/RM3_seaState_Power.mat', 'RM3');  % Time series output of proteus WEC model for a single point absorber with no dock attached.
                    
                    % extract seaState to use for power calc.
                    switch mi.depVar
                        case 'AUV Model'
                            seaState = mi.rdVars.seaState;
                        case 'WEC Power Gen / Wave Resource'
                            seaState = depVarCount;
                    end
                    
                    % calc
                    wec.reshapePowerGen(RM3(seaState).Power, RM3(seaState).Time, mi.simHrs, mi.dt);  
                    
                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;  % save power gen
                    modOut.seaState(depVarCount) = seaState;  % save sea state

                case 2 % Time series of power gen
                    % load data
                    dataTime = load(mi.rdVars.dataFiles{depVarCount}, mi.rdVars.tVarName);
                    pGen_dataTime = load(mi.rdVars.dataFiles{depVarCount}, mi.rdVars.pwrVarName);

                    % reshape power gen. to fit simulation timestep
                    wec.reshapePowerGen(pGen_dataTime, dataTime, mi.simHrs, mi.dt);
                     
                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;
                    if isempty(modOut.dataIn)
                        modOut.dataIn = mi.rdVars.dataFiles;
                    end

                case 3  % Value of mean power
                    % save
                    if isempty(modOut.meanPowerGen)
                        modOut.meanPowerGen = mi.rdVars.meanPower;
                    end
    
                    % save for use in current sim. iteration
                    wec.meanPowerGen = modOut.meanPowerGen(depVarCount);
                    wec.lowPowerGen = 0.75 * wec.meanPowerGen; 
                    wec.powerGenMeans = ones(size(mi.simTime)) * wec.meanPowerGen; 

                case 4 % Time series of wave specs (Hs, Te, Tp)
                    % load data
                    mi.rdVars.dataTable = readtable(mi.rdVars.dataFiles{depVarCount},'VariableNamingRule','modify');

                    % calculate
                    wec.calcPowerGen(mi.rdVars.dataTable,'meanPwr', mi.simTime, 0, 0); %% set windowOverrideIndx to 403 (and use Oregon dataset) to replicate paper results  

                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;
                    if isempty(modOut.dataIn)
                        modOut.dataIn = mi.rdVars.dataFiles; 
                    end

                case 5 % Mean Wave specs (Hs, Te, Tp)
                    % calculate
                    wec.calcPowerGen(mi.rdVars.waveSpecTable(depVarCount,:), [], mi.simTime, 0, 0);
                    wec.lowPowerGen = 0.75*wec.meanPowerGen; 
                    wec.powerGenMeans = ones(size(mi.simTime)) * wec.meanPowerGen;

                    % save
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;
                    if isempty(modOut.dataIn)
                        modOut.dataIn = mi.rdVars.waveSpecTable;
                    end

                case 6 % Power matrix & Hs, Te time series
                    % load data
                    pMat = load(mi.rdVars.dataFiles{depVarCount}, mi.rdVars.pMatVarName);  pMat = pMat.pMat;
                    dataTime = load(mi.rdVars.dataFiles{depVarCount}, mi.rdVars.tVarName);  dataTime = dataTime.dataTime;
                    Hs = load(mi.rdVars.dataFiles{depVarCount}, mi.rdVars.HsVarName);  Hs = Hs.Hs;
                    Te = load(mi.rdVars.dataFiles{depVarCount}, mi.rdVars.TeVarName);  Te = Te.Te;
                     
                    % Calculate power from grid interpolation
                    pMatTableFn = griddedInterpolant( repmat(pMat(2:end,1), 1, size(pMat, 2)-1), repmat(pMat(1,2:end), size(pMat, 1)-1, 1), pMat(2:end, 2:end),'linear','nearest');  % X1 grid, X2 grid, valueMatrix, interpolation method, extrapolation method
                    pGen_dataTime = pMatTableFn(Hs, Te); 

                     % reshape power gen
                    wec.reshapePowerGen(pGen_dataTime, dataTime, mi.simHrs, mi.dt);
                    modOut.meanPowerGen(depVarCount) = wec.meanPowerGen;
                    if isempty(modOut.dataIn)
                        modOut.dataIn = dataFiles;
                    end
            end

        end
        


    end  % methods
end