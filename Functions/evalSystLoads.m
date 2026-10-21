% Ama Hartman - 2026

function [maxHotelLoad, maxChargingLoad] = evalSystLoads(modIn, auv, fleetSize, energyStorage)
% evalSystLoads evaluates the maximum operational loads of a given system 
% on the central battery. 
% 
% INPUTS: 
% - modIn:  ModelInput object containing 'simGoal' and 'incorpStagger'
%           values. 
% - auv:    AUV object array
% - fleetSize:      Scalar size of AUV fleet in simulation
% - energyStorage:  EnergyStorage object
% 
% OUTPUTS: 
% - maxHotelLoad:    Maximum aggregate hotel load for the given system
% - maxChargingLoad: Maximum load on the system's central battery during
%                    standard operation, not accounting for any power input

switch modIn.simGoal

    % calc. AUV loads
    case 1
        % fleetSize # of identical 'auv' AUVs
        auvMaxHotel = fleetSize * (auv.hotelLoad) / (auv.n_powerTransfer * auv.n_battery * energyStorage.n_battery);

        if modIn.incorpStagger == 1 
            numSimulRecharge = ceil(auv.chargeTime * fleetSize / (auv.missionTime + auv.chargeTime));  % ceil(t_recharge / t_stagger), t_stagger = (t_mission + t_recharge) / N_auv
        else
            numSimulRecharge = fleetSize;
        end
        auvMaxChargingLoad = numSimulRecharge * (auv.chargeRate + auv.hotelLoad) / (auv.n_powerTransfer * auv.n_battery * energyStorage.n_battery);

    case 2
        % fleetSize = numel(auv) # of DIFFERENT AUVs ...
        auvMaxHotel = sum( [auv.hotelLoad] ./ ([auv.n_powerTransfer] .* [auv.n_battery] .* energyStorage.n_battery) );
        
        auvChargingLoads = ([auv.chargeRate] + [auv.hotelLoad]) ./ ([auv.n_powerTransfer] .* [auv.n_battery] .* energyStorage.n_battery);
        if modIn.incorpStagger == 1
            estStaggerTime = min( ([auv.missionTime]+[auv.chargeTime]) ./ length(auv) );
            numSimulRecharge = ceil(min([auv.chargeTime]) / estStaggerTime);  % ceil(t_recharge(smallest) / t_stagger), t_stagger = min((t_mission + t_recharge) / N_auv)

            auvMaxChargingLoad = max( movsum(auvChargingLoads, numSimulRecharge, "Endpoints", "Discard") );  % max load from moving sum
        else
            auvMaxChargingLoad = sum(auvChargingLoads);
        end
end

esHotel = (energyStorage.baseHotelLoad + energyStorage.dockHotelLoad*fleetSize) / (energyStorage.n_battery * energyStorage.n_powerTrnsfr);

maxHotelLoad = esHotel + auvMaxHotel;
maxChargingLoad = esHotel + auvMaxChargingLoad; 

end