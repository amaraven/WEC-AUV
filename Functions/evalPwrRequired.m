% Ama Hartman 2026

function powerRequired = evalPwrRequired(auv, wec, energyStorage)
% evalPwrRequired calculates the mean power needed to support a WEC-AUV
% system composed of the AUV fleet saved in the auv object array and
% component specifications saved in wec and energyStorage. 
arguments
    auv AUV
    wec WEC
    energyStorage EnergyStorage
end

t_c = [auv.missionTime] + [auv.chargeTime]; % time required to complete an AUV mission + recharge cycle

pwrToAUVs = [auv.chargeLoad]./(energyStorage.n_battery*t_c);  % averaged AUV draw rate throughout mission (no draw) + recharge (draw) cycle

eStorageHotel = energyStorage.baseHotelLoad + energyStorage.dockHotelLoad*(numel(auv));

rechargePower = (sum(pwrToAUVs) + wec.hotelLoad/(wec.n_battery^2) + eStorageHotel/energyStorage.n_battery/energyStorage.n_powerTrnsfr)/(energyStorage.n_wecPwrTrnsfr*energyStorage.n_battery);  % sum(pwrToAUVs) is the same as mean(pwrToAUVs)*numel(auv)

hotelPower = (sum([auv.hotelLoad]./([auv.n_powerTransfer].*[auv.n_battery].*energyStorage.n_battery)) + wec.hotelLoad/(wec.n_battery^2) + eStorageHotel/energyStorage.n_battery/energyStorage.n_powerTrnsfr)/(energyStorage.n_wecPwrTrnsfr*energyStorage.n_battery);

powerRequired = max([rechargePower, hotelPower]);
end