% modelStartup adds all relevant subfolders and files to the working path
% and opens the model GUI. Run to initiate. 
% 
% Ama Hartman

function modelStartup
    % Add subfolders to path
    addpath('Component Classes');
    addpath('Functions');
    addpath('InputData');
    addpath('GUI');

    % run GUI
    powerModelGUI;
end