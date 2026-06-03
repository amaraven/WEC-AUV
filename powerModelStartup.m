% Ama Hartman

function powerModelStartup
    % Add subfolders to path
    addpath('Component Classes');
    addpath('Functions');
    addpath('InputData');
    addpath('GUI');

    % run GUI
    powerModelGUI;
end