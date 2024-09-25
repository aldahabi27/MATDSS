function MATDSSManageSimDataApp_LoadSimDataFile(SimDataApp, MATDSSApp)
% MATDSSManageSimDataApp_LoadSimDataFile(SimDataApp, MATDSSApp)
% This function loads simulation data into the MATDSS application and saves
% it in a variable called 'LoadedSimData'. The MATDSS application will
% check the available data, compare it with the currently loaded tables,
% and decide whether to use the existing data or regenerate new data.
%
% Parameters:
%   - SimDataApp: The SimData management application instance.
%   - MATDSSApp: The MATDSS application instance.
%
% Last Update for this function was on MATDSS App Ver 0.96 (20 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% If two arguments are provided, determine which simulation data file is selected
if nargin > 1
    % Find the index of the selected simulation data file
    iSimData = find(strcmpi([SimDataApp.SimDataListBox.Items], SimDataApp.SimDataListBox.Value));

    % Load the selected simulation data into MATDSSApp
    MATDSSApp.LoadedSimData = SimDataApp.SimDataFilesData(iSimData);

    % Save the full address of the loaded simulation data file
    MATDSSApp.LoadedSimData.FileAddress = SimDataApp.SimDataFilesFullAddress(iSimData);
else
    % If only one argument is provided, use SimDataApp as MATDSSApp
    MATDSSApp = SimDataApp;
end

% Enable buttons to allow further interaction after loading simulation data
MATDSSApp.LoadConfigFromSavedSimDataButton.Enable = "on";
MATDSSApp.UseSimDataButton.Enable = "on";
MATDSSApp.ReloadSimDataButton.Enable = "on";

% Pre-load tables in the MATDSS application using the loaded simulation data
confirmation = 'continue';

% Check if the loaded simulation data corresponds to a different OpenDSS file
if ~strcmpi(MATDSSApp.OpenDSSFilesListBox.Value, MATDSSApp.LoadedSimData.Name)
    % Prompt the user to confirm if they want to continue with the mismatched file
    confirmation = questdlg('Loaded SimData are for a different OpenDSS filename. Do you want to continue?', ...
        'Confirmation', 'Continue', 'Cancel', 'Cancel');
end

% Proceed only if the user confirms they want to continue
if strcmpi(confirmation, 'continue')
    % Loop through each tab in the CircuitConfigurationsTabGroup and load simulation data
    for i = 1:length(MATDSSApp.CircuitConfigurationsTabGroup.Children)
        MATDSSApp.CircuitConfigurationsTabGroup.SelectedTab = MATDSSApp.CircuitConfigurationsTabGroup.Children(i);

        % Load the simulation data for the selected tab
        MATDSSApp_LoadSimData(MATDSSApp, 1);
    end

    % Store the current data of the ControlAreasTable
    MATDSSApp.ControlAreasTableData = MATDSSApp.ControlAreasTable.Data;

    % Load control parameters into the app from the loaded simulation data
    MATDSSApp.GeneralGainEditField.Value = MATDSSApp.LoadedSimData.ControlParameters.GeneralGain;
    MATDSSApp.RoUEditField.Value = MATDSSApp.LoadedSimData.ControlParameters.RoU;
    MATDSSApp.LoopsizeEditField.Value = MATDSSApp.LoadedSimData.ControlParameters.LoopSize;

    % Load time settings from the simulation data
    MATDSSApp.DurationEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.Duration;
    MATDSSApp.SimTimeStepEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.SimTimeStep;
    MATDSSApp.ContTimeStepEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.ContTimeStep;
    MATDSSApp.MeasTimeStepEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.MeasTimeStep;
    MATDSSApp.MeasDelayEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.MeasDelay;
    MATDSSApp.StabilizationTimeEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.StabilizationTime;

    % Load P0Set function settings
    iP0Set = MATDSS_StrComp({MATDSSApp.P0SetFunctionButtonGroup.Children.Text}, MATDSSApp.LoadedSimData.P0SetFunctionSettings.Type);
    MATDSSApp.P0SetFunctionButtonGroup.SelectedObject = MATDSSApp.P0SetFunctionButtonGroup.Children(iP0Set);
    MATDSSApp.ChangekWEditField.Value = MATDSSApp.LoadedSimData.P0SetFunctionSettings.kWVal;
    MATDSSApp.TrackP0Button.Value = MATDSSApp.LoadedSimData.P0SetFunctionSettings.TrackP0Toggle;

    % Introduce a brief pause to allow UI updates
    pause(2);
end
end


%% Old function code
%
%{

function MATDSSManageSimDataApp_LoadSimDataFile(SimDataApp,MATDSSApp)
% This function Loads simulation data into MATDSSApplication and saves them
% in a variable called LoadedSimData. MATDSS application will check the
% available data and compare it to the loaded tables and then decide to
% either use or regenerate this data.

if nargin > 1
    % Which data file is selected?
    iSimData = find(strcmpi([SimDataApp.SimDataListBox.Items],SimDataApp.SimDataListBox.Value));

    MATDSSApp.LoadedSimData = SimDataApp.SimDataFilesData(iSimData);
    MATDSSApp.LoadedSimData.FileAddress = SimDataApp.SimDataFilesFullAddress(iSimData);

else
    MATDSSApp = SimDataApp;
end


MATDSSApp.LoadConfigFromSavedSimDataButton.Enable = "on";
MATDSSApp.UseSimDataButton.Enable = "on";
MATDSSApp.ReloadSimDataButton.Enable = "on";


% Pre-load tables in MATDSSApplication with SimData (always confirm if the
% loaded data is for different OpenDSS file
confirmation = 'continue';
if ~strcmpi(MATDSSApp.OpenDSSFilesListBox.Value, MATDSSApp.LoadedSimData.Name) % Check if loaded sim data are for different opendss file
    confirmation = questdlg('Loaded SimData are for different OpenDSS filename. Do you want to continue?', ...
        'Confirmation', 'Continue', 'Cancel', 'Cancel');
end

if strcmpi(confirmation, 'continue')
    % First load tables in MATDSSApplication
    for i = 1:length(MATDSSApp.CircuitConfigurationsTabGroup.Children)
        MATDSSApp.CircuitConfigurationsTabGroup.SelectedTab = MATDSSApp.CircuitConfigurationsTabGroup.Children(i);
        MATDSSApp_LoadSimData(MATDSSApp,1);
    end
    MATDSSApp.ControlAreasTableData = MATDSSApp.ControlAreasTable.Data;

    % Loading all control parameters data
    MATDSSApp.GeneralGainEditField.Value = MATDSSApp.LoadedSimData.ControlParameters.GeneralGain;
    MATDSSApp.RoUEditField.Value = MATDSSApp.LoadedSimData.ControlParameters.RoU;
    MATDSSApp.LoopsizeEditField.Value = MATDSSApp.LoadedSimData.ControlParameters.LoopSize;

    % Loading Time Settings
    MATDSSApp.DurationEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.Duration;
    MATDSSApp.SimTimeStepEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.SimTimeStep;
    MATDSSApp.ContTimeStepEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.ContTimeStep;
    MATDSSApp.MeasTimeStepEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.MeasTimeStep;
    MATDSSApp.MeasDelayEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.MeasDelay;
    MATDSSApp.StabilizationTimeEditField.Value = MATDSSApp.LoadedSimData.TimeSettings.StabilizationTime;


    % Loading Simulation Settings
    iP0Set = MATDSS_StrComp({MATDSSApp.P0SetFunctionButtonGroup.Children.Text},MATDSSApp.LoadedSimData.P0SetFunctionSettings.Type);
    MATDSSApp.P0SetFunctionButtonGroup.SelectedObject = MATDSSApp.P0SetFunctionButtonGroup.Children(iP0Set);
    MATDSSApp.ChangekWEditField.Value = MATDSSApp.LoadedSimData.P0SetFunctionSettings.kWVal;
    MATDSSApp.TrackP0Button.Value = MATDSSApp.LoadedSimData.P0SetFunctionSettings.TrackP0Toggle;
    pause(2);
end
end


%}