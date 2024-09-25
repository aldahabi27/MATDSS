function MATDSSApp_Defaults(app)
% MATDSSApp_Defaults(app)
% This function resets the MATDSS Application (Main Window) GUI to
% its default values/settings. It will re-enable all buttons that are
% disabled and reset timing parameters.
%
% Last Update for this function was on MATDSS App Ver 0.96
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Add necessary paths for MATDSS functions and files
addpath DSS_Files;
addpath Functions;
addpath Functions\App;
addpath Functions\MATDSS;
addpath Exports;
addpath Config; % Location where parameters configurations can be specified (defining DERs, setting L2C parameters, etc.)
addpath Config\SimSetupMats;

% Update application status to indicate initialization
MATDSSApp_Status(app, 'Initializing (Checking MATDSS Circuit Files)', 0.5);
pause(0.001);

% Get list of *.dss files to select from in the drop-down list
FilesList = dir([pwd '\DSS_Files\MATDSS_*.dss']); % Look for *.dss files starting with "MATDSS_"
app.OpenDSSFilesListBox.Tag = [pwd '\DSS_Files\']; % Save folder location in the tag for easy access

% Update the file list in the GUI
FilesListNames = {FilesList.name};
app.OpenDSSFilesListBox.Items = FilesListNames;

% Reset the UI
MATDSSApp_Status(app, 'Initializing (reset the UI)', 0.75);
pause(0.001);

% Set default timing values
app.DurationEditField.Value = '10'; % Duration in seconds
app.SimTimeStepEditField.Value = '10'; % Time step for simulation in milliseconds
app.ContTimeStepEditField.Value = '100'; % Control time step in milliseconds
app.MeasTimeStepEditField.Value = app.ContTimeStepEditField.Value; % Measurement time step in milliseconds
app.MeasDelayEditField.Value = '0'; % Measurement delay in milliseconds
app.StabilizationTimeEditField.Value = '0'; % Stabilization time in seconds

% Configure P0Set function buttons and text according to selected option
MATDSSApp_P0SetFunButtonGroup(app);

% Reset Progress Bar settings and enable all buttons and text boxes
app.ProgressBar1.IconAlignment = 'center';
WaitBar1 = permute(repmat(app.ProgressBar1.BackgroundColor, 150, 1, 1500), [1, 3, 2]);
app.ProgressBar1.Icon = WaitBar1;
% MATDSSApp_GUIChanges(app, 'Enable_GUI'); % Reset GUI

% Display the first message as if the app has just started
MATDSSApp_Details(app, {['MATDSS Application V. ' app.AppVer]; ['Program is Initializing, Date: ' datestr(datetime('now'))]; '******************************************************'; ''});
app.MATDSSApplicationUIFigure.Tag = app.AppVer; % Reset the tag

% Check if there are any OpenDSS files and update the list
pause(0.5);
MATDSSApp_Status(app, 'Please wait!', 0.99);
pause(0.5);
if ~isempty(FilesListNames)
    MATDSSApp_OpenDSSFilesListBoxValueChangedFcn(app); % Update options in the Plot Properties Panel
end

% Initialize the Disturbance Table if not already done
if strcmpi(app.DisturbanceTable.Tag, 'uninitialized')
    DisturbanceTableHeaders = {'ID', [char(39) 'Bus name' char(39)], 'Phases ([1,2,3])', 'nPhases', ['Conn. Type (' char(39) 'Y' char(39) ' or ' char(39) 'D' char(39) ')'], 'P (kW)', 'PF', '[t1,t2] (s)', 'Enabled'};
    DisturbanceTableVarTypes = repmat({'string'}, 1, length(DisturbanceTableHeaders));
    DisturbanceTable = table('Size', [0, length(DisturbanceTableHeaders)], 'VariableTypes', DisturbanceTableVarTypes, 'VariableNames', DisturbanceTableHeaders);
    app.DisturbanceTable.Data = DisturbanceTable;
    app.DisturbanceTable.ColumnName = DisturbanceTableHeaders;
    app.DisturbanceTable.Tag = 'initialized';
    app.DeleteDisturbanceButton.Enable = false;
end

% Add default disturbance values if the table is empty
if ~size(app.DisturbanceTable.Data, 1)
    i = size(app.DisturbanceTable.Data, 1) + 1;
    NewDisturbanceTableValues = {num2str(i), "''", "[1, 2, 3]", "3", "Y", "20", "0.9", "[1, inf]", "T"};
    app.DisturbanceTable.Data = [app.DisturbanceTable.Data; NewDisturbanceTableValues];
    app.DisturbanceTableData = app.DisturbanceTable.Data;
end




% Configure the application based on the selected tab
MATDSSApp_ConfigurationTabSelect(app);

end