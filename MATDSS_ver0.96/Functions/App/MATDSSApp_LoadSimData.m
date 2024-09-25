function MATDSSApp_LoadSimData(app, Flag_DontAskForConfirmation)
% MATDSSApp_LoadSimData(app, Flag_DontAskForConfirmation)
% This function loads the tables/lists from the loaded simulation data
% handle into the MATDSS application. It checks if the loaded data is
% associated with the currently selected OpenDSS file and prompts
% the user for confirmation if necessary.
%
% Parameters:
%   - app: The MATDSS application instance.
%   - Flag_DontAskForConfirmation: A flag to skip confirmation for
%     loading data from a different OpenDSS file (default is 0).
%
% Last Update for this function was on MATDSS App Ver 0.96 (20 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

if nargin < 2
    Flag_DontAskForConfirmation = 0; % Default value for confirmation flag if not provided
end

% Initialize confirmation variable
confirmation = 'continue';

% Check if LoadedSimData is non-empty
if isempty(app.LoadedSimData)
    confirmation = 'cancel'; % If no data is loaded, cancel operation
elseif ~strcmpi(app.OpenDSSFilesListBox.Value, app.LoadedSimData.Name) && ~Flag_DontAskForConfirmation
    % If the loaded simulation data is from a different OpenDSS file
    confirmation = questdlg('Loaded SimData are for different OpenDSS filename. Do you want to continue?', ...
        'Confirmation', 'Continue', 'Cancel', 'Cancel'); % Prompt user for confirmation
end

% If the user confirms loading the data
if strcmpi(confirmation, 'continue')
    selectedTab = app.CircuitConfigurationsTabGroup.SelectedTab; % Get currently selected tab

    % Check which tab is selected and update tables accordingly
    switch selectedTab.Title
        case 'DERs'
            % Load DER data into the corresponding table
            app.DERsTable.Data = app.LoadedSimData.Table.DER;
            app.DERsTable.ColumnName = app.LoadedSimData.Table.DER.Properties.VariableNames; % Set column names
            app.DERsTable.ColumnWidth = 'auto'; % Auto-adjust column width
            app.DERsTab.Tag = 'Loaded'; % Mark tab as loaded
        case 'V & I'
            % Load voltage and current tracking data
            app.VTrackingListBox.Items = app.LoadedSimData.MATDSS.Sim.Meas.AllNodesNames(3 + 1:end); % Load all node names for voltage tracking
            app.ITrackingListBox.Items = app.LoadedSimData.MATDSS.Sim.Meas.AllBranchesPhasesNames; % Load all branch phase names for current tracking

            % Set the selected values based on loaded data
            app.VTrackingListBox.Value = app.VTrackingListBox.Items(MATDSS_StrComp(app.VTrackingListBox.Items, app.LoadedSimData.Table.VTrackingListBoxValues));
            app.ITrackingListBox.Value = app.ITrackingListBox.Items(MATDSS_StrComp(app.ITrackingListBox.Items, app.LoadedSimData.Table.ITrackingListBoxValues));
        case 'Control Areas'
            % Load control area settings data into respective tables
            app.ControlAreasTable.Data = app.LoadedSimData.Table.CASettings; % Load control areas data
            app.ControlAreasTable.ColumnName = app.LoadedSimData.Table.CASettings.Properties.VariableNames; % Set column names
            app.ControlAreasTable.ColumnWidth = 'auto'; % Auto-adjust column width

            % Load controllers settings data
            app.ControllersSettingsTable.Data = app.LoadedSimData.Table.ControllersSettings;
            app.ControllersSettingsTable.ColumnName = app.LoadedSimData.Table.ControllersSettings.Properties.VariableNames; % Set column names
            app.ControllersSettingsTable.ColumnWidth = 'auto'; % Auto-adjust column width

            % Load VDER settings data
            app.VDERsTable.Data = app.LoadedSimData.Table.VDERSettings;
            app.VDERsTable.ColumnName = app.LoadedSimData.Table.VDERSettings.Properties.VariableNames; % Set column names
            app.VDERsTable.ColumnWidth = 'auto'; % Auto-adjust column width

            % Load disturbance data
            app.DisturbanceTable.Data = app.LoadedSimData.Table.Disturbance;
            app.DisturbanceTable.ColumnName = app.LoadedSimData.Table.Disturbance.Properties.VariableNames; % Set column names
    end
end
end


%% Old function code
%{

function MATDSSApp_LoadSimData(app,Flag_DontAskForConfirmation)
% This function will load the tables/lists from loaded simData handle to
% MATDSSApplication

if nargin < 2
    Flag_DontAskForConfirmation = 0;
end

% Check if LoadedSimData is none-empty
confirmation = 'continue';
if isempty(app.LoadedSimData) % Check if LoadedSimData is none-empty
    confirmation = 'cancel';
elseif ~strcmpi(app.OpenDSSFilesListBox.Value, app.LoadedSimData.Name) && ~Flag_DontAskForConfirmation % Check if loaded sim data are for different opendss file
    confirmation = questdlg('Loaded SimData are for different OpenDSS filename. Do you want to continue?', ...
        'Confirmation', 'Continue', 'Cancel', 'Cancel');
end

% If user wants to load it anyway, then just load it!
if strcmpi(confirmation,'continue')
    selectedTab = app.CircuitConfigurationsTabGroup.SelectedTab;

    % check which tab is selected and update tables
    % accordingly.
    switch selectedTab.Title
        case 'DERs'
            app.DERsTable.Data = app.LoadedSimData.Table.DER;
            app.DERsTable.ColumnName = app.LoadedSimData.Table.DER.Properties.VariableNames;
            app.DERsTable.ColumnWidth = 'auto';
            app.DERsTab.Tag = 'Loaded';
        case 'V & I'
            app.VTrackingListBox.Items = app.LoadedSimData.MATDSS.Sim.Meas.AllNodesNames(3+1:end);
            app.ITrackingListBox.Items = app.LoadedSimData.MATDSS.Sim.Meas.AllBranchesPhasesNames;
            app.VTrackingListBox.Value = app.VTrackingListBox.Items(MATDSS_StrComp(app.VTrackingListBox.Items,app.LoadedSimData.Table.VTrackingListBoxValues));
            app.ITrackingListBox.Value = app.ITrackingListBox.Items(MATDSS_StrComp(app.ITrackingListBox.Items,app.LoadedSimData.Table.ITrackingListBoxValues));
        case 'Control Areas'
            app.ControlAreasTable.Data = app.LoadedSimData.Table.CASettings;
            app.ControlAreasTable.ColumnName = app.LoadedSimData.Table.CASettings.Properties.VariableNames;
            app.ControlAreasTable.ColumnWidth = 'auto';

            app.ControllersSettingsTable.Data = app.LoadedSimData.Table.ControllersSettings;
            app.ControllersSettingsTable.ColumnName = app.LoadedSimData.Table.ControllersSettings.Properties.VariableNames;
            app.ControllersSettingsTable.ColumnWidth = 'auto';

            app.VDERsTable.Data = app.LoadedSimData.Table.VDERSettings;
            app.VDERsTable.ColumnName = app.LoadedSimData.Table.VDERSettings.Properties.VariableNames;
            app.VDERsTable.ColumnWidth = 'auto';
            
            app.DisturbanceTable.Data = app.LoadedSimData.Table.Disturbance;
            app.DisturbanceTable.ColumnName = app.LoadedSimData.Table.Disturbance.Properties.VariableNames;
    end
end
end

%}