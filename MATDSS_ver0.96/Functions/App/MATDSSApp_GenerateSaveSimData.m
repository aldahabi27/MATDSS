function MATDSSApp_GenerateSaveSimData(app)
% MATDSSApp_GenerateSaveSimData saves the current simulation setup data to a .mat file.
%
% The saved data includes:
%    - Text-based copies of DERs, Voltage, Current, CASettings, Controllers
%      Settings, VDERs, Control Parameters, Time settings, and Simulation
%      Settings.
%    - Properly formatted data for reloading using MATDSSApp_LoadSimData(app).
%
% This function performs the following steps:
%    1. Updates the application status and disables the UseSimData button.
%    2. Initializes and populates the SimSetupData structure with data from UI components.
%    3. Generates additional data (e.g., VI Tables) and initializes simulation parameters.
%    4. Creates DER devices and disturbance loads, and initializes system matrices.
%    5. Saves the simulation data to a .mat file and restores the UseSimData button state.
%
% 
% This function saves various simulation data into a .mat file located in the Config folder.
% The data includes configurations for DERs, voltage, current, control areas settings, 
% controllers settings, VDER settings, control parameters, time settings, and simulation settings.
% It also initializes the simulation, generates required variables and matrices, and performs 
% various operations to set up and prepare the simulation environment.
%
% Parameters:
%   - app: The application object containing UI components and data to be saved.
%
% Last Update: MATDSS App Ver 0.96 (19 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com




% Update status and briefly pause to notify the user that data generation has started
MATDSSApp_Status(app, 'Generating SimData...', 'Generating Simulation Data Started.');
pause(0.1);

% Disable the UseSimData button to prevent conflicts during data generation
State_UseSimDataButton = app.UseSimDataButton.Value;
app.UseSimDataButton.Value = 0;

% Initialize the SimSetupData structure with the name from the app
SimSetupData = struct;
SimSetupData.Name = app.MATDSSExcelSheet;

% Copy data from tables and list boxes into SimSetupData
MATDSSApp_Details(app, 'Copying Tables to SimData');
pause(0.1);

% Initialize the Table field in SimSetupData
SimSetupData.Table = struct;
SimSetupData.Table.DER = app.DERsTable.Data; % DERs Table data
SimSetupData.Table.VTrackingListBoxValues = app.VTrackingListBox.Value; % Voltage tracking list box values
SimSetupData.Table.ITrackingListBoxValues = app.ITrackingListBox.Value; % Current tracking list box values
SimSetupData.Table.VTrackingListBoxItems = app.VTrackingListBox.Items; % Voltage tracking list box items
SimSetupData.Table.ITrackingListBoxItems = app.ITrackingListBox.Items; % Current tracking list box items
SimSetupData.Table.CASettings = app.ControlAreasTable.Data; % Control Areas settings
SimSetupData.Table.ControllersSettings = app.ControllersSettingsTable.Data; % Controllers settings
SimSetupData.Table.VDERSettings = app.VDERsTable.Data; % VDERs settings
SimSetupData.Table.Disturbance = app.DisturbanceTable.Data; % Disturbance table data

% Copy table data (variables) and convert to lowercase for consistency
MATDSSApp_Details(app, 'Copying Tables-Data to SimData');
pause(0.1);

% Initialize the TableData field in SimSetupData
SimSetupData.TableData = struct;
SimSetupData.TableData.DER = app.DERsTable.Data.Variables; % DERs table variables
SimSetupData.TableData.CASettings = app.ControlAreasTable.Data.Variables; % Control Areas settings variables
SimSetupData.TableData.ControllersSettings = app.ControllersSettingsTable.Data.Variables; % Controllers settings variables
SimSetupData.TableData.VDERSettings = app.VDERsTable.Data.Variables; % VDERs settings variables
SimSetupData.TableData.Disturbance = app.DisturbanceTable.Data.Variables; % Disturbance table variables

% Convert all controller settings and control area settings to lowercase
SimSetupData.TableData.ControllersSettings = cellfun(@lower, SimSetupData.TableData.ControllersSettings, 'UniformOutput', false);
SimSetupData.TableData.CASettings = cellfun(@lower, SimSetupData.TableData.CASettings, 'UniformOutput', false);

% Generate and add VI Tables to SimSetupData
MATDSS = struct;
MATDSS.Table.DER = SimSetupData.Table.DER;
MATDSS.TableData.DER = SimSetupData.TableData.DER;
SimSetupData.Table.VI = MATDSSApp_VITableFun(app, MATDSS); % Generate VI Tables
SimSetupData.TableData.VI = SimSetupData.Table.VI.Variables; % Add VI Tables data
clear MATDSS; % Clear temporary MATDSS variable

% Copy control parameters from UI components
MATDSSApp_Details(app, 'Copying Control Parameters to SimData');
pause(0.1);

SimSetupData.ControlParameters = struct;
SimSetupData.ControlParameters.GeneralGain = app.GeneralGainEditField.Value; % General Gain value
SimSetupData.ControlParameters.RoU = app.RoUEditField.Value; % RoU value
SimSetupData.ControlParameters.LoopSize = app.LoopsizeEditField.Value; % Loop size value
SimSetupData.ControlParameters.DisturbanceTable = app.DisturbanceTable.Data; % Disturbance table data

% Copy time settings from UI components
MATDSSApp_Details(app, 'Copying Time Parameters to SimData');
pause(0.1);

SimSetupData.TimeSettings = struct;
SimSetupData.TimeSettings.Duration = app.DurationEditField.Value; % Simulation duration
SimSetupData.TimeSettings.SimTimeStep = app.SimTimeStepEditField.Value; % Simulation time step
SimSetupData.TimeSettings.ContTimeStep = app.ContTimeStepEditField.Value; % Continuous time step
SimSetupData.TimeSettings.MeasTimeStep = app.MeasTimeStepEditField.Value; % Measurement time step
SimSetupData.TimeSettings.MeasDelay = app.MeasDelayEditField.Value; % Measurement delay
SimSetupData.TimeSettings.StabilizationTime = app.StabilizationTimeEditField.Value; % Stabilization time

% Copy simulation settings from UI components
MATDSSApp_Details(app, 'Copying Simulation Settings to SimData');
pause(0.1);

SimSetupData.P0SetFunctionSettings = struct;
SimSetupData.P0SetFunctionSettings.Type = app.P0SetFunctionButtonGroup.SelectedObject.Text; % P0 Set Function type
SimSetupData.P0SetFunctionSettings.kWVal = app.ChangekWEditField.Value; % kW value
SimSetupData.P0SetFunctionSettings.TrackP0Toggle = app.TrackP0Button.Value; % Track P0 toggle

% Initialize simulation variables and matrices
GlobalGain = str2double(app.GeneralGainEditField.Value); % Convert General Gain to double

MATDSSApp_Details(app, 'Generating Simulation Matrices...');
pause(0.1);

MATDSSApp_Details(app, 'Initializing Simulation...');
pause(0.1);
[MATDSS] = MATDSSApp_Initialization(app, GlobalGain); % Initialize simulation with Global Gain

MATDSS.Table = SimSetupData.Table;
MATDSS.TableData = SimSetupData.TableData;

% Check if DSS is loaded and linked successfully
if MATDSS.Sim.DSSStartOk
    % Generate DER devices
    MATDSSApp_Details(app, 'Generating DERs');
    pause(0.1);
    [MATDSS, DER] = MATDSS_DER(MATDSS);

    % Generate disturbance loads
    MATDSSApp_Details(app, 'Generating Disturbances');
    pause(0.1);
    MATDSS = MATDSS_Disturbance(MATDSS);
    for i = 1:20
        MATDSS.Sim.DSSSolution.Solve; % Solve the DSS solution 20 times
    end

    % Generate YBus matrices
    MATDSSApp_Details(app, 'Generating Ybus matrices');
    pause(0.1);
    MATDSS = MATDSS_YBusMatrix(MATDSS, 1);

    % Initialize time index
    t = 1; % Initial time
    MATDSS.Sim.t = t;
    at = MATDSS.Time.Sim.TimeSpan(t); % Actual time in seconds
    MATDSS.Sim.at = at;

    % Initialize measurement vectors
    MATDSSApp_Details(app, 'Generating Measurement vectors');
    pause(0.1);
    MATDSS = MATDSS_Measurement(MATDSS);
    MATDSS.ControlSignals.P0Set = ones(length(MATDSS.Sim.Meas.at), 1) .* sum(MATDSS.Meas.P0(:, t));

    % Define controllers and control areas
    MATDSSApp_Details(app, 'Defining controllers (this might take a while...)');
    pause(0.1);
    [MATDSS, DER] = MATDSS_Controllers(MATDSS, DER);

    % Initialize controllers
    MATDSS.Sim.RoU = 0;
    [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS, DER);

    % Initialize DERs
    MATDSSApp_Details(app, 'Initializing DERs');
    pause(0.1);
    for i = 1:MATDSS.Sim.nDER
        [MATDSS, DER] = MATDSS_DERUpdate(MATDSS, DER, i);
    end

    % Initialize disturbance loads
    MATDSS = MATDSS_DisturbanceUpdate(MATDSS);

    % Update live measurements (commented out as it might be unnecessary for saving)
    % MATDSSApp_LiveMeasurements(MATDSS, DER);
end

% Save simulation data to .mat file
MATDSSApp_Details(app, 'Packing and Saving Simulation Data file');
MATDSS.Sim.DSSObj = [];
MATDSS.Sim.DSSText = [];
MATDSS.Sim.DSSCircuit = [];
MATDSS.Sim.DSSSolution = [];
SimSetupData.MATDSS = MATDSS;
SimSetupData.DER = DER;

% Create filename with timestamp
currenttime = datetime;
currenttime.Format = 'MMddyyHHmm';
currenttime = char(currenttime);

SimSetupData.currenttime = currenttime;

% Create directory if it does not exist
if exist(fullfile([pwd '\Config\SimSetupMats']), 'dir') == 0
    MATDSSApp_Details(app, 'SimData folder is not found, creating it now.');
    mkdir([pwd '\Config'], 'SimSetupMats');
    MATDSSApp_Details(app, 'Folder created, saving SimData file now');
end

% Save SimSetupData to .mat file
save([pwd '\Config\SimSetupMats\' app.MATDSSExcelSheet '_' currenttime '.mat'], 'SimSetupData', '-v7.3');

% Restore the UseSimData button state
app.UseSimDataButton.Value = State_UseSimDataButton;

% Update status to indicate successful export
MATDSSApp_Status(app, 'Simulation Data Exported Successfully', ['Exporting completed to ' app.MATDSSExcelSheet '_' currenttime '.mat']);
end


%% Old Function Code

%{

function MATDSSApp_GenerateSaveSimData(app)
% This function will save the following in SimSetupData.mat file in Config folder
%    * Text based copy of DERs, Voltage, Current, CASettings, Controllers
%    Settings, VDERs, Control Parameters, Time settings, Simulation
%    Settings.
%
%    * A proper type for the data to be able to reload them in the program
%    using MATDSSApp_LoadSimData(app)
MATDSSApp_Status(app,['Generating SimData...'],['Generating Simulation Data Started.']);
pause(0.1)


State_UseSimDataButton = app.UseSimDataButton.Value;
app.UseSimDataButton.Value = 0;


% Now, we define our new entry to the SimSetupData
SimSetupData = struct;
SimSetupData.Name = app.MATDSSExcelSheet;

MATDSSApp_Details(app,'Copying Tables to SimData')
pause(0.1)
% Copying all tables
SimSetupData.Table = struct;
SimSetupData.Table.DER = app.DERsTable.Data;
SimSetupData.Table.VTrackingListBoxValues = app.VTrackingListBox.Value;
SimSetupData.Table.ITrackingListBoxValues = app.ITrackingListBox.Value;
SimSetupData.Table.VTrackingListBoxItems = app.VTrackingListBox.Items;
SimSetupData.Table.ITrackingListBoxItems = app.ITrackingListBox.Items;
SimSetupData.Table.CASettings = app.ControlAreasTable.Data;  
SimSetupData.Table.ControllersSettings = app.ControllersSettingsTable.Data;
SimSetupData.Table.VDERSettings = app.VDERsTable.Data;
SimSetupData.Table.Disturbance = app.DisturbanceTable.Data;



MATDSSApp_Details(app,'Copying Tables-Data to SimData')
pause(0.1)
% Copying all tables data
SimSetupData.TableData = struct;
SimSetupData.TableData.DER = app.DERsTable.Data.Variables;
SimSetupData.TableData.CASettings = app.ControlAreasTable.Data.Variables;  
SimSetupData.TableData.ControllersSettings = app.ControllersSettingsTable.Data.Variables;
SimSetupData.TableData.VDERSettings = app.VDERsTable.Data.Variables;
SimSetupData.TableData.Disturbance = app.DisturbanceTable.Data.Variables;
SimSetupData.TableData.ControllersSettings = cellfun(@lower,SimSetupData.TableData.ControllersSettings,'UniformOutput',false);
SimSetupData.TableData.CASettings = cellfun(@lower,SimSetupData.TableData.CASettings,'UniformOutput',false);


% Generating IV Tables
MATDSS = struct;
MATDSS.Table.DER = SimSetupData.Table.DER;
MATDSS.TableData.DER = SimSetupData.TableData.DER;
SimSetupData.Table.VI = MATDSSApp_VITableFun(app,MATDSS);
SimSetupData.TableData.VI = SimSetupData.Table.VI.Variables;
clear MATDSS;

MATDSSApp_Details(app,'Copying Control Parameters to SimData')
pause(0.1)
% Copying all control parameters data
SimSetupData.ControlParameters = struct;
SimSetupData.ControlParameters.GeneralGain = app.GeneralGainEditField.Value;
SimSetupData.ControlParameters.RoU = app.RoUEditField.Value;
SimSetupData.ControlParameters.LoopSize = app.LoopsizeEditField.Value;
SimSetupData.ControlParameters.DisturbanceTable = app.DisturbanceTable.Data;

MATDSSApp_Details(app,'Copying Time Parameters to SimData')
pause(0.1)
% Copying Time Settings
SimSetupData.TimeSettings = struct;
SimSetupData.TimeSettings.Duration = app.DurationEditField.Value;
SimSetupData.TimeSettings.SimTimeStep = app.SimTimeStepEditField.Value;
SimSetupData.TimeSettings.ContTimeStep = app.ContTimeStepEditField.Value;
SimSetupData.TimeSettings.MeasTimeStep = app.MeasTimeStepEditField.Value;
SimSetupData.TimeSettings.MeasDelay = app.MeasDelayEditField.Value;
SimSetupData.TimeSettings.StabilizationTime = app.StabilizationTimeEditField.Value;

MATDSSApp_Details(app,'Copying Simulation Settings to SimData')
pause(0.1)
% Copying Simulation Settings
SimSetupData.P0SetFunctionSettings = struct;
SimSetupData.P0SetFunctionSettings.Type = app.P0SetFunctionButtonGroup.SelectedObject.Text;
SimSetupData.P0SetFunctionSettings.kWVal = app.ChangekWEditField.Value;
SimSetupData.P0SetFunctionSettings.TrackP0Toggle = app.TrackP0Button.Value;


% Run Initiallization and Generate all needed variables/matrices for the
% simulation
GlobalGain = str2double(app.GeneralGainEditField.Value);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               System Initialization               %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSSApp_Details(app,'Generating Simulation Matrices...')
pause(0.1)

MATDSSApp_Details(app,'Initializing Simulation...')
pause(0.1)
[MATDSS] = MATDSSApp_Initialization(app,GlobalGain);

MATDSS.Table = SimSetupData.Table;
MATDSS.TableData = SimSetupData.TableData;
if MATDSS.Sim.DSSStartOk % if DSS is loaded and linked successfully


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               Creating DER Devices                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSSApp_Details(app,'Generating DERs')
    pause(0.1)
    [MATDSS, DER] = MATDSS_DER(MATDSS);

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            Creating Disturbance Loads             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSSApp_Details(app,'Generating Disturbances')
    pause(0.1)
    MATDSS = MATDSS_Disturbance(MATDSS);
    for i = 1:20
        MATDSS.Sim.DSSSolution.Solve;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                    YBus Matrix                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSSApp_Details(app,'Generating Ybus matrices')
    pause(0.1)
    MATDSS = MATDSS_YBusMatrix(MATDSS,1);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            Time Index Initiallization             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = 1; %time = 0 s, we do all initiallization before the for loop.
    % MATDSSApp_ProgressBar(0,app.ProgressBar1,'t = 0 s', app.ProgressBarStatus);
    MATDSS.Sim.t = t;
    at = MATDSS.Time.Sim.TimeSpan(t); %actual time in seconds
    MATDSS.Sim.at = at;
    % app.MATDSSRunVariables.t = t;
    % app.MATDSSRunVariables.at = at;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%           Initializing the Measurements           %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSSApp_Details(app,'Generating Measurement vectors')
    pause(0.1)
    MATDSS = MATDSS_Measurement(MATDSS);
    MATDSS.ControlSignals.P0Set = ones(length(MATDSS.Sim.Meas.at),1).*sum(MATDSS.Meas.P0(:,t));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      Defining Controllers and Control Areas       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSSApp_Details(app,'Defining controllers (this might take a while...)')
    pause(0.1)
    [MATDSS, DER] = MATDSS_Controllers(MATDSS,DER); %Initialize the controllers


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            Initializing the Controlles            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MATDSS.Sim.RoU = 0;
    [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS,DER);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 DER Initialization                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSSApp_Details(app,'Initiallizing DERs')
    pause(0.1)
    for i = 1:MATDSS.Sim.nDER
        [MATDSS, DER] = MATDSS_DERUpdate(MATDSS,DER,i);
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%          Initializing Disturbance Loads           %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSS = MATDSS_DisturbanceUpdate(MATDSS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%             Live Measurements Updates             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATDSSApp_LiveMeasurements(MATDSS,DER);

end

MATDSSApp_Details(app,'Packing and Saving Simulation Data file')
MATDSS.Sim.DSSObj = [];
MATDSS.Sim.DSSText = [];
MATDSS.Sim.DSSCircuit = [];
MATDSS.Sim.DSSSolution = [];
SimSetupData.MATDSS = MATDSS;
SimSetupData.DER = DER;

currenttime = datetime;
currenttime.Format = 'MMddyyHHmm';
currenttime = char(currenttime);

SimSetupData.currenttime = currenttime;

if exist(fullfile([pwd '\Config\SimSetupMats']),'dir') == 0
    MATDSSApp_Details(app,'SimData folder is not found, creating it now.');
    mkdir([pwd '\Config'], 'SimSetupMats');
    MATDSSApp_Details(app, 'Folder created, saving SimData file now');
end
save([pwd '\Config\SimSetupMats\' app.MATDSSExcelSheet '_' currenttime '.mat'], "SimSetupData", "-v7.3");
    

app.UseSimDataButton.Value = State_UseSimDataButton;
MATDSSApp_Status(app,'Simulation Data Exported Successfully', ['Exporting completed to '  app.MATDSSExcelSheet '_' currenttime '.mat']);




end

%}