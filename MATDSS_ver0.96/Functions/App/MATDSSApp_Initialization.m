function MATDSS = MATDSSApp_Initialization(app, Gain)
% MATDSSApp_Initialization(app, Gain)
% This function initializes the MATDSS structure and sets up the 
% simulation environment based on the configurations provided in the 
% MATDSS Application GUI. It prepares the necessary data structures 
% (MATDSS, DER) that are used for simulations and removes any direct 
% dependencies from MATDSS_Sim, allowing it to run independently 
% (e.g., from Simulink).
%
% This initialization sets up time parameters, controller settings, and 
% measurement structures. It also handles various simulation-specific 
% configurations such as the time step, stabilization time, and the 
% input/output profiles needed during the simulation.
%
% The Gain parameter is optional and may be used for additional settings 
% not explicitly defined in the function.
%
% Parameters:
%   - app: The MATDSS application instance, containing GUI elements and 
%          configurations.
%   - Gain: A parameter related to controller gains (optional).
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com


% Check if output is required, otherwise update status and pause for UI feedback
if nargout == 0
    MATDSSApp_Status(app, 'Initializing...', 'Initializing the simulation...');
    pause(0.1); % Small pause for user to observe status update in GUI
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                Main MATDSS Struct                 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the main MATDSS struct that will handle all variables needed for 
% simulation. This struct will contain time, measurements, simulation 
% parameters, plot configurations, and controller information.
MATDSS = struct;

% Deprecated L2C struct replaced by Cont struct for controller details.
% MATDSS.L2C = struct; 

MATDSS.Time = struct;   % Contains time-related data such as simulation steps and duration
MATDSS.Meas = struct;   % Stores measurement data used by the controller
MATDSS.Sim = struct;    % Stores simulation-specific parameters (e.g., loads, changes)
MATDSS.Plot = struct;   % Handles variables and settings related to plots
MATDSS.Cont = struct;   % Stores controller information (new struct replacing L2C)
MATDSS.Disturbance = struct; % Handles load disturbances during the simulation
MATDSS.UseSimDataFlag = app.UseSimDataButton.Value; % Flag to check if 'Use SimData' is enabled

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  Time Parameters                  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Retrieve time-related parameters from the GUI fields for the simulation
TDuration = str2double(app.DurationEditField.Value); % Total simulation duration (seconds)
StabilizationTime = str2double(app.StabilizationTimeEditField.Value); 
% Extra time given for stabilization before and after the main event (e.g., load change)

SimTimeStep = str2double(app.SimTimeStepEditField.Value) * 1e-3; 
% Time step for the simulation, defined in milliseconds and converted to seconds

SimTimeEnd = TDuration + 1 * StabilizationTime; 
% The total simulation end time, including stabilization period

SimTimeSpan = 0:SimTimeStep:SimTimeEnd; 
% The full time span for the simulation, from 0 to the end time, using the defined time step

% Ensure these time steps align with the main simulation time step (SimTimeStep)
ContTimeStep = str2double(app.ContTimeStepEditField.Value) * 1e-3; 
% Controller update time step (in seconds)

MeasTimeStep = str2double(app.MeasTimeStepEditField.Value) * 1e-3; 
% Measurement update time step (in seconds)

MeasDelay = str2double(app.MeasDelayEditField.Value) * 1e-3; 
% Delay in the measurements (in seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 MATDSS.Time Struct                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Populate the MATDSS.Time struct with time-related data
MATDSS.Time.Sim.TimeSpan = SimTimeSpan;  % Simulation time span
MATDSS.Time.Sim.TimeStep = SimTimeStep;  % Main simulation time step
MATDSS.Time.Sim.ST = StabilizationTime;  % Stabilization time before and after events
MATDSS.Time.Cont.TimeStep = ContTimeStep; % Controller time step
MATDSS.Time.Meas.TimeStep = MeasTimeStep; % Measurement update time step
MATDSS.Time.Meas.Delay = MeasDelay;       % Measurement delay
% MATDSS.Time.Meas.QueueIndex = []; % Optionally, store indices for measurement queues (if needed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  MATDSS.Sim Struct                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define simulation-specific parameters in the MATDSS.Sim struct
DSSFileName = [app.OpenDSSFilesListBox.Tag app.OpenDSSFilesListBox.Value]; 
% Retrieve OpenDSS file name from the application GUI

MATDSS.Sim.Bus0 = 1:3; 
% The first three nodes (Bus 0) are reserved for the interface bus in the controller algorithm

MATDSS.Sim.DSSFileName = DSSFileName; 
% Store the DSS file name for use in the simulation

MATDSS.Sim.nDER = 0; 
% Initially, there are 0 DERs (Distributed Energy Resources) defined

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                   P0Set Function                  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the P0SetFunction struct to handle power setpoints (e.g., for DER)
MATDSS.Sim.P0SetFunction = struct;

% Retrieve the type of P0SetFunction from the GUI (UnitStep or Ramp)
P0SetFunctionType = app.P0SetFunctionButtonGroup.SelectedObject.Text; 
MATDSS.Sim.P0SetFunction.FunctionType = P0SetFunctionType; 

MATDSS.Sim.P0SetFunction.Flag = 1; 
% Flag to activate the P0 set function

MATDSS.Sim.P0SetFunction.ChangekW = str2double(app.ChangekWEditField.Value) * 1e3; 
% Convert kW change value from GUI input into Watts

% Deprecated hardcoded function types (left as comments for future reference)
% MATDSS.Sim.P0SetFunction.FunctionType = 'UnitStep'; % UnitStep or Ramp
% MATDSS.Sim.P0SetFunction.FunctionType = 'Ramp'; % UnitStep or Ramp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 MATDSS.Meas Struct                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize measurement variables in the MATDSS.Meas struct
MATDSS.Meas.at = [];            % Placeholder for time instances at which measurements are taken
MATDSS.Meas.VMagProfile = [];   % Voltage magnitude profile over time
MATDSS.Meas.VMagProfilePu = []; % Voltage magnitude profile in per-unit (p.u.)
MATDSS.Meas.VProfile = [];      % Full voltage profile
MATDSS.Meas.VProfileMag = [];   % Voltage magnitude data
MATDSS.Meas.VProfileAng = [];   % Voltage angle data
MATDSS.Meas.IinjProfile = [];   % Current injection profile
MATDSS.Meas.SProfile = [];      % Complex power (S) profile
MATDSS.Meas.PProfile = [];      % Active power (P) profile
MATDSS.Meas.QProfile = [];      % Reactive power (Q) profile
MATDSS.Meas.P0 = [];            % Initial active power setpoint (P0)
MATDSS.Meas.Q0 = [];            % Initial reactive power setpoint (Q0)
MATDSS.Meas.IProfile = [];      % Current profile over time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%               MATDSS.Sim.Meas Struct              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize empty arrays for simulation measurement variables
MATDSS.Sim.Meas.at = [];              % Time vector for measurements
MATDSS.Sim.Meas.VMagProfile = [];     % Voltage magnitude profile
MATDSS.Sim.Meas.VMagProfilePu = [];   % Voltage magnitude profile in per unit
MATDSS.Sim.Meas.VProfile = [];        % Voltage profile (complex numbers)
MATDSS.Sim.Meas.VProfileMag = [];     % Voltage profile magnitudes
MATDSS.Sim.Meas.VProfileAng = [];     % Voltage profile angles
MATDSS.Sim.Meas.IinjProfile = [];     % Injected current profile
MATDSS.Sim.Meas.SProfile = [];        % Apparent power profile
MATDSS.Sim.Meas.PProfile = [];        % Active power profile
MATDSS.Sim.Meas.QProfile = [];        % Reactive power profile
MATDSS.Sim.Meas.P0 = [];              % Initial active power
MATDSS.Sim.Meas.Q0 = [];              % Initial reactive power
MATDSS.Sim.Meas.IProfile = [];        % Current profile

% Define measurement time vectors
MATDSS.Meas.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Meas.TimeStep)==0)); % Measurement times based on time step
MATDSS.Sim.Meas.at = MATDSS.Time.Sim.TimeSpan; % Copy the simulation time span

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%               MATDSS.Disturbance Struct            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize disturbance-related variables
MATDSS.Disturbance.DSSName = [];      % Name of the disturbance
MATDSS.Disturbance.ID = [];           % ID of the disturbance in the table
MATDSS.Disturbance.Bus = [];          % Bus associated with the disturbance
MATDSS.Disturbance.Phases = [];       % Phases involved in the disturbance
MATDSS.Disturbance.Nphase = [];       % Number of phases involved
MATDSS.Disturbance.ConnType = [];     % Connection type (wye or delta)
MATDSS.Disturbance.P = [];            % Power in kW
MATDSS.Disturbance.pf = [];           % Power factor
MATDSS.Disturbance.t = [];            % Time of disturbance [start_time, end_time]
MATDSS.Disturbance.Enabled = [];      % Enable flag for the disturbance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 MATDSS.Cont Gain                  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Controller gain and settings
MATDSS.Cont.Gain = Gain;  % Set the main controller gain
MATDSS.Cont.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Cont.TimeStep)==0)); % Controller time based on step size
MATDSS.Cont.RoU = str2double(app.RoUEditField.Value);  % RoU value from app input
MATDSS.Cont.LoopSize = str2double(app.LoopsizeEditField.Value); % Loop size value from app input
MATDSS.Cont.DisableControllersFlag = app.DisableControllersButton.Value; % Disable controller flag
MATDSS.Cont.s = app.TrackP0Button.Value; % Signal for tracking P0 (active power)

% LowPass Filter Settings
MATDSS.Cont.LPF = struct;  % Initialize LowPass Filter structure
MATDSS.Cont.LPF.EnableLPFFlag = app.EnableLPFCheckBox.Value; % Enable/disable LPF based on app checkbox
MATDSS.Cont.LPF.Tc = 0; % Initialize filter cutoff time constant

% If LowPass Filter is enabled
if MATDSS.Cont.LPF.EnableLPFFlag
    % Read the filter cutoff time constant from the app
    MATDSS.Cont.LPF.Tc = app.TccutoffperiodinsEditField.Value;
    MATDSS.Cont.LPF.Ts = MATDSS.Time.Cont.TimeStep;  % Set the filter sampling time
    MATDSS.Cont.LPF.fc = 1/MATDSS.Cont.LPF.Tc;      % Calculate cutoff frequency
    MATDSS.Cont.LPF.wc = 2*pi*MATDSS.Cont.LPF.fc;   % Calculate cutoff angular frequency
    MATDSS.Cont.LPF.TF = tf([1],[1/MATDSS.Cont.LPF.wc, 1]); % Create transfer function
    MATDSS.Cont.LPF.DTF = c2d(MATDSS.Cont.LPF.TF,MATDSS.Cont.LPF.Ts,'tustin'); % Discretize transfer function

    % Discretized transfer function setup (used for filtering signals)
    MATDSS.Cont.LPF.Num = MATDSS.Cont.LPF.DTF.Numerator{:};  % Get numerator coefficients
    MATDSS.Cont.LPF.Den = MATDSS.Cont.LPF.DTF.Denominator{:}; % Get denominator coefficients

    % Calculate filter coefficients for real-time filtering
    MATDSS.Cont.LPF.d1d0 = MATDSS.Cont.LPF.Den(2)/MATDSS.Cont.LPF.Den(1);  % Denominator ratio
    MATDSS.Cont.LPF.n0d0 = MATDSS.Cont.LPF.Num(1)/MATDSS.Cont.LPF.Den(1);  % Numerator ratio
    MATDSS.Cont.LPF.n1d0 = MATDSS.Cont.LPF.Num(2)/MATDSS.Cont.LPF.Den(1);  % Numerator ratio
end

% Read PID Controller Settings
MATDSS.Cont.PID = struct;  % Initialize PID controller structure
MATDSS.Cont.PID.EnablePFlag = app.EnablePControlCheckBox.Value;  % Enable P control flag
MATDSS.Cont.PID.EnableDFlag = app.EnableDControlCheckBox.Value;  % Enable D control flag
MATDSS.Cont.PID.DControlAppliedto = 'none';  % Default: no D control applied
MATDSS.Cont.PID.PControlAppliedto = 'none';  % Default: no P control applied

% If P control is enabled
if MATDSS.Cont.PID.EnablePFlag
    MATDSS.Cont.PID.Kp = app.KpEditField.Value;  % Proportional gain
    % P control kappa values (user-defined parameters)
    MATDSS.Cont.PID.kappa_p.lambda = app.kappa_plambdaEditField.Value;
    MATDSS.Cont.PID.kappa_p.mu = app.kappa_pmuEditField.Value;
    MATDSS.Cont.PID.kappa_p.eta = app.kappa_petaEditField.Value;
    MATDSS.Cont.PID.kappa_p.psi = app.kappa_ppsiEditField.Value;
    MATDSS.Cont.PID.kappa_p.gamma = app.kappa_pgammaEditField.Value;
    MATDSS.Cont.PID.kappa_p.nu = app.kappa_pnuEditField.Value;
    MATDSS.Cont.PID.kappa_p.zeta = app.kappa_pzetaEditField.Value;
    MATDSS.Cont.PID.PControlAppliedto = app.PControlAppliedtoButtonGroup.SelectedObject.Text; % Set control application type
end

% If D control is enabled
if MATDSS.Cont.PID.EnableDFlag
    MATDSS.Cont.PID.Kd = app.KdEditField.Value;  % Derivative gain
    % D control kappa values (user-defined parameters)
    MATDSS.Cont.PID.kappa_d.lambda = app.kappa_dlambdaEditField.Value;
    MATDSS.Cont.PID.kappa_d.mu = app.kappa_dmuEditField.Value;
    MATDSS.Cont.PID.kappa_d.eta = app.kappa_detaEditField.Value;
    MATDSS.Cont.PID.kappa_d.psi = app.kappa_dpsiEditField.Value;
    MATDSS.Cont.PID.kappa_d.gamma = app.kappa_dgammaEditField.Value;
    MATDSS.Cont.PID.kappa_d.nu = app.kappa_dnuEditField.Value;
    MATDSS.Cont.PID.kappa_d.zeta = app.kappa_dzetaEditField.Value;

    % Choose the D signal type (Error or measurement-based)
    MATDSS.Cont.PID.Dsignal = app.DSignalButtonGroup.SelectedObject.Text;
    MATDSS.Cont.PID.DControlAppliedto = app.DControlAppliedtoButtonGroup.SelectedObject.Text; % Set where D control is applied
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                     DER Struct                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DER Struct
% DER = struct; % details of DER struct are given in MATDSS_DERNew function.
% The Distributed Energy Resources (DER) struct will store the DER 
% information. This is initialized in a separate function called MATDSS_DERNew.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                Use Loaded Sim Data                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MATDSS.UseSimDataFlag
    % If the flag for using simulation data is set to true, the following
    % block will execute. It indicates that pre-loaded simulation data
    % from the MATDSS application will be used during initialization.

    % Notify the user that DER information will not be imported from
    % SimData, as DER definitions depend on OpenDSS environment links.
    MATDSSApp_Details(app,'Checking saved Data in "SimData"');
    
    % Warning message for the user: The SimData will only be used to speed up
    % the generation of large matrices such as Ybus and ABMH. The actual
    % DER and controller data used in the simulation will come from the
    % MATDSS application window (unless explicitly reloaded from SimData).
    MATDSSApp_Details(app,'Warning - Simulation Data loaded will only be used for Ybus, ABMH and large matrices generation to speed up initiallization time. This means you should make sure that all changes made for the tables that might affect the shape/values of these matrices are already done before you generate SimData. The properties of DERs (excluding bus, type, or # of phases), Controllers (excluding monitored phases, defined control areas), Disturbances, can be modifed after loading SimData. The valeus lsited in the tables within MATDSSApp will be used in the simulation. To restore SimData values, click on "Reload SimData" or "Load Config. From Saved SimData" to reload them to the tables.');
    
    % Assigning loaded simulation data to the MATDSS struct for future use.
    MATDSS.LoadedSimData = app.LoadedSimData;
    
    % A brief pause is added here to allow the system to process the loading.
    pause(0.1)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      Reading the Configurations and Settings      %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section handles reading and setting the configurations for the simulation
% environment based on the MATDSS app's tables.

% Notify the user that the app is reading configurations and settings.
MATDSSApp_Details(app,'Reading Configurations and Settings');
pause(0.1)  % Small delay to process the message.

% Read DER Table from the app's UI (MATDSS_DER.xlsx or equivalent source).
MATDSS.Table.DER = app.DERsTable.Data; % Extract the data from the DERs table.
MATDSS.TableData.DER = MATDSS.Table.DER.Variables; % Store the DER table's variable data.

% Generate the Voltage and Current (VI) table, a key component for 
% tracking voltage and current levels in the simulation environment.
VITable = MATDSSApp_VITableFun(app,MATDSS); 
MATDSS.Table.VI = VITable; % Store the VI table in MATDSS structure.
MATDSS.TableData.VI = VITable.Variables; % Extract the variables from VI table.

% % % Legacy code for reading VI table from an Excel file. Commented out for now.
% % % If the MATDSS_VI.xlsx file exists and contains data for the selected sheet,
% % % it will be loaded and used for voltage phase tracking.
% % if exist([pwd '\Config\MATDSS_VI.xlsx'], 'file')
% %     MATDSSVISheetnames =  cellstr(sheetnames('MATDSS_VI.xlsx')); % Get list of sheet names.
% %     MATDSSVIFlag = MATDSS_StrComp(MATDSSVISheetnames,app.OpenDSSFilesListBox.Value); % Verify target sheet's presence.
% %     
% %     if MATDSSVIFlag > 0
% %         VIReadOptions = detectImportOptions('MATDSS_VI.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % Set import options.
% %         VIReadOptions = setvartype(VIReadOptions,VIReadOptions.VariableNames(:),'string'); % Set variable types to string.
% %         MATDSS.Table.VI = readtable([pwd '\Config\MATDSS_VI.xlsx'],VIReadOptions); % Load the data into the table.
% %         MATDSS.TableData.VI = MATDSS.Table.VI.Variables; % Extract the variables from the table.
% %         clear VIReadOptions  % Clear the options to save memory.
% %     end
% % end

% Read the Control Areas (CA) table for managing the system's control areas.
MATDSS.Table.CASettings = app.ControlAreasTable.Data; % Read the Control Areas settings from the UI table.
MATDSS.TableData.CASettings = MATDSS.Table.CASettings.Variables; % Extract the variables.
MATDSS.TableData.CASettings = cellfun(@lower,MATDSS.TableData.CASettings,'UniformOutput',false); % Convert the data to lowercase for consistency.

% % % Legacy code for reading Control Areas from an Excel file. Commented out for now.
% % % If the MATDSS_ControlAreas.xlsx file exists, the system will attempt to load 
% % % the settings for control areas from the specified sheet.
% % if exist([pwd '\Config\MATDSS_ControlAreas.xlsx'], 'file')
% %     MATDSSCASettingsSheetnames =  cellstr(sheetnames('MATDSS_ControlAreas.xlsx')); % Get list of sheet names.
% %     MATDSSCASettingsFlag = MATDSS_StrComp(MATDSSCASettingsSheetnames,app.OpenDSSFilesListBox.Value); % Check if target sheet exists.
% %     if MATDSSCASettingsFlag > 0
% %         CASettingsReadOptions = detectImportOptions('MATDSS_ControlAreas.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % Set read options.
% %         CASettingsReadOptions = setvartype(CASettingsReadOptions,CASettingsReadOptions.VariableNames(:),'string'); % Set variable types to string.
% %         MATDSS.Table.CASettings = readtable([pwd '\Config\MATDSS_ControlAreas.xlsx'],CASettingsReadOptions); % Load data into table.
% %         MATDSS.TableData.CASettings= MATDSS.Table.CASettings.Variables; % Extract the variables.
% %         MATDSS.TableData.CASettings = cellfun(@lower,MATDSS.TableData.CASettings,'UniformOutput',false); % Convert to lowercase.
% %     end
% % end

% Read the Controllers Settings from the UI table in MATDSS.
MATDSS.Table.ControllersSettings = app.ControllersSettingsTable.Data; % Read the Controllers settings from the UI table.
MATDSS.TableData.ControllersSettings= MATDSS.Table.ControllersSettings.Variables; % Extract the variable data.
MATDSS.TableData.ControllersSettings = cellfun(@lower,MATDSS.TableData.ControllersSettings,'UniformOutput',false); % Convert to lowercase for consistency.

% % % Legacy code for reading Controller Settings from an Excel file. Commented out for now.
% % % If the MATDSS_ControllersSettings.xlsx file exists, the system will try
% % % to load the settings for controllers from the selected sheet.
% % if exist([pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'file')
% %     MATDSSControllersSettingsSheetnames =  cellstr(sheetnames('MATDSS_ControllersSettings.xlsx')); % Get list of sheet names.
% %     MATDSSControllersSettingsFlag = MATDSS_StrComp(MATDSSControllersSettingsSheetnames,app.OpenDSSFilesListBox.Value); % Check if target sheet exists.
% %     if MATDSSControllersSettingsFlag > 0
% %         ControllersSettingsReadOptions = detectImportOptions('MATDSS_ControllersSettings.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % Set import options.
% %         ControllersSettingsReadOptions = setvartype(ControllersSettingsReadOptions,ControllersSettingsReadOptions.VariableNames(:),'string'); % Set variable types to string.
% %         MATDSS.Table.ControllersSettings = readtable([pwd '\Config\MATDSS_ControllersSettings.xlsx'],ControllersSettingsReadOptions); % Load data into table.
% %         MATDSS.TableData.ControllersSettings= MATDSS.Table.ControllersSettings.Variables; % Extract the variables.
% %         MATDSS.TableData.ControllersSettings = cellfun(@lower,MATDSS.TableData.ControllersSettings,'UniformOutput',false); % Convert to lowercase.
% %     end
% % end

% Try reading MATDSS_VDERs.xlsx file
MATDSS.Table.VDERSettings = app.VDERsTable.Data; % Load VDER settings from the table in the MATDSS application
MATDSS.TableData.VDERSettings = MATDSS.Table.VDERSettings.Variables; % Store the VDER settings data as variables
% MATDSS.TableData.VDERSettings = cellfun(@lower,MATDSS.TableData.VDERSettings,'UniformOutput',false); % Optional code to convert to lowercase, currently commented out

% Try reading MATDSS_VDERs.xlsx file
MATDSS.Table.Disturbance = app.DisturbanceTable.Data; % Load disturbance data from the table in the MATDSS application
MATDSS.TableData.Disturbance = MATDSS.Table.Disturbance.Variables; % Store the disturbance data as variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              Compile OpenDSS Circuit              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 0
    MATDSSApp_Details(app,'Initiating OpenDSS Engine'); % Display initialization status
    pause(0.1) % Pause briefly to ensure smooth initialization
end

% The compile command "text" that will be used to run/execute the DSS file later
MATDSS.Sim.DSSCompile = ['Compile "' DSSFileName '"']; % Store the DSS compile command with the file name

% Initiate DSS link
MATDSS.Sim.DSSObj = actxserver('OpenDSSEngine.DSS'); % Create an ActiveX server for the OpenDSS engine

% Flag if link is established
MATDSS.Sim.DSSStartOk = MATDSS.Sim.DSSObj.Start(0); % Attempt to start the DSS engine (returns 1 if successful)

% If DSS link was successfully established
if MATDSS.Sim.DSSStartOk
    % Define the Text interface (to run commands from MATLAB)
    MATDSS.Sim.DSSText = MATDSS.Sim.DSSObj.Text; % Initialize the Text interface for command execution

    % Compile the DSS file to load the circuit
    MATDSS.Sim.DSSText.command = MATDSS.Sim.DSSCompile; % Compile the DSS file using the command stored earlier

    % Set up the interface variables DSSCircuit and DSSSolution
    MATDSS.Sim.DSSCircuit = MATDSS.Sim.DSSObj.ActiveCircuit; % Get the active DSS circuit
    MATDSS.Sim.DSSSolution = MATDSS.Sim.DSSCircuit.Solution; % Get the solution object for the active circuit
    % MATDSS.Sim.DSSLines % Optional line, currently commented out
else
    % If DSS engine failed to start, show an error message
    msgbox('DSSEngine Failed to start!'); % Display a failure message
    % Log error details
    MATDSSApp_Details(app,{'!****************************************!';...
        '!* Could not initiate link with OpenDSS *!';...
        '!****************************************!';...
        ' ';'Simulation Aborted!';''});
    return % Abort the process and exit the function
end

% This flag enables logging of all V, I, P, and Q for all buses
% Data is stored in MATDSS.Sim.Meas for plotting purposes (not accessible to controllers)
MATDSS.Sim.Meas.Flag = 1; % Enable the measurement logging flag

% If there are no output arguments, update the MATDSS application instance
if nargout == 0
   app.Main.MATDSS = MATDSS; % Save the MATDSS data back into the app's main structure
end

% Indicate that the initialization is complete
MATDSSApp_Details(app,'Initialization Completed!'); % Display a message indicating successful initialization
pause(0.1) % Pause briefly

end % End of the MATDSSApp_Initialization function




%% Old Function Code
%{
function MATDSS = MATDSSApp_Initialization(app,Gain)
% this function will initiallize the MATDSS and DER variables according to
% the configurations set in MATDSS Application. We will move out any
% dependce of MATDSS_Sim to the application, which will allow us to run
% that function later from Simulink for example!
if nargout == 0
    MATDSSApp_Status(app,['Initializing...'],['Initializing the simulation...']);
    pause(0.1)
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                Main MATDSS Struct                 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main MATDSS Struct that will handle all variables in our run.
MATDSS = struct;
% MATDSS.L2C = struct; This is the old L2C struct that is deprecated
MATDSS.Time = struct; % Time struct that will contain time information for simulation, measurements and L2C
MATDSS.Meas = struct; % Measurement Struct that will contain the measurements for our controller to use
MATDSS.Sim = struct; % All modeling related parameters that may not be available for
% our controller (like loads information and load changes and other
% details that we might not have access to in the controller)
MATDSS.Plot = struct; % This is the struct that will handle the plot related variables and features.
MATDSS.Cont = struct; % Controllers struct that will contain all information about the controllers that we have in the system (This is the new struct replacing L2C in previous versions)
MATDSS.Disturbance = struct; % Struct that will handle load disturbances
MATDSS.UseSimDataFlag = app.UseSimDataButton.Value; % Check if Use SimData is enabled




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  Time Parameters                  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Main Simulation Time Parameters (read the parameters from the MATDSS
% Application GUI
% Start time is always t=0
TDuration = str2double(app.DurationEditField.Value);
StabilizationTime = str2double(app.StabilizationTimeEditField.Value); % Extra time that will be given for the system before and after the actual P0SetFunction is active. This is used to help stabilize and remove any fluctuations at the beginning.
SimTimeStep = str2double(app.SimTimeStepEditField.Value)*1e-3; % Time step of the main simulation, always recommended to be small steps
SimTimeEnd = TDuration + 1*StabilizationTime; % End time of simulation
SimTimeSpan = 0:SimTimeStep:SimTimeEnd; % simulation time window in seconds

% For the following time steps, make sure that these are reachable (mod =
% 0) uing the SimTimeStep
ContTimeStep = str2double(app.ContTimeStepEditField.Value)*1e-3;%*SimTimeStep;%0.1; % Level 2 Controller time step in seconds
MeasTimeStep = str2double(app.MeasTimeStepEditField.Value)*1e-3;%SimTimeStep;%0.25; % Measurement updates time step in seconds
MeasDelay = str2double(app.MeasDelayEditField.Value)*1e-3; % delay in the measurements in seconds


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 MATDSS.Time Struct                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Time Struct
MATDSS.Time.Sim.TimeSpan = SimTimeSpan; % Simulation Main time span
MATDSS.Time.Sim.TimeStep = SimTimeStep; % Simuation Main time step
MATDSS.Time.Sim.ST = StabilizationTime; % Simulation stabilization time
MATDSS.Time.Cont.TimeStep = ContTimeStep; % Level 2 Controller time step
MATDSS.Time.Meas.TimeStep = MeasTimeStep;
MATDSS.Time.Meas.Delay = MeasDelay;
% MATDSS.Time.Meas.QueueIndex = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  MATDSS.Sim Struct                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sim Struct
DSSFileName = [app.OpenDSSFilesListBox.Tag app.OpenDSSFilesListBox.Value];
MATDSS.Sim.Bus0 = 1:3; %the first three nodes are for the interface bus (bus 0 in our controller algorithm)
MATDSS.Sim.DSSFileName = DSSFileName;
MATDSS.Sim.nDER = 0; % initially we have 0 DERs defined!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                   P0Set Function                  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P0Set Function
MATDSS.Sim.P0SetFunction = struct;
P0SetFunctionType = app.P0SetFunctionButtonGroup.SelectedObject.Text;
MATDSS.Sim.P0SetFunction.FunctionType = P0SetFunctionType; % UnitStep or Ramp
MATDSS.Sim.P0SetFunction.Flag = 1;
MATDSS.Sim.P0SetFunction.ChangekW = str2double(app.ChangekWEditField.Value)*1e3;
% MATDSS.Sim.P0SetFunction.FunctionType = 'UnitStep'; % UnitStep or Ramp
% MATDSS.Sim.P0SetFunction.FunctionType = 'Ramp'; % UnitStep or Ramp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 MATDSS.Meas Struct                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSS.Meas.at = [];
MATDSS.Meas.VMagProfile = [];
MATDSS.Meas.VMagProfilePu = [];
MATDSS.Meas.VProfile = [];
MATDSS.Meas.VProfileMag = [];
MATDSS.Meas.VProfileAng = [];
MATDSS.Meas.IinjProfile = [];
MATDSS.Meas.SProfile = [];
MATDSS.Meas.PProfile = [];
MATDSS.Meas.QProfile = [];
MATDSS.Meas.P0 = [];
MATDSS.Meas.Q0 = [];
MATDSS.Meas.IProfile = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               MATDSS.Sim.Meas Struct              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSS.Sim.Meas.at = [];
MATDSS.Sim.Meas.VMagProfile = [];
MATDSS.Sim.Meas.VMagProfilePu = [];
MATDSS.Sim.Meas.VProfile = [];
MATDSS.Sim.Meas.VProfileMag = [];
MATDSS.Sim.Meas.VProfileAng = [];
MATDSS.Sim.Meas.IinjProfile = [];
MATDSS.Sim.Meas.SProfile = [];
MATDSS.Sim.Meas.PProfile = [];
MATDSS.Sim.Meas.QProfile = [];
MATDSS.Sim.Meas.P0 = [];
MATDSS.Sim.Meas.Q0 = [];
MATDSS.Sim.Meas.IProfile = [];


MATDSS.Meas.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Meas.TimeStep)==0));
MATDSS.Sim.Meas.at = MATDSS.Time.Sim.TimeSpan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               MATDSS.Sim.Meas Struct              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSS.Disturbance.DSSName = [];
MATDSS.Disturbance.ID = []; % ID Number in the table
MATDSS.Disturbance.Bus = []; % Busname
MATDSS.Disturbance.Phases = []; % Phases
MATDSS.Disturbance.Nphase = []; % number of phases
MATDSS.Disturbance.ConnType = []; % Connection type, wye or delta
MATDSS.Disturbance.P = []; % in kW
MATDSS.Disturbance.pf = []; % power factor
MATDSS.Disturbance.t = []; % in s [start_time, end_time (could be inf.)]
MATDSS.Disturbance.Enabled = []; % this to set if this disturbance is to be enabled in this run or not.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 MATDSS.Cont Gain                  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Controller Gain
MATDSS.Cont.Gain = Gain;
MATDSS.Cont.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Cont.TimeStep)==0));
MATDSS.Cont.RoU = str2double(app.RoUEditField.Value);
MATDSS.Cont.LoopSize = str2double(app.LoopsizeEditField.Value);
MATDSS.Cont.DisableControllersFlag = app.DisableControllersButton.Value;
MATDSS.Cont.s = app.TrackP0Button.Value;



% Read LowPass filter Settings
MATDSS.Cont.LPF = struct;
MATDSS.Cont.LPF.EnableLPFFlag = app.EnableLPFCheckBox.Value;
MATDSS.Cont.LPF.Tc = 0;
% Read LowPass Filter Cutoff time


% LowPass Filter Setup
if MATDSS.Cont.LPF.EnableLPFFlag
    MATDSS.Cont.LPF.Tc = app.TccutoffperiodinsEditField.Value; % Read as a number!
    MATDSS.Cont.LPF.Ts = MATDSS.Time.Cont.TimeStep;
    MATDSS.Cont.LPF.fc = 1/MATDSS.Cont.LPF.Tc;
    MATDSS.Cont.LPF.wc = 2*pi*MATDSS.Cont.LPF.fc;
    MATDSS.Cont.LPF.TF = tf([1],[1/MATDSS.Cont.LPF.wc, 1]);
    MATDSS.Cont.LPF.DTF = c2d(MATDSS.Cont.LPF.TF,MATDSS.Cont.LPF.Ts,'tustin');

    %{
        g_out/ g_in = num/den;
        g_out * den = g_in * num
        
        den = [d0 d1];
        num = [n0 n1];
        
        yd0 + y1d1 = xn0 + x1n1;

        y = y1 -(d1/d0) + x n0/d0 + x1 n1/d0;
        
        d1d0 = d1/d0;
        n0d0 = n0/d0;
        n1d0 = n1/d0;
        
        
    %}

    MATDSS.Cont.LPF.Num = MATDSS.Cont.LPF.DTF.Numerator{:};
    MATDSS.Cont.LPF.Den = MATDSS.Cont.LPF.DTF.Denominator{:};
    
    MATDSS.Cont.LPF.d1d0 = MATDSS.Cont.LPF.Den(2)/MATDSS.Cont.LPF.Den(1);
    MATDSS.Cont.LPF.n0d0 = MATDSS.Cont.LPF.Num(1)/MATDSS.Cont.LPF.Den(1);
    MATDSS.Cont.LPF.n1d0 = MATDSS.Cont.LPF.Num(2)/MATDSS.Cont.LPF.Den(1);
end


% Read PID Controller Settings
MATDSS.Cont.PID = struct;
MATDSS.Cont.PID.EnablePFlag = app.EnablePControlCheckBox.Value;
MATDSS.Cont.PID.EnableDFlag = app.EnableDControlCheckBox.Value;
MATDSS.Cont.PID.DControlAppliedto = 'none'; %default none
MATDSS.Cont.PID.PControlAppliedto = 'none'; %default none

% LowPass Filter Setup
if MATDSS.Cont.PID.EnablePFlag
    MATDSS.Cont.PID.Kp = app.KpEditField.Value;
    MATDSS.Cont.PID.kappa_p.lambda = app.kappa_plambdaEditField.Value;
    MATDSS.Cont.PID.kappa_p.mu = app.kappa_pmuEditField.Value;
    MATDSS.Cont.PID.kappa_p.eta = app.kappa_petaEditField.Value;
    MATDSS.Cont.PID.kappa_p.psi = app.kappa_ppsiEditField.Value;
    MATDSS.Cont.PID.kappa_p.gamma = app.kappa_pgammaEditField.Value;
    MATDSS.Cont.PID.kappa_p.nu = app.kappa_pnuEditField.Value;
    MATDSS.Cont.PID.kappa_p.zeta = app.kappa_pzetaEditField.Value;
    MATDSS.Cont.PID.PControlAppliedto = app.PControlAppliedtoButtonGroup.SelectedObject.Text;
end

if MATDSS.Cont.PID.EnableDFlag
    MATDSS.Cont.PID.Kd = app.KdEditField.Value;
    MATDSS.Cont.PID.kappa_d.lambda = app.kappa_dlambdaEditField.Value;
    MATDSS.Cont.PID.kappa_d.mu = app.kappa_dmuEditField.Value;
    MATDSS.Cont.PID.kappa_d.eta = app.kappa_detaEditField.Value;
    MATDSS.Cont.PID.kappa_d.psi = app.kappa_dpsiEditField.Value;
    MATDSS.Cont.PID.kappa_d.gamma = app.kappa_dgammaEditField.Value;
    MATDSS.Cont.PID.kappa_d.nu = app.kappa_dnuEditField.Value;
    MATDSS.Cont.PID.kappa_d.zeta = app.kappa_dzetaEditField.Value;
    MATDSS.Cont.PID.Dsignal = app.DSignalButtonGroup.SelectedObject.Text; % Error (default) or Y (to use measurements)
    MATDSS.Cont.PID.DControlAppliedto = app.DControlAppliedtoButtonGroup.SelectedObject.Text;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                     DER Struct                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DER Struct
% DER = struct; % details of DER struct are given in MATDSS_DERNew function.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                Use Loaded Sim Data                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MATDSS.UseSimDataFlag
    % Only MATDSS is going to have information from SimData that is loaded.
    % DER variable will be setup in MATDSS_DER function. This cannot be
    % done by loading the data from "SimData" as it has links to actual
    % loads that need to be defined in OpenDSS environment. Therefore, we
    % don't copy or use that information here.

    % We will copy only useful information related to initialization. All
    % gains are going to be used from MATDSS application window. Therefore,
    % if you want to use the data from "SimData", make sure you load them
    % to MATDSS Application windows first by clicking on "Reload SimData".
    MATDSSApp_Details(app,'Checking saved Data in "SimData"');
    MATDSSApp_Details(app,'Warning - Simulation Data loaded will only be used for Ybus, ABMH and large matrices generation to speed up initiallization time. This means you should make sure that all changes made for the tables that might affect the shape/values of these matrices are already done before you generate SimData. The properties of DERs (excluding bus, type, or # of phases), Controllers (excluding monitored phases, defined control areas), Disturbances, can be modifed after loading SimData. The valeus lsited in the tables within MATDSSApp will be used in the simulation. To restore SimData values, click on "Reload SimData" or "Load Config. From Saved SimData" to reload them to the tables.');
    MATDSS.LoadedSimData = app.LoadedSimData;
    pause(0.1)

    

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      Reading the Configurations and Settings      %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading the configurations from within MATDSS App Tables.
MATDSSApp_Details(app,'Reading Configurations and Settings');
pause(0.1)


% Read DER Table
% Reading the MATDSS_DER.xlsx file to define the DERs
MATDSS.Table.DER = app.DERsTable.Data;
MATDSS.TableData.DER = MATDSS.Table.DER.Variables;

VITable = MATDSSApp_VITableFun(app,MATDSS);

MATDSS.Table.VI = VITable;
MATDSS.TableData.VI = VITable.Variables;
%
%
% % Read VI Table
% % Getting Voltage phases indecices for tracking
% if exist([pwd '\Config\MATDSS_VI.xlsx'], 'file')
%     MATDSSVISheetnames =  cellstr(sheetnames('MATDSS_VI.xlsx')); % List of all sheets in the Excel file
%     MATDSSVIFlag = MATDSS_StrComp(MATDSSVISheetnames,app.OpenDSSFilesListBox.Value); % Check if teh target sheet is available in the Excel file
%
%     if MATDSSVIFlag > 0
%         VIReadOptions = detectImportOptions('MATDSS_VI.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % set the import settings
%         VIReadOptions = setvartype(VIReadOptions,VIReadOptions.VariableNames(:),'string');
%         MATDSS.Table.VI = readtable([pwd '\Config\MATDSS_VI.xlsx'],VIReadOptions);%,'Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve');
%         MATDSS.TableData.VI = MATDSS.Table.VI.Variables;
%         clear VIReadOptions
%     end
% end




% Try reading MATDSS_ControlAreas.xlsx file
MATDSS.Table.CASettings = app.ControlAreasTable.Data;
MATDSS.TableData.CASettings= MATDSS.Table.CASettings.Variables;
MATDSS.TableData.CASettings = cellfun(@lower,MATDSS.TableData.CASettings,'UniformOutput',false);

%
% if exist([pwd '\Config\MATDSS_ControlAreas.xlsx'], 'file')
%     MATDSSCASettingsSheetnames =  cellstr(sheetnames('MATDSS_ControlAreas.xlsx')); % List of all sheets in the Excel file
%     MATDSSCASettingsFlag = MATDSS_StrComp(MATDSSCASettingsSheetnames,app.OpenDSSFilesListBox.Value); % Check if teh target sheet is available in the Excel file
%     if MATDSSCASettingsFlag > 0
%         % Read Control Areas Settings
%         CASettingsReadOptions = detectImportOptions('MATDSS_ControlAreas.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % set the import settings
%         CASettingsReadOptions = setvartype(CASettingsReadOptions,CASettingsReadOptions.VariableNames(:),'string');
%         MATDSS.Table.CASettings = readtable([pwd '\Config\MATDSS_ControlAreas.xlsx'],CASettingsReadOptions);
%         MATDSS.TableData.CASettings= MATDSS.Table.CASettings.Variables;
%         MATDSS.TableData.CASettings = cellfun(@lower,MATDSS.TableData.CASettings,'UniformOutput',false);
%     end
% end






% Try reading MATDSS_ControllersSettings.xlsx file
MATDSS.Table.ControllersSettings = app.ControllersSettingsTable.Data;
MATDSS.TableData.ControllersSettings= MATDSS.Table.ControllersSettings.Variables;
MATDSS.TableData.ControllersSettings = cellfun(@lower,MATDSS.TableData.ControllersSettings,'UniformOutput',false);
% if exist([pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'file')
%     MATDSSControllersSettingsSheetnames =  cellstr(sheetnames('MATDSS_ControllersSettings.xlsx')); % List of all sheets in the Excel file
%     MATDSSControllersSettingsFlag = MATDSS_StrComp(MATDSSControllersSettingsSheetnames,app.OpenDSSFilesListBox.Value); % Check if teh target sheet is available in the Excel file
%     if MATDSSControllersSettingsFlag > 0
%         % Read Control Areas Settings
%         ControllersSettingsReadOptions = detectImportOptions('MATDSS_ControllersSettings.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % set the import settings
%         ControllersSettingsReadOptions = setvartype(ControllersSettingsReadOptions,ControllersSettingsReadOptions.VariableNames(:),'string');
%         MATDSS.Table.ControllersSettings = readtable([pwd '\Config\MATDSS_ControllersSettings.xlsx'],ControllersSettingsReadOptions);
%         MATDSS.TableData.ControllersSettings= MATDSS.Table.ControllersSettings.Variables;
%         MATDSS.TableData.ControllersSettings = cellfun(@lower,MATDSS.TableData.ControllersSettings,'UniformOutput',false);
%     end
% end



% Try reading MATDSS_VDERs.xlsx file
MATDSS.Table.VDERSettings = app.VDERsTable.Data;
MATDSS.TableData.VDERSettings = MATDSS.Table.VDERSettings.Variables;
% MATDSS.TableData.VDERSettings = cellfun(@lower,MATDSS.TableData.VDERSettings,'UniformOutput',false);



% Try reading MATDSS_VDERs.xlsx file
MATDSS.Table.Disturbance = app.DisturbanceTable.Data;
MATDSS.TableData.Disturbance = MATDSS.Table.Disturbance.Variables;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              Compile OpenDSS Circuit              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 0
    MATDSSApp_Details(app,'Initiating OpenDSS Engine');
    pause(0.1)
end
% The compile command "text" that will be used to run/excute the DSS file
% later
MATDSS.Sim.DSSCompile = ['Compile "' DSSFileName '"'];

% Initiate DSS link
MATDSS.Sim.DSSObj = actxserver('OpenDSSEngine.DSS');

% Flag if link is established
MATDSS.Sim.DSSStartOk = MATDSS.Sim.DSSObj.Start(0);

if MATDSS.Sim.DSSStartOk
    % Define the Text interface (to run commands from MATLAB)
    MATDSS.Sim.DSSText = MATDSS.Sim.DSSObj.Text;

    % Compile the DSS file first to load the circuit
    MATDSS.Sim.DSSText.command = MATDSS.Sim.DSSCompile;

    % Set up the interface variables DSSCircuit and DSSSolution
    MATDSS.Sim.DSSCircuit = MATDSS.Sim.DSSObj.ActiveCircuit;
    MATDSS.Sim.DSSSolution = MATDSS.Sim.DSSCircuit.Solution;
    %     MATDSS.Sim.DSSLines
else
    msgbox('DSSEngine Failed to start!');
    MATDSSApp_Details(app,{'!****************************************!';...
        '!* Could not initiate link with OpenDSS *!';...
        '!****************************************!';...
        ' ';'Simulation Aborted!';''});
    return
end



% This flag is to enable logging of all V, I, P and Q of all buses and save
% their data in MATDSS.Sim.Meas for plotting purposes. This data is not
% accessible to the controllers
MATDSS.Sim.Meas.Flag = 1;


if nargout == 0
   app.Main.MATDSS = MATDSS;
end

MATDSSApp_Details(app,'Initialization Completed!');
pause(0.1)
end
%}