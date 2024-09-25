function [MATDSS] = MATDSS_Initialize(app, Gain,Flag_UseSimData)
if nargin < 3
    Flag_UseSimData = 0;
end
if Flag_UseSimData
    MATDSS = app.MyRun.MATDSS;
    DER = app.MyRun.DER;
    DSSFileName = MATDSS.Sim.DSSFileName;
    MATDSS.Disturbance = [];
    MATDSS.Disturbance = struct;
    
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

    MATDSS.Meas.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Meas.TimeStep)==0));
    MATDSS.Sim.Meas.at = MATDSS.Time.Sim.TimeSpan;
    MATDSS.Cont.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Cont.TimeStep)==0));
    
    ControllersSettingsTableData = app.ControllersSettingsTable.Data.Variables;
    FixFields = {'Area', 'alpha', 'rp', 'rbard', 'E', 'vul', 'vll', 'iul', 'arho',...
        'asigma', 'alambda', 'amu', 'aeta', 'apsi', 'agamma', 'anu', 'azeta', 'crho',...
        'csigma', 'clambda', 'cmu', 'ceta', 'cpsi', 'cgamma', 'cnu', 'czeta'}';
    
    MATDSS.Cont.Gain = Gain;

    for i = 1:length(MATDSS.Cont.CA)
        ControllerFields = fieldnames(MATDSS.Cont.CA(i));
        
        [i_field] = MATDSS_StrComp(ControllerFields, FixFields);
        [...
            MATDSS.Cont.CA(i).Area,...
            MATDSS.Cont.CA(i).Type,...
            MATDSS.Cont.CA(i).alpha,...
            MATDSS.Cont.CA(i).rp,...
            MATDSS.Cont.CA(i).rbard,...
            MATDSS.Cont.CA(i).E,...
            MATDSS.Cont.CA(i).vul,...
            MATDSS.Cont.CA(i).vll,...
            MATDSS.Cont.CA(i).iul,...
            MATDSS.Cont.CA(i).arho,...
            MATDSS.Cont.CA(i).asigma,...
            MATDSS.Cont.CA(i).alambda,...
            MATDSS.Cont.CA(i).amu,...
            MATDSS.Cont.CA(i).aeta,...
            MATDSS.Cont.CA(i).apsi,...
            MATDSS.Cont.CA(i).agamma,...
            MATDSS.Cont.CA(i).anu,...
            MATDSS.Cont.CA(i).azeta,...
            MATDSS.Cont.CA(i).crho,...
            MATDSS.Cont.CA(i).csigma,...
            MATDSS.Cont.CA(i).clambda,...
            MATDSS.Cont.CA(i).cmu,...
            MATDSS.Cont.CA(i).ceta,...
            MATDSS.Cont.CA(i).cpsi,...
            MATDSS.Cont.CA(i).cgamma,...
            MATDSS.Cont.CA(i).cnu,...
            MATDSS.Cont.CA(i).czeta] = deal(ControllersSettingsTableData{i,:});


        % Changing strings to numbers in the corresponding locations.
        for j = 1:length(i_field)
            if j ~= [2,8]
                MATDSS.Cont.CA(i).(ControllerFields{i_field(j)}) = str2double(MATDSS.Cont.CA(i).(ControllerFields{i_field(j)}));
            end

            if j == 2 && ~strcmp(MATDSS.Cont.CA(i).(ControllerFields{i_field(j)}),'auto')
                MATDSS.Cont.CA(i).(ControllerFields{i_field(j)}) = str2double(MATDSS.Cont.CA(i).(ControllerFields{i_field(j)}));
            end

        end
        if strcmpi(MATDSS.Cont.CA(i).alpha, 'auto')
            MATDSS.Cont.CA(i).alpha = str2double(app.GeneralGainEditField.Value).*str2double(app.ContTimeStepEditField.Value)./1e3;
        else
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.CA(i).alpha.*str2double(app.GeneralGainEditField.Value).*str2double(app.ContTimeStepEditField.Value)./1e3;
        end


    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                   P0Set Function                  %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % P0Set Function
    P0SetFunctionType = app.P0SetFunctionButtonGroup.SelectedObject.Text;
    MATDSS.Sim.P0SetFunction.FunctionType = P0SetFunctionType; % UnitStep or Ramp
    MATDSS.Sim.P0SetFunction.Flag = 1;
    % MATDSS.Sim.P0SetFunction.FunctionType = 'UnitStep'; % UnitStep or Ramp
    % MATDSS.Sim.P0SetFunction.FunctionType = 'Ramp'; % UnitStep or Ramp



     % Read DER Table
    % Reading the MATDSS_DER.xlsx file to define the DERs
    % MATDSS.Table.DER = readtable([pwd '\Config\MATDSS_DER.xlsx'],'Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve');
    % MATDSS.TableData.DER = MATDSS.Table.DER.Variables;

    
else
    % MATDSS = MATDSS_Initialize
    % This function will initialize the simulation environment, set time
    % parameters, and initiate DSS Object. The function also defines the DER
    % struct that will be used in our circuits.
    % Defining Main Structs!

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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                     DER Struct                    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DER Struct
    % DER = struct; % details of DER struct are given in MATDSS_DERNew function.





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%              Compile OpenDSS Circuit              %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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







    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%      Reading the Configurations and Settings      %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reading All Excel Tables and save their values in MATDSS.Tables



    % Read DER Table
    % Reading the MATDSS_DER.xlsx file to define the DERs
    % MATDSS.Table.DER = readtable([pwd '\Config\MATDSS_DER.xlsx'],'Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve');
    % MATDSS.TableData.DER = MATDSS.Table.DER.Variables;





%     % Read VI Table
%     % Getting Voltage phases indecices for tracking
%     if exist([pwd '\Config\MATDSS_VI.xlsx'], 'file')
%         MATDSSVISheetnames =  cellstr(sheetnames('MATDSS_VI.xlsx')); % List of all sheets in the Excel file
%         MATDSSVIFlag = MATDSS_StrComp(MATDSSVISheetnames,app.OpenDSSFilesListBox.Value); % Check if teh target sheet is available in the Excel file
% 
%         if MATDSSVIFlag > 0
%             VIReadOptions = detectImportOptions('MATDSS_VI.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % set the import settings
%             VIReadOptions = setvartype(VIReadOptions,VIReadOptions.VariableNames(:),'string');
%             MATDSS.Table.VI = readtable([pwd '\Config\MATDSS_VI.xlsx'],VIReadOptions);%,'Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve');
%             MATDSS.TableData.VI = MATDSS.Table.VI.Variables;
%             clear VIReadOptions
%         end
%     end
%     %{
% Deprecated file configration
% % Read L2C Table
% % Read MATDSS_L2C.xlsx to get L2C Parameters
% L2CReadOptions = detectImportOptions('MATDSS_L2C.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % set the import settings
% L2CReadOptions = setvartype(L2CReadOptions,L2CReadOptions.VariableNames(:),'string');
% MATDSS.Table.L2C = readtable([pwd '\Config\MATDSS_L2C.xlsx'],L2CReadOptions);%'Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve',);
% MATDSS.TableData.L2C= MATDSS.Table.L2C.Variables;
% MATDSS.TableData.L2C = cellfun(@lower,MATDSS.TableData.L2C,'UniformOutput',false);
% 
%     %}
% 
% 
% 
% 
%     % Try reading MATDSS_ControlAreas.xlsx file
%     if exist([pwd '\Config\MATDSS_ControlAreas.xlsx'], 'file')
%         MATDSSCASettingsSheetnames =  cellstr(sheetnames('MATDSS_ControlAreas.xlsx')); % List of all sheets in the Excel file
%         MATDSSCASettingsFlag = MATDSS_StrComp(MATDSSCASettingsSheetnames,app.OpenDSSFilesListBox.Value); % Check if teh target sheet is available in the Excel file
%         if MATDSSCASettingsFlag > 0
%             % Read Control Areas Settings
%             CASettingsReadOptions = detectImportOptions('MATDSS_ControlAreas.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % set the import settings
%             CASettingsReadOptions = setvartype(CASettingsReadOptions,CASettingsReadOptions.VariableNames(:),'string');
%             MATDSS.Table.CASettings = readtable([pwd '\Config\MATDSS_ControlAreas.xlsx'],CASettingsReadOptions);
%             MATDSS.TableData.CASettings= MATDSS.Table.CASettings.Variables;
%             MATDSS.TableData.CASettings = cellfun(@lower,MATDSS.TableData.CASettings,'UniformOutput',false);
%         end
%     end
% 
% 
% 
% 
% 
% 
%     % Try reading MATDSS_ControllersSettings.xlsx file
%     if exist([pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'file')
%         MATDSSControllersSettingsSheetnames =  cellstr(sheetnames('MATDSS_ControllersSettings.xlsx')); % List of all sheets in the Excel file
%         MATDSSControllersSettingsFlag = MATDSS_StrComp(MATDSSControllersSettingsSheetnames,app.OpenDSSFilesListBox.Value); % Check if teh target sheet is available in the Excel file
%         if MATDSSControllersSettingsFlag > 0
%             % Read Control Areas Settings
%             ControllersSettingsReadOptions = detectImportOptions('MATDSS_ControllersSettings.xlsx','Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve'); % set the import settings
%             ControllersSettingsReadOptions = setvartype(ControllersSettingsReadOptions,ControllersSettingsReadOptions.VariableNames(:),'string');
%             MATDSS.Table.ControllersSettings = readtable([pwd '\Config\MATDSS_ControllersSettings.xlsx'],ControllersSettingsReadOptions);
%             MATDSS.TableData.ControllersSettings= MATDSS.Table.ControllersSettings.Variables;
%             MATDSS.TableData.ControllersSettings = cellfun(@lower,MATDSS.TableData.ControllersSettings,'UniformOutput',false);
%         end
%     end
end


MATDSS.Cont.Gain = Gain;
MATDSS.Cont.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Cont.TimeStep)==0));
MATDSS.Cont.RoU = str2double(app.RoUEditField.Value);
MATDSS.Cont.LoopSize = str2double(app.LoopsizeEditField.Value);
MATDSS.Cont.DisableControllersFlag = app.DisableControllersButton.Value;
MATDSS.Cont.s = app.TrackP0Button.Value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              Compile OpenDSS Circuit              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% app.StatusLabel.Text = 'MATDSS Initiallization Completed Successfully';


MATDSS.UseSimDataFlag = app.UseSimDataButton.Value;
end %function end


