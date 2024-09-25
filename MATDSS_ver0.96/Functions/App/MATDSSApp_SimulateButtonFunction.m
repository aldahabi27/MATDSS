function MATDSSApp_SimulateButtonFunction(app)
% MATDSSApp_SimulateButtonFunction(app)
% This function serves as the core of the MATDSS application.
% It initiates the simulation process while locking the GUI to prevent
% any modifications to the configurations during runtime. The function 
% prepares the Run variable, which contains all simulation outputs 
% used in the editable plots.
%
% Error handling is implemented to catch any issues that may occur 
% during the simulation. If an error is encountered, it will be 
% displayed in the application details text area, allowing for 
% easy troubleshooting.
%
% The primary responsibility of this function is to call the 
% MATDSS_Sim function, which executes the simulation based on the 
% current settings. Additionally, it updates the progress bar for 
% multiple runs, although this feature may be deprecated in future 
% versions (in which case this function may be combined with 
% MATDSS_Sim).
%
% Parameters:
%   - app: The MATDSS application instance containing properties and 
%          methods related to the simulation and GUI management.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com

% Switch to the first tab in the TabGroup and DetailsDERSTabGroup
app.TabGroup.SelectedTab = app.TabGroup.Children(1);
app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(1);

% Clear any saved data in MyRun variables
app.MyRun.Initialization = [];
app.MyRun.SimulationCompleted = false;
app.MyRun.MATDSS = [];
app.MyRun.DER = [];

% Lock the GUI and reset the progress bars
MATDSSApp_ProgressBar(0, app.ProgressBar1, 't = 0 s', app.ProgressBarStatus);
app.StopButton.Tag = '0';
app.PauseButton.Tag = '0';

% Update the details area to indicate simulation start
MATDSSApp_Details(app, {''; ['Simulation started for "' app.OpenDSSFilesListBox.Value '". Please wait!']});

% Use try-catch to handle potential errors during simulation
try
    % Prepare the Gains array from the General Gain edit field
    if app.GeneralGainEditField.Value(end) == ';'
        app.GeneralGainEditField.Value(end) = []; 
    end
    Gains = cellfun(@str2double, strsplit(app.GeneralGainEditField.Value, ';'));

    % Loop over the Gains to perform multiple simulations
    for iGain = 1:length(Gains) 
        Gain = Gains(iGain);  % Set the current gain value for simulation
        
        % Initialize the application for the current gain value
        MATDSSApp_Initialization(app, Gain);

        % Call the MATDSS_Sim function to run the simulation with the given Gain value
        [MATDSS_iGain, DER_iGain] = MATDSS_Sim(app, app.Main.MATDSS);

        % Save the simulation output in MyRun handle for future reference
        app.MyRun.MATDSS = [MATDSS_iGain];
        app.MyRun.DER = [DER_iGain];
        app.MyRun.SimulationCompleted = true;

        % Check if the Stop button was pressed to break the loop
        if app.StopButton.Tag ~= '0'
            break;
        end
    end

% Catch any errors and display the details in the MATDSS Details tab
catch MATDSS_Error
    ErrorString = [];
    ErrorStack = MATDSS_Error.stack;
    for i = 1:length(ErrorStack)
        ErrorFiles = ErrorStack(i).file;
        ErrorLines = ErrorStack(i).line;
        ErrorString = [ErrorString; strcat(ErrorFiles, " -> Line ", num2str(ErrorLines))];
    end
    ErrorMsg = sprintf('%s\n\n%s\n%s\n', ...
                       MATDSS_Error.identifier, ...
                       MATDSS_Error.message, ...
                       strjoin(ErrorString, '\n'));
    
    % Display error message in the details area
    MATDSSApp_Details(app, {ErrorMsg; ' '; '******************************************************'; ' '});
    assignin('base', 'MATDSS_Error', MATDSS_Error);
    app.StatusLabel.Text = 'Error occurred!';
    
    % Update details area with error information
    MATDSSApp_Details(app, {'Error in MATDSSApp_SimulateButtonFunction occurred. Check the details above for more information'; 'The program will stop the simulation.'});
    app.SimulateButton.Enable = "on";  % Re-enable the simulate button
end

% Update status label if no error occurred
if ~strcmp(app.StatusLabel.Text, 'Error occurred!')
    app.StatusLabel.Text = 'Ready!';
    MATDSSApp_Details(app, {'Simulation completed successfully'; ''});
end

% Clear variables except for the app instance
clearvars -except app

% Pass the MyRun variable to the workspace after simulations
if exist('MATDSS_iGain', 'var')
    assignin('base', 'MyRun', app.MyRun);
end

% In case of an error, we can still check the details by examining the "app" variable
assignin('base', 'app', app);


%% Old function code

%{
function MATDSSApp_SimulateButtonFunction(app)
% MATDSSApp_SimulateButtonFunction(app)
% This function is the core of MATDSS Application.
% The function will lock the GUI to prevent any changes in the
% configurations while the simulation is running. The function will also
% prepare Run variable that will contain all simulation outputs that will
% be used in the editable plots.
%
%
% The function also will handle all errors that might occure during the
% simulation and will display the errors in the application details text
% area.
%
% The main role of this function is call MATDSS_Sim to run the simulation
% for the current settings. The function also updates the progress bar
% related to multiple runs. This feature might be deprecated in future
% versions. (in that case, this function will be combined with MATDSS_Sim)
app.TabGroup.SelectedTab = app.TabGroup.Children(1);
app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(1);


% Clear any saved data in MyRun variables.
app.MyRun.Initialization = [];
app.MyRun.SimulationCompleted = false;
app.MyRun.MATDSS = [];
app.MyRun.DER = [];
% app.StatusLabel.FontColor = [0,0,0];


% Lock the GUI and reset the progress bars, stop and paused buttons.
% MATDSSApp_GUIChanges(app, 'Disable_GUI');
MATDSSApp_ProgressBar(0,app.ProgressBar1,'t = 0 s',app.ProgressBarStatus);
app.StopButton.Tag = '0';
app.PauseButton.Tag = '0';
% app.SaveButton.Enable = 'off';
% app.StepButton.Enable = "off";
% app.ContStButton.Enable = "off";
% app.ConfigurationButton.Enable = 'off';
MATDSSApp_Details(app,{'';['Simulation started for "' app.OpenDSSFilesListBox.Value '". Please wait!']});

% I use try and catch to direct error messages to the MATDSS details tab.
% This is important if we are to keep the variables and all program
% excusions active and not brake the flow.
try
    % This part is to detect how many Genera gain values are given to run
    % multiple simulations as a result.
    if app.GeneralGainEditField.Value(end) == ';'; app.GeneralGainEditField.Value(end) = []; end
    Gains = cellfun(@str2double,strsplit(app.GeneralGainEditField.Value,';'));

    for iGain = 1:length(Gains) % loop over the gains
        % Multiple runs progress bar settings
%         app.StatusLabel.Text = ['Gain = ' num2str(Gains(iGain))];
        
        
        % Set the current gain value to run the simulation for
        Gain = Gains(iGain);
        
        MATDSSApp_Initialization(app,Gain);

        [MATDSS_iGain, DER_iGain] = MATDSS_Sim(app,app.Main.MATDSS); % Call MATDSS_Sim Function to start the simulation with the given Gain value


        % Save the simulation output in MyRun handle for reference if
        % needed later.
        app.MyRun.MATDSS = [MATDSS_iGain];
        app.MyRun.DER = [DER_iGain];
        app.MyRun.SimulationCompleted = true;
        % app.StepButton.Enable = "off";
        % app.ContStButton.Enable = "off";
        if app.StopButton.Tag ~= '0'
            break;
        end
    end
    

    % disp('hi');

    
% Catch the erro and display all the details in the MATDSS Details tab
catch MATDSS_Error
    ErrorString = [];
    ErrorStack = MATDSS_Error.stack;
    for i = 1:length(ErrorStack)
        ErrorFiles = ErrorStack(i).file;
        ErrorLines = ErrorStack(i).line;

        ErrorString = [ErrorString;strcat(ErrorFiles," -> Line ", num2str(ErrorLines))];
    end
    ErrorMsg = sprintf('%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', MATDSS_Error.identifier,MATDSS_Error.message,ErrorString);
    %                 app.DetailsTextArea.Value = [app.DetailsTextArea.Value; ErrorMsg;' ';'******************************************************';' '];
    MATDSSApp_Details(app,{ErrorMsg;' ';'******************************************************';' '});
    assignin('base','MATDSS_Error',MATDSS_Error);
    app.StatusLabel.Text = 'Error occured!';
    % app.StatusLabel.FontColor = app.MyRun.Plot.Colors(12,:);
    MATDSSApp_Details(app,{'Error in MATDSSApp_SimulateButtonFunction occured. Check the details above for more information';'The program will stop the simulation.'});
    app.SimulateButton.Enable = "on";
end
if ~strcmp(app.StatusLabel.Text,'Error occured!')
    app.StatusLabel.Text = 'Ready!';
    MATDSSApp_Details(app,{'Simulation completed successfully';''});
end
clearvars -except app
% Re-enable the GUI
% MATDSSApp_GUIChanges(app, 'Enable_GUI');

% Pass the the MyRun variable to the workspace at the end of the
% Simulations.
if exist('MATDSS_iGain','var')
    assignin('base','MyRun',app.MyRun);
end

% In case error happend, we still can check the details by looking at the
% "app" variable which would contain the data.
assignin('base','app',app);



%}