function [MATDSS, DER] = MATDSS_Sim(app, MATDSS)
% MATDSS_Sim(app, MATDSS)
% This function performs the main simulation process for the MATDSS
% application. It establishes a link with OpenDSS and sets up the
% environment for time-domain dynamic simulation, including:
%
% 1. Initializing OpenDSS simulation engine.
% 2. Creating DER devices and disturbance loads.
% 3. Running the power flow solution in OpenDSS.
% 4. Extracting the YBus matrix from the simulated circuit.
% 5. Initializing time-based indices and preparing measurement data.
% 6. Defining and initializing controllers and updating control signals.
%
% This function manages the simulation flow from the setup phase to
% updating the DERs and applying disturbances. Live measurements and
% visual feedback can be updated in the app during the simulation.
%
% Parameters:
%   - app: The MATDSS application instance, containing GUI components.
%   - MATDSS: The structure holding simulation settings, circuit data,
%             and results.
%
% Returns:
%   - MATDSS: Updated MATDSS structure with results from the simulation.
%   - DER: Updated DER structure with dynamic responses.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com




if nargin < 2
    MATDSS = app; % If MATDSS is not provided, use 'app' as the main structure
    app = []; % Set 'app' to empty if it's not provided
end

% Check if the 'app' instance exists, and update the status to indicate
% that the simulation has started
if ~isempty(app)
    MATDSSApp_Status(app,'Simulation in Progress','Simulation Started');
    pause(0.1); % Slight pause for status update visibility
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               System Initialization               %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the simulation and set up the environment by calling
% MATDSS_Initialize function (currently commented out).

if MATDSS.Sim.DSSStartOk % If DSS is successfully loaded and linked

    % Update the status to indicate the creation of DER devices
    if ~isempty(app)
        MATDSSApp_Details(app, 'Creating DERs in OpenDSS');
        pause(0.1); % Slight pause for status update
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%               Creating DER Devices                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call the function to create DER devices and set them in the OpenDSS
    [MATDSS, DER] = MATDSS_DER(MATDSS);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%            Creating Disturbance Loads             %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate the creation of disturbance loads
    if ~isempty(app)
        MATDSSApp_Details(app, 'Creating Disturbance Loads in OpenDSS');
        pause(0.1); % Slight pause for status update
    end

    % Create disturbance loads in the OpenDSS circuit
    MATDSS = MATDSS_Disturbance(MATDSS);

    % Solve the power flow equations for 20 iterations to initialize
    for i = 1:20
        MATDSS.Sim.DSSSolution.Solve;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                    YBus Matrix                    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate that the Ybus matrix is being generated
    if ~isempty(app)
        MATDSSApp_Details(app, 'Generating Ybus Matrix');
        pause(0.1); % Slight pause for status update
    end

    % Generate the Ybus matrix for the current OpenDSS circuit
    MATDSS = MATDSS_YBusMatrix(MATDSS);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%            Time Index Initialization             %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize time index and setup for the simulation time span
    t = 1; % Initialize the time index at t = 1 (corresponding to 0 s)
    MATDSS.Sim.t = t; % Save the time index in the MATDSS structure
    at = MATDSS.Time.Sim.TimeSpan(t); % Get the actual time (in seconds)
    MATDSS.Sim.at = at; % Save the actual time in the MATDSS structure






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%           Initializing the Measurements           %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate measurement setup
    if ~isempty(app)
        MATDSSApp_Details(app, 'Setting up measurements');
        pause(0.1); % Slight pause for status update
    end

    % Initialize measurements for voltage, current, power, etc.
    MATDSS = MATDSS_Measurement(MATDSS);
    % Set initial control signals for power (P) and reactive power (Q)
    MATDSS.ControlSignals.P0Set = ones(length(MATDSS.Sim.Meas.at),1) .* sum(MATDSS.Meas.P0(:,t));
    MATDSS.ControlSignals.Q0Set = ones(length(MATDSS.Sim.Meas.at),1) .* sum(MATDSS.Meas.Q0(:,t));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%      Defining Controllers and Control Areas       %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate controller setup
    if ~isempty(app)
        MATDSSApp_Details(app, 'Setting up the controllers');
        pause(0.1); % Slight pause for status update
    end

    % Initialize the controllers and set control areas for DERs
    [MATDSS, DER] = MATDSS_Controllers(MATDSS, DER);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% STOPPED HERE STOPPED HERE %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Update the status to indicate initialization of controllers and DERs
    if ~isempty(app)
        MATDSSApp_Details(app, 'Initializing controllers and DERs');
        pause(0.1); % Slight pause for status update
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%            Initializing the Controllers            %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSS.Sim.RoU = 0; % Initialize the Rate of Update for the simulation (deprecated!)
    % Update controller setpoints and states
    [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS, DER);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                 DER Initialization                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update each DER device with its initialization values
    for i = 1:MATDSS.Sim.nDER
        [MATDSS, DER] = MATDSS_DERUpdate(MATDSS, DER, i);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%          Initializing Disturbance Loads           %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the disturbance loads to reflect initial conditions
    MATDSS = MATDSS_DisturbanceUpdate(MATDSS);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%             Live Measurements Updates             %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update live measurements in the app (if GUI is active)
    if ~isempty(app)
        MATDSSApp_LiveMeasurements(app, MATDSS, DER);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                  Simulation Start                 %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the count for continuous step button clicks based on time steps
    ContStButtonClickedCounter = MATDSS.Time.Cont.TimeStep / MATDSS.Time.Sim.TimeStep;

    % Initialize flag for continuous step button click details
    ContStButtonClickedDetailsFlag = 0;

    % Check if the application instance is available
    if ~isempty(app)
        % Update GUI status to indicate that the simulation is running
        MATDSSApp_Status(app, 'Simulation is Running!', 'Simulation is Running!')

        % Set the selected tab to the first tab in the main tab group
        app.TabGroup.SelectedTab = app.TabGroup.Children(1);

        % Set the selected tab to the second tab in the details tab group
        app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(2);
    end

    % Loop through each time step in the simulation time span
    for t = 2:length(MATDSS.Time.Sim.TimeSpan)
        % Check if the application instance is available
        if ~isempty(app)
            % Check if the stop button has been pressed
            if app.StopButton.Tag == '1'
                return; % Exit the simulation loop if stop is clicked
            end
        end

        % Increment the Run of Update (RoU) counter (deprecated)
        MATDSS.Sim.RoU = MATDSS.Sim.RoU + 1;

        % Update the time index for the simulation
        MATDSS.Sim.t = t;
        % Get the current simulation time from the time span
        at = MATDSS.Time.Sim.TimeSpan(t);
        % Store the current time in the MATDSS structure
        MATDSS.Sim.at = at;

        % Pause briefly to allow the GUI to process interactions
        pause(0.001);

        % Update the progress bar based on the current time step
        if (round(100 * t / length(MATDSS.Time.Sim.TimeSpan)) ~= round(100 * (t - 1) / length(MATDSS.Time.Sim.TimeSpan))) && ~isempty(app)
            % Update the progress bar display with the current time
            MATDSSApp_ProgressBar(t / length(MATDSS.Time.Sim.TimeSpan), app.ProgressBar1, ['t = ' num2str(round(at, 1)) ' s'], app.ProgressBarStatus);
            pause(0.001); % Brief pause to allow GUI update
        end

        % If the pause button is activated, keep looping until it is released
        if ~isempty(app)
            while app.PauseButton.Value
                app.StepButton.Enable = "on"; % Enable the step button
                app.ContStButton.Enable = "on"; % Enable the continuous step button
                pause(0.1); % Pause briefly to avoid busy-waiting

                % Check if the step button has been clicked
                if strcmpi(app.StepButton.Tag, "clicked")
                    app.StepButton.Tag = "not clicked"; % Reset step button tag
                    % Log the step event with current time information
                    MATDSSApp_Details(app, {'Step has been clicked.'; ['t = ' num2str(t - 1) ' (= ' num2str(MATDSS.Time.Sim.TimeSpan(t - 1)) 's)']; ' '});
                    break; % Exit the loop to continue simulation
                end

                % Check if the continuous step button has been clicked
                if strcmpi(app.ContStButton.Tag, "clicked")
                    ContStButtonClickedCounter = 0; % Reset the counter for continuous steps
                    ContStButtonClickedDetailsFlag = 1; % Set flag for detail logging
                end

                % Check if the continuous step button has been clicked within the allowed count
                if (ContStButtonClickedCounter < MATDSS.Time.Cont.TimeStep / MATDSS.Time.Sim.TimeStep)
                    ContStButtonClickedCounter = ContStButtonClickedCounter + 1; % Increment the counter
                    app.ContStButton.Tag = "not clicked"; % Reset continuous button tag
                    break; % Exit the loop to continue simulation
                else
                    % If the continuous step button was clicked, log the event
                    if ContStButtonClickedDetailsFlag
                        MATDSSApp_Details(app, {'Controller Step has been clicked.'; ['t = ' num2str(t - 1) ' (= ' num2str(MATDSS.Time.Sim.TimeSpan(t - 1)) 's)']; ' '});
                        ContStButtonClickedDetailsFlag = 0; % Reset the detail flag
                    end
                end
            end
        end

        % Disable step and continuous step buttons if pause is not active
        if ~isempty(app)
            if ~app.PauseButton.Value
                app.StepButton.Enable = "off"; % Disable step button
                app.ContStButton.Enable = "off"; % Disable continuous step button
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%              Update the Measurements              %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MATDSS = MATDSS_Measurement(MATDSS); % Update measurement data within MATDSS

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%                   Update P0_Set                   %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MATDSS = MATDSS_P0SetFunction(MATDSS); % Update P0 settings in MATDSS

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%               Update the Controllers               %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS, DER); % Update controllers based on the latest MATDSS data

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%               Show Live Measurements              %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(app)
            MATDSSApp_LiveMeasurements(app, MATDSS, DER); % Update the live measurement display in the GUI
            MATDSSApp_DERStatus(app, MATDSS, DER); % Update the DER status in the GUI
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%                  Update the DERs                  %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:MATDSS.Sim.nDER
            [MATDSS, DER] = MATDSS_DERUpdate(MATDSS, DER, i); % Update each DER based on MATDSS data
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%             Update Disturbance Loads              %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MATDSS = MATDSS_DisturbanceUpdate(MATDSS); % Update disturbance loads in MATDSS

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%          Debugging Code during run time           %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Commented out debugging code for troubleshooting purposes
        % if round(at-1,2) == 0 && false
        %     my_circut = MATDSS.Sim.DSSCircuit; % handle of DSS_Circuit
        %     myloads = my_circut.Loads; % handle of loads in DSS_Circuit
        %     myloads.Name = 'IlyasLoad'; % Set load name for debugging
        %     myloads.kW = 200; % Set load power for debugging
        % end

        % if mod(MATDSS.Sim.at,1) == 0 || (round(at-1.4,2) == 0)
        %     disp('hi'); % Debugging output
        % end
    end

end % End of the simulation loop
end




%% Old function code
%{
%{
% MATDSS_Sim(app, MATDSS)
% This function performs the main simulation process for the MATDSS 
% application. It establishes a link with OpenDSS and sets up the 
% environment for time-domain dynamic simulation, including:
%
% 1. Initializing OpenDSS simulation engine.
% 2. Creating DER devices and disturbance loads.
% 3. Running the power flow solution in OpenDSS.
% 4. Extracting the YBus matrix from the simulated circuit.
% 5. Initializing time-based indices and preparing measurement data.
% 6. Defining and initializing controllers and updating control signals.
%
% This function manages the simulation flow from the setup phase to 
% updating the DERs and applying disturbances. Live measurements and 
% visual feedback can be updated in the app during the simulation.
%
% Parameters:
%   - app: The MATDSS application instance, containing GUI components.
%   - MATDSS: The structure holding simulation settings, circuit data, 
%             and results.
%
% Returns:
%   - MATDSS: Updated MATDSS structure with results from the simulation.
%   - DER: Updated DER structure with dynamic responses.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com

function [MATDSS, DER] = MATDSS_Sim(app, MATDSS)

if nargin < 2
    MATDSS = app; % If MATDSS is not provided, use 'app' as the main structure
    app = []; % Set 'app' to empty if it's not provided
end

% Check if the 'app' instance exists, and update the status to indicate 
% that the simulation has started
if ~isempty(app)
    MATDSSApp_Status(app,'Simulation in Progress','Simulation Started');
    pause(0.1); % Slight pause for status update visibility
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               System Initialization               %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the simulation and set up the environment by calling
% MATDSS_Initialize function (currently commented out).

if MATDSS.Sim.DSSStartOk % If DSS is successfully loaded and linked

    % Update the status to indicate the creation of DER devices
    if ~isempty(app)
        MATDSSApp_Details(app, 'Creating DERs in OpenDSS');
        pause(0.1); % Slight pause for status update
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               Creating DER Devices                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call the function to create DER devices and set them in the OpenDSS
    [MATDSS, DER] = MATDSS_DER(MATDSS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            Creating Disturbance Loads             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate the creation of disturbance loads
    if ~isempty(app)
        MATDSSApp_Details(app, 'Creating Disturbance Loads in OpenDSS');
        pause(0.1); % Slight pause for status update
    end

    % Create disturbance loads in the OpenDSS circuit
    MATDSS = MATDSS_Disturbance(MATDSS);
    
    % Solve the power flow equations for 20 iterations to initialize
    for i = 1:20
        MATDSS.Sim.DSSSolution.Solve;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                    YBus Matrix                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate that the Ybus matrix is being generated
    if ~isempty(app)
        MATDSSApp_Details(app, 'Generating Ybus Matrix');
        pause(0.1); % Slight pause for status update
    end

    % Generate the Ybus matrix for the current OpenDSS circuit
    MATDSS = MATDSS_YBusMatrix(MATDSS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            Time Index Initialization             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize time index and setup for the simulation time span
    t = 1; % Initialize the time index at t = 1 (corresponding to 0 s)
    MATDSS.Sim.t = t; % Save the time index in the MATDSS structure
    at = MATDSS.Time.Sim.TimeSpan(t); % Get the actual time (in seconds)
    MATDSS.Sim.at = at; % Save the actual time in the MATDSS structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%           Initializing the Measurements           %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate measurement setup
    if ~isempty(app)
        MATDSSApp_Details(app, 'Setting up measurements');
        pause(0.1); % Slight pause for status update
    end

    % Initialize measurements for voltage, current, power, etc.
    MATDSS = MATDSS_Measurement(MATDSS);
    % Set initial control signals for power (P) and reactive power (Q)
    MATDSS.ControlSignals.P0Set = ones(length(MATDSS.Sim.Meas.at),1) .* sum(MATDSS.Meas.P0(:,t));
    MATDSS.ControlSignals.Q0Set = ones(length(MATDSS.Sim.Meas.at),1) .* sum(MATDSS.Meas.Q0(:,t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      Defining Controllers and Control Areas       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the status to indicate controller setup
    if ~isempty(app)
        MATDSSApp_Details(app, 'Setting up the controllers');
        pause(0.1); % Slight pause for status update
    end

    % Initialize the controllers and set control areas for DERs
    [MATDSS, DER] = MATDSS_Controllers(MATDSS, DER);

    % Update the status to indicate initialization of controllers and DERs
    if ~isempty(app)
        MATDSSApp_Details(app, 'Initializing controllers and DERs');
        pause(0.1); % Slight pause for status update
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%            Initializing the Controllers            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSS.Sim.RoU = 0; % Initialize the Rate of Update for the simulation
    % Update controller setpoints and states
    [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS, DER);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                 DER Initialization                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update each DER device with its initialization values
    for i = 1:MATDSS.Sim.nDER
        [MATDSS, DER] = MATDSS_DERUpdate(MATDSS, DER, i);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%          Initializing Disturbance Loads           %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the disturbance loads to reflect initial conditions
    MATDSS = MATDSS_DisturbanceUpdate(MATDSS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%             Live Measurements Updates             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update live measurements in the app (if GUI is active)
    if ~isempty(app)    
        MATDSSApp_LiveMeasurements(app, MATDSS, DER);
    end

end % End of DSSStartOk check


%}

function [MATDSS, DER] = MATDSS_Sim(app, MATDSS)
% [MATDSS, DER] = MATDSS_Sim(app, Gain)
% MATDSS_Sim function is the main function that run a time-based dynamic
% simulation. The function roles are:
%
% 1. Initialize the link with OpenDSS
% 2. Read Excel files and save their data in the corresponding table
%    variables in MATDSS.
% 3. Define OpenDSS loads to represent the DERs in the system.
% 4. Run a time domain simulation by calling OpenDSS engine to solve teh
% power flow equations, and update the output of DER devices following the
% set-points provided by the Level 2 Controller (L2C).
%
%

if nargin < 2
    MATDSS = app;
    app = [];
end

if ~isempty(app)
    MATDSSApp_Status(app,'Simulation in Progress','simulation Started')
    pause(0.1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%               System Initialization               %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [MATDSS] = MATDSS_Initialize(app,Gain, app.UseSimDataButton.Value);
%


if MATDSS.Sim.DSSStartOk % if DSS is loaded and linked successfully

    if ~isempty(app)
        MATDSSApp_Details(app, 'Creating DERs in OpenDSS')
        pause(0.1)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%               Creating DER Devices                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [MATDSS, DER] = MATDSS_DER(MATDSS);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%            Creating Disturbance Loads             %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(app)
        MATDSSApp_Details(app, 'Creating Disturbance Loads in OpenDSS')
        pause(0.1)
    end
    MATDSS = MATDSS_Disturbance(MATDSS);

    for i = 1:20
        MATDSS.Sim.DSSSolution.Solve;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                    YBus Matrix                    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(app)
        MATDSSApp_Details(app, 'Generating Ybus Matrix')
        pause(0.1)
    end
    MATDSS = MATDSS_YBusMatrix(MATDSS);


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
    if ~isempty(app)
        MATDSSApp_Details(app, 'Setting up measurements')
        pause(0.1)
    end
    MATDSS = MATDSS_Measurement(MATDSS);
    MATDSS.ControlSignals.P0Set = ones(length(MATDSS.Sim.Meas.at),1).*sum(MATDSS.Meas.P0(:,t));
    MATDSS.ControlSignals.Q0Set = ones(length(MATDSS.Sim.Meas.at),1).*sum(MATDSS.Meas.Q0(:,t));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%      Defining Controllers and Control Areas       %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(app)
        MATDSSApp_Details(app, 'Setting up the controllers')
        pause(0.1)
    end
    [MATDSS, DER] = MATDSS_Controllers(MATDSS,DER); %Initialize the controllers



    if ~isempty(app)
        MATDSSApp_Details(app, 'Initiallizing controllers and DERs')
        pause(0.1)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%            Initializing the Controlles            %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MATDSS.Sim.RoU = 0;
    [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS,DER);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                 DER Initialization                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if ~isempty(app)
        MATDSSApp_LiveMeasurements(app,MATDSS,DER);
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                  Simulation Start                 %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ContStButtonClickedCounter = MATDSS.Time.Cont.TimeStep/MATDSS.Time.Sim.TimeStep;
    ContStButtonClickedDetailsFlag = 0;
    if ~isempty(app)
        MATDSSApp_Status(app, 'Simulation is Running!', 'Simulation is Running!')
        app.TabGroup.SelectedTab = app.TabGroup.Children(1);
        app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(2);

    end
    % MATDSSApp_Status(app,'Simulation is Running!')
    for t = 2:length(MATDSS.Time.Sim.TimeSpan)
        % If stop is clicked, stop the simulation
        if ~isempty(app)
            if app.StopButton.Tag == '1'
                return;
            end
        end

        MATDSS.Sim.RoU = MATDSS.Sim.RoU + 1;

        % disp('1')



        % Update time index
        MATDSS.Sim.t = t;
        at = MATDSS.Time.Sim.TimeSpan(t);
        MATDSS.Sim.at = at;
        % app.MATDSSRunVariables.t = t;
        % app.MATDSSRunVariables.at = at;




        pause(0.001); % This is needed to allow the GUI to interrupt the code if a button status changes (stop, pause, live meas. etc.)
        % Progress bar update
        if (round(100*t/length(MATDSS.Time.Sim.TimeSpan)) ~= round(100*(t-1)/length(MATDSS.Time.Sim.TimeSpan))) && ~isempty(app)
            MATDSSApp_ProgressBar(t/length(MATDSS.Time.Sim.TimeSpan),app.ProgressBar1,['t = ' num2str(round(at,1)) ' s'], app.ProgressBarStatus);
            pause(0.001);
        end
        %
        % if at <= 0.1 || (at >= 5 && at <= 5.5)
        %     [MATDSS.Sim.Meas.PProfile(1:3,1), MATDSS.Sim.Meas.PProfile(1:3,end-1)]
        %     MATDSS.Sim.Meas.PProfile(1:3,1) - MATDSS.Sim.Meas.PProfile(1:3,end-1)
        %     sum(MATDSS.Sim.Meas.PProfile(1:3,1) - MATDSS.Sim.Meas.PProfile(1:3,end-1))./1e3
        %     disp('time time time');
        % end
        %

        % If pause is "ON", then keep looping here :D
        if ~isempty(app)
            while app.PauseButton.Value
                app.StepButton.Enable = "on";
                app.ContStButton.Enable = "on";
                pause(0.1);
                if strcmpi(app.StepButton.Tag,"clicked")
                    app.StepButton.Tag = "not clicked";
                    MATDSSApp_Details(app,{'Step has been clicked.';['t = ' num2str(t-1) ' (= ' num2str(MATDSS.Time.Sim.TimeSpan(t-1)) 's)'];' '});
                    break;
                end
                if strcmpi(app.ContStButton.Tag, "clicked")
                    ContStButtonClickedCounter = 0;
                    ContStButtonClickedDetailsFlag = 1;
                end

                if (ContStButtonClickedCounter < MATDSS.Time.Cont.TimeStep/MATDSS.Time.Sim.TimeStep)
                    ContStButtonClickedCounter = ContStButtonClickedCounter + 1;
                    app.ContStButton.Tag = "not clicked";
                    break;
                else
                    if ContStButtonClickedDetailsFlag
                        MATDSSApp_Details(app,{'Controller Step has been clicked.';['t = ' num2str(t-1) ' (= ' num2str(MATDSS.Time.Sim.TimeSpan(t-1)) 's)'];' '});
                        ContStButtonClickedDetailsFlag = 0;
                    end
                end
            end
        end
        %
        if ~isempty(app)
            if ~app.PauseButton.Value
                app.StepButton.Enable = "off";
                app.ContStButton.Enable = "off";
            end
        end
        % disp('2')



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%              Update the Measurements              %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MATDSS = MATDSS_Measurement(MATDSS);
        % disp('3')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%                   Update P0_Set                   %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MATDSS = MATDSS_P0SetFunction (MATDSS);
        % disp('4')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%               Update the Controlles               %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS,DER);
        % disp('5')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%               Show Live Measurements              %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(app)
            MATDSSApp_LiveMeasurements(app,MATDSS,DER);
            % disp('6')
            MATDSSApp_DERStatus(app,MATDSS,DER);
            % disp('7')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%                  Update the DERs                  %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:MATDSS.Sim.nDER
            [MATDSS, DER] = MATDSS_DERUpdate(MATDSS,DER,i);
        end
        % disp('8')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%             Update Disturbance Loads              %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MATDSS = MATDSS_DisturbanceUpdate(MATDSS);

        % disp('9')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%          Debugging Code during run time           %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         if round(at-1,2) == 0 && false
        %             my_circut = MATDSS.Sim.DSSCircuit; % handle of DSS_Circuit
        %             myloads = my_circut.Loads; % handle of loads in DSS_Circuit
        %             myloads.Name = 'IlyasLoad';
        %             myloads.kW = 200;
        %         end
        %
        %         if mod(MATDSS.Sim.at,1) == 0 || (round(at-1.4,2) == 0)
        % %             disp('hi');
        %         end
        %         [DER(1).kWSetpoint(end), DER(1).kW(end), abs(DER(1).kWSetpoint(end)-DER(1).kWSetpoint(end-1));...
        %          DER(1).kvarSetpoint(end), DER(1).kvar(end),abs(DER(1).kvarSetpoint(end)-DER(1).kvarSetpoint(end-1))]


        % disp('10')
    end

end





end

%}