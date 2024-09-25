function MATDSS = MATDSS_DisturbanceUpdate(MATDSS)
% MATDSS_DisturbanceUpdate(MATDSS)
% This function updates the status of disturbances within the MATDSS simulation 
% based on their configuration and the current simulation time.
% It interacts with the OpenDSS circuit to adjust load characteristics dynamically 
% during the simulation run.
%
% The function loops through all defined disturbances and adjusts their 
% power output if they are enabled and within their active time window. 
% If a disturbance is active, its specified power (kW) is applied to the 
% corresponding load in OpenDSS. Otherwise, the load's power is set to 0.
%
% Parameters:
%   - MATDSS: The MATDSS application instance containing simulation properties, 
%             including the disturbance definitions, time, and OpenDSS circuit.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com
%




% This function updates the status of disturbances according to their 
% configurations and the current simulation time. It interacts with the 
% OpenDSS circuit by updating the corresponding loads' power output 
% (in kilowatts) based on the disturbance schedule and settings.
%

% Setting local variables for simulation time
t = MATDSS.Sim.t;   % Current simulation time step
at = MATDSS.Sim.at; % Absolute time in the simulation

% Retrieving the number of disturbances to update
nDis = MATDSS.Sim.nDis; % Number of active disturbances in the simulation

% Loop over each disturbance to check its status and apply updates
for i = 1:nDis
    
    % Obtain a handle to the OpenDSS circuit and load elements
    my_circut = MATDSS.Sim.DSSCircuit; % Handle for the DSS_Circuit object
    myloads = my_circut.Loads; % Handle for the DSS loads object
    
    % Select the specific load associated with the disturbance using its DSS name
    myloads.Name = MATDSS.Disturbance(i).DSSName; % Name of the load (disturbance)
    
    % Check if the disturbance is active based on its time window and enabled status
    if ((at - MATDSS.Time.Sim.ST) >= MATDSS.Disturbance(i).t(1) && (at - MATDSS.Time.Sim.ST) <= MATDSS.Disturbance(i).t(2)) ...
            && MATDSS.Disturbance(i).Enabled
        
        % If the disturbance is active and enabled, set the power (kW) of the load
        myloads.kW = MATDSS.Disturbance(i).P; % Apply the disturbance's active power setting
        
        % (Optional) Code to modify power factor can be uncommented if needed
        % myloads.PF = MATDSS.Disturbance(i).pf; % Set the power factor (if defined)
    else
        % If the disturbance is outside its active time window or disabled, set the load to zero power
        myloads.kW = 0; % Set power to 0 kW when disturbance is not active
    end
    
end

end


%% Old function code

%{

function MATDSS = MATDSS_DisturbanceUpdate(MATDSS)
% MATDSS = MATDSS_DisturbanceUpdate(MATDSS)
%
% This function will update the status of the disturbance according to its
% configurations.
%
%


% Setting some easy callback variables
t = MATDSS.Sim.t;
at = MATDSS.Sim.at;


nDis = MATDSS.Sim.nDis; % Number of Disturbances

for i = 1:nDis
    my_circut = MATDSS.Sim.DSSCircuit; % handle of DSS_Circuit
    myloads = my_circut.Loads; % handle of loads in DSS_Circuit
    myloads.Name = MATDSS.Disturbance(i).DSSName; %Select the DER to be updated

    if ((at-MATDSS.Time.Sim.ST) >= MATDSS.Disturbance(i).t(1) && (at-MATDSS.Time.Sim.ST) <= MATDSS.Disturbance(i).t(2)) && MATDSS.Disturbance(i).Enabled
        myloads.kW = MATDSS.Disturbance(i).P;
%         myloads.PF = MATDSS.Disturbance(i).pf;
    else
        myloads.kW = 0;
    end
end


end

%}