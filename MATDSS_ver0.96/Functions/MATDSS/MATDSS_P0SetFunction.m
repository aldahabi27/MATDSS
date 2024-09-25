function MATDSS = MATDSS_P0SetFunction(MATDSS)
% MATDSS_P0SetFunction adjusts the P0 set points based on the configured 
% function type for the simulation.
%
% This function modifies the P0 set points in the MATDSS application 
% according to the specified function type, which can include ramp 
% changes or unit step changes. The function takes into account the 
% simulation time and adjusts the control signals accordingly.
%
% The function supports multiple function types:
% - 'Ramp': Gradually increases or decreases the set points based on 
%           the configured change per DER.
% - 'Unit Step': Changes the set points by a fixed amount when triggered.
%
% The function updates the ControlSignals.P0Set array with the new 
% values for each time step during the simulation.
%
% Parameters:
%   - MATDSS: A structure containing all necessary data and settings 
%             for the simulation, including time parameters and control 
%             signals.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

t = MATDSS.Sim.t; % Current time index for the simulation
at = MATDSS.Sim.at; % Actual time for the simulation

% Check if the P0SetFunction type is defined; if not, set it to 'Ramp'
if ~isfield(MATDSS.Sim.P0SetFunction, 'FunctionType')
    MATDSS.Sim.P0SetFunction.FunctionType = 'Ramp'; % Default function type
end

% Determine the action based on the function type specified
switch convertStringsToChars(lower(MATDSS.Sim.P0SetFunction.FunctionType))
    case 'ramp'
        P0SetChangePerDER = MATDSS.Sim.P0SetFunction.ChangekW; % Change rate per DER
        P0SetChange = P0SetChangePerDER; % Set change to the defined rate
        % Check conditions based on simulation time
        if at > MATDSS.Time.Sim.ST && mod(at, 1) == 0
            % Determine if the change should increase or decrease
            if at - MATDSS.Time.Sim.ST > 1 * (MATDSS.Time.Sim.TimeSpan(end) - MATDSS.Time.Sim.ST) / 2 && ...
               at - MATDSS.Time.Sim.ST < MATDSS.Time.Sim.TimeSpan(end) - MATDSS.Time.Sim.ST
                MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1) - P0SetChange; % Decrease
            elseif at - MATDSS.Time.Sim.ST < (MATDSS.Time.Sim.TimeSpan(end) - MATDSS.Time.Sim.ST) / 2
                MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1) + P0SetChange; % Increase
            else
                MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1); % No change
            end
        else
            MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1); % No change if conditions are not met
        end
    case 'unit step'
        P0SetChange = MATDSS.Sim.P0SetFunction.ChangekW; % Fixed change amount for unit step
        % Check if the conditions for the unit step are met
        if at >= MATDSS.Time.Sim.ST && MATDSS.Sim.P0SetFunction.Flag && mod(at, MATDSS.Time.Cont.TimeStep) == 0
            MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1) + P0SetChange; % Apply change
            MATDSS.Sim.P0SetFunction.Flag = 0; % Reset flag after applying change
        else
            MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1); % No change if conditions are not met
        end
    otherwise
        MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1); % Default case, no change
end 

end


%% Old function code

%{

function MATDSS = MATDSS_P0SetFunction (MATDSS)
t = MATDSS.Sim.t;
at = MATDSS.Sim.at;

if ~isfield(MATDSS.Sim.P0SetFunction,'FunctionType')
    MATDSS.Sim.P0SetFunction.FunctionType = 'Ramp';
end

switch convertStringsToChars(lower(MATDSS.Sim.P0SetFunction.FunctionType))
    case 'ramp'
        P0SetChangePerDER = MATDSS.Sim.P0SetFunction.ChangekW;
        P0SetChange = P0SetChangePerDER;% * MATDSS.Sim.nDER;
        if at > MATDSS.Time.Sim.ST && mod(at,1) == 0
            if at-MATDSS.Time.Sim.ST > 1*(MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST)/2 && at-MATDSS.Time.Sim.ST < MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST
                MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1) - P0SetChange;
            elseif at-MATDSS.Time.Sim.ST < (MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST)/2
                MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1) + P0SetChange;
            else
                MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1);
            end
        else
            MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1);
        end
    case 'unit step'
        P0SetChange = MATDSS.Sim.P0SetFunction.ChangekW; % 2MW
        if at >= MATDSS.Time.Sim.ST && MATDSS.Sim.P0SetFunction.Flag && mod(at,MATDSS.Time.Cont.TimeStep) == 0
            MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1) + P0SetChange;
            MATDSS.Sim.P0SetFunction.Flag = 0;
        else
            MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1);% - 100;
        end
    otherwise
        MATDSS.ControlSignals.P0Set(t) = MATDSS.ControlSignals.P0Set(t-1);
end

%}