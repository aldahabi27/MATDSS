function MATDSS = MATDSS_Measurement(MATDSS)
% MATDSS_Measurement(MATDSS)
% This function is responsible for updating the measurements in the MATDSS
% simulation. It calls the measurement functions VProfile and IPQ, which 
% update the voltage profile and active/reactive power profiles, respectively.
% Additionally, the function is designed to handle the delay index for 
% measurements, although the delay mechanism is currently hardcoded to 0.
%
% When the simulation time step 't' equals 1, the function performs a power 
% flow solution twice in OpenDSS. The voltage and power profiles are updated
% on each time step. The function may also display progress (e.g., using disp),
% but those statements are currently commented out.
%
% Parameters:
%   - MATDSS: The MATDSS application instance, containing properties related
%             to the simulation, including time step and solution management.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Solve the power flow at time step 1
if MATDSS.Sim.t == 1
    % Solve the OpenDSS power flow twice for initialization purposes
    MATDSS.Sim.DSSSolution.Solve;
    MATDSS.Sim.DSSSolution.Solve;
end

% Solve the power flow for the current time step
MATDSS.Sim.DSSSolution.Solve;

% Obtain the voltage profile (VProfile) using the MATDSS_VProfile function
MATDSS = MATDSS_VProfile(MATDSS);

% Obtain the power profiles (P, Q, P0, Q0) using the MATDSS_IPQ function
MATDSS = MATDSS_IPQ(MATDSS);

% Uncommented hardcoded delay value (currently set to 0 for no delay)
% MeasDelay = 0; % in index delay for now (multiples of MeasTimeStep)

% Uncommented code for handling delay in measurement queue index
% MATDSS.Time.Meas.QueueIndex = MATDSS.Sim.t - ...
%     MeasDelay*MATDSS.Time.Meas.TimeStep/MATDSS.Time.SimTimeStep;

end


%% Old function Code
% 
%{

function MATDSS = MATDSS_Measurement(MATDSS)
% MATDSS = MATDSS_Measurement(MATDSS)
% This function will call the measurement functions VProfile and IPQ to
% update the measurements. In addition, the function will update the delay
% index for measurements (updating the index pointing to the queue)

% Solve the power flow
if MATDSS.Sim.t == 1
    MATDSS.Sim.DSSSolution.Solve;
    MATDSS.Sim.DSSSolution.Solve;
end
% if MATDSS.Sim.t >= 14
%     disp('h');
% end
    MATDSS.Sim.DSSSolution.Solve;
    % disp('a')
MATDSS = MATDSS_VProfile(MATDSS); %obtain the voltage profile
% disp('b')
MATDSS = MATDSS_IPQ(MATDSS); % obtain the P, Q, P0 and Q0 profiles
% disp('c')
% Hardcoding the delay now
% MeasDelay = 0; % in index delay for now (multiples of MeasTimeStep)

% MATDSS.Time.Meas.QueueIndex = MATDSS.Sim.t -...
%     MeasDelay*MATDSS.Time.Meas.TimeStep/MATDSS.Time.SimTimeStep;



end

%}