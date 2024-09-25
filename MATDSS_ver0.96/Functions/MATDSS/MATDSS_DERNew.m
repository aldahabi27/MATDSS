function [MATDSS, DER] = MATDSS_DERNew(MATDSS, DER, DERInd, Name, Bus, Tau, DERType, Nodes, ConnType, Nphase, Mode, Px, Qx, Pmin, Pmax, Qmin, Qmax, ax, cx)
% MATDSS_DERNew Function
% This function defines a new Distributed Energy Resource (DER) in the MATDSS
% application. It interfaces with OpenDSS to create the DER load based on
% the provided parameters and manages dynamic updates within the DER
% structure.
%
% Input:
%   MATDSS    - Main MATDSS structure containing simulation variables and configurations.
%   DER       - Structure storing information about created DERs.
%   DERInd    - Index for tracking different DER devices.
%   Name      - Name of the DER device.
%   Bus       - Bus name or number where the DER is defined.
%   Tau       - Time constant for DER dynamics.
%   DERType   - Type of DER (e.g., VDER, conventional).
%   Nodes     - Nodes where the DER is connected.
%   ConnType  - Connection type (Wye or Delta).
%   Nphase    - Number of phases.
%   Mode      - Type of OpenDSS representation (currently supports load type).
%   Px        - Coefficients of P(x) in cost function.
%   Qx        - Coefficients of Q(x) in cost function.
%   Pmin      - Minimum active power limit of the DER.
%   Pmax      - Maximum active power limit of the DER.
%   Qmin      - Minimum reactive power limit of the DER.
%   Qmax      - Maximum reactive power limit of the DER.
%   ax, cx    - Scaling constants for step size in optimization.
%
% Output:
%   MATDSS    - Updated MATDSS structure.
%   DER       - Updated DER structure with new DER information.
%
% Notes:
% - If the DERType is 'VDER', no OpenDSS load definition occurs, and DER dynamics are not initialized.
% - The function initializes DER parameters, defines the load in OpenDSS, and updates the DER structure.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com



% Reading and converting the input variables from "string" to corresponding data type.
DERInd = str2double(DERInd); % Convert DERInd to double
Name = convertStringsToChars(Name); % Convert Name to character array
Bus = convertStringsToChars(Bus); % Convert Bus to character array
Tau = str2double(Tau); % Convert Tau to double
DERType = convertStringsToChars(DERType); % Convert DERType to character array
Nodes = str2num(Nodes); % Convert Nodes to numeric array
ConnType = convertStringsToChars(ConnType); % Convert ConnType to character array
Nphase = str2double(Nphase); % Convert Nphase to double
Mode = convertStringsToChars(Mode); % Convert Mode to character array
Px = str2num(Px); % Convert Px to numeric array
Qx = str2num(Qx); % Convert Qx to numeric array
Pmin = str2double(Pmin); % Convert Pmin to double
Pmax = str2double(Pmax); % Convert Pmax to double
Qmin = str2double(Qmin); % Convert Qmin to double
Qmax = str2double(Qmax); % Convert Qmax to double
ax = str2double(ax); % Convert ax to double
cx = str2double(cx); % Convert cx to double

timespan = MATDSS.Sim.Meas.at; % Time span for measurements
setpoint_timespan = MATDSS.Cont.at; % Setpoint time span

% Getting bus number (check if bus name is given)
if Bus(1) == "'" % Check if Bus is a string of bus name
    BusNum = MATDSS_StrComp(MATDSS.Sim.DSSCircuit.AllBusNames, Bus(2:end-1)); % Find bus number
else % Bus is a string of bus number
    BusNum = str2double(Bus); % Convert Bus to double
end

% Select the bus in OpenDSS to define the DER OpenDSS load
myBus = MATDSS.Sim.DSSCircuit.Buses(int16(BusNum-1)); % Select bus by index

% Prepare NodesText for node connection string
NodesText = [];
myNodes = Nodes;
while ~isempty(myNodes)
    NodesText = [NodesText, '.', num2str(myNodes(1))]; % Concatenate nodes to NodesText
    myNodes = myNodes(2:end); % Move to the next node
end

% Check if DERType is 'vder'; if so, skip OpenDSS load definition and DER dynamics
if ~strcmpi(DERType, 'vder')
    % Define DER in OpenDSS
    DSSNewDERCommand = ['New Load.', Name, '  Bus1=', myBus.Name, NodesText, '  Phases=', num2str(Nphase), ...
        '  Conn=', ConnType, '  Model=1 kV=', num2str(myBus.kVBase * sqrt(3)), '  kW=0  kvar=0 Vminpu=0.1']; % OpenDSS command

    % Prepare the DER struct to handle DER dynamics
    DER(DERInd).x.ODETime = []; % Initialize ODE time
    DER(DERInd).x.ODEx = []; % Initialize ODE state
    MATDSS.Sim.DSSText.Command = DSSNewDERCommand; % Send command to OpenDSS
    DER(DERInd).Tau = Tau; % Set DER time constant
end

% Save all values in the DER struct
DER(DERInd).DER_Ind = DERInd; % DER index
DER(DERInd).DSSName = Name; % DER name in OpenDSS
DER(DERInd).BusNum = BusNum; % Bus number
DER(DERInd).BusName = myBus.Name; % Bus name in OpenDSS
DER(DERInd).Type = DERType; % DER type
DER(DERInd).Nodes = Nodes; % Nodes connected to DER
DER(DERInd).ConnType = ConnType; % Connection type
DER(DERInd).NPhases = Nphase; % Number of phases
DER(DERInd).Mode = Mode; % OpenDSS mode representation
DER(DERInd).Px = Px; % Coefficients of P(x) in cost function
DER(DERInd).Qx = Qx; % Coefficients of Q(x) in cost function
DER(DERInd).Pmin = Pmin; % Minimum active power limit
DER(DERInd).Pmax = Pmax; % Maximum active power limit
DER(DERInd).Qmin = Qmin; % Minimum reactive power limit
DER(DERInd).Qmax = Qmax; % Maximum reactive power limit
DER(DERInd).ax = ax; % Scaling constant ax
DER(DERInd).cx = cx; % Scaling constant cx
DER(DERInd).MultiPhaseRatio = 1/Nphase; % Multi-phase ratio
DER(DERInd).W = [timespan', zeros(length(timespan), 1)]; % Initialize W array
DER(DERInd).Var = [timespan', zeros(length(timespan), 1)]; % Initialize Var array
DER(DERInd).i_WVar = 1; % Initialize index for W and Var
DER(DERInd).WSetpoint = [setpoint_timespan', zeros(length(setpoint_timespan), 1)]; % Initialize WSetpoint
DER(DERInd).VarSetpoint = [setpoint_timespan', zeros(length(setpoint_timespan), 1)]; % Initialize VarSetpoint
DER(DERInd).WUnfilteredSetpoint = zeros(length(setpoint_timespan), 1); % Initialize WUnfilteredSetpoint
DER(DERInd).VarUnvfilteredSetpoint = zeros(length(setpoint_timespan), 1); % Initialize VarUnvfilteredSetpoint
DER(DERInd).i_setpoint = 1; % Initialize index for setpoint
DER(DERInd).UpdateStatus = 1; % Update status flag
DER(DERInd).LLCMoveSetPoint = 0; % Flag for moving setpoint due to constraints

% Update the number of DERs in the Sim struct
MATDSS.Sim.nDER = MATDSS.Sim.nDER + 1; % Increment number of DERs

% clearvars -except app MATDSS DER % Clear unnecessary variables

end



%% Old function code
% 
%{ 

function [MATDSS, DER] = MATDSS_DERNew(MATDSS,DER,DERInd,Name,Bus,Tau,DERType,Nodes,ConnType,Nphase,Mode,Px,Qx,Pmin,Pmax,Qmin,Qmax,ax,cx)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Update documentation   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


This function will define a new load in the circuit following the
properties passed as an input.

Input: 

MATDSS      -> The main MATDSS struct which contains all variables related to
               the current simulation. It is important to establish OpenDSS link prior to
               calling MATDSS_DERNew function.

DER         -> DER struct where the new DER information will be saved.

DERInd      -> DERInd is an index that the user can give to the DER device.
               It will be used to track different DER devices later and map
               them to their corresponding control areas.

Name        -> Name of the DER device.

Bus         -> Bus name or number where the new load is to be defined.

Tau         -> Time constant used in modeling the dynamics of DER

DER_Type    -> The list of possible types of DERs (could affect the cost function
               and other related parameters, refer to the types written in the
               MATDSS_Main) (*Note: Currently this is not used to change
               the corresponding cost function)

Nodes       -> The nodes where the load will be connected

Conn_type   -> Wye or Delta connection

Nphase      -> number of phases

Mode        -> Define the type of OpenDSS representation. Currently only supports
               load type

Px          -> Coefficients of P(x) in the cost function

Qx          -> Coefficients of Q(x) in the cost function

Pmin, Pmax, Qmin, Qmax -> are the P and Q limits of the DER.

ax, cx      -> are the scaling constants of the general step size alpha and
               c

%}


% Fix bus voltage later, check the following in opendss excel
%{
Property kV As Double
    Member of OpenDSSengine.Loads
    Set kV rating for active Load. For 2 or more phases set Line-Line kV. Else actual kV across terminals.
%}

% Reading and converting the input variables from "string" to corresponding
% data type.
DERInd = str2double(DERInd);
Name = convertStringsToChars(Name);
Bus = convertStringsToChars(Bus);
Tau = str2double(Tau);
DERType = convertStringsToChars(DERType);
Nodes = str2num(Nodes);
ConnType = convertStringsToChars(ConnType);
Nphase = str2double(Nphase);
Mode = convertStringsToChars(Mode);
Px = str2num(Px);
Qx = str2num(Qx);
Pmin = str2double(Pmin);
Pmax = str2double(Pmax);
Qmin = str2double(Qmin);
Qmax = str2double(Qmax);
ax = str2double(ax);
cx = str2double(cx);

timespan = MATDSS.Sim.Meas.at;
setpoint_timespan = MATDSS.Cont.at;

% Getting bus number (check if bus name is given)
if Bus(1) == "'" %string of bus name
    BusNum = MATDSS_StrComp(MATDSS.Sim.DSSCircuit.AllBusNames,Bus(2:end-1));
else %string of bus number
    BusNum = str2double(Bus);
end

%Select the bus in OpenDSS to define the DER OpenDSS load
myBus = MATDSS.Sim.DSSCircuit.Buses(int16(BusNum-1));


% Nphase = length(Nodes);
NodesText = [];
myNodes = Nodes;
while ~isempty(myNodes)
    NodesText = [NodesText,'.',num2str(myNodes(1))];
    myNodes = myNodes(2:end);
end




% check if this DER is a VDER, then don't define it in OpenDSS and there is
% no need for DER dynamics.
if ~strcmpi(DERType,'vder')
% if this DER is an actual DER device, then define it in OpenDSS
% DSSCommand = "New Load.S37a  Bus1=37.1   Phases=1 Conn=Wye   Model=2 kV=2.4   kW=40.0  kvar=20.0"
DSSNewDERCommand = ['New Load.' Name '  Bus1=' myBus.Name NodesText   '  Phases=' num2str(Nphase) ...
    '  Conn=' ConnType  '  Model=1 kV=' num2str(myBus.kVBase*sqrt(3)) '  kW=0  kvar=0 Vminpu=0.1'];% pf=1'];

% Prepare the DER struct to handle DER dynamics
DER(DERInd).x.ODETime = [];
DER(DERInd).x.ODEx = [];
MATDSS.Sim.DSSText.Command = DSSNewDERCommand;
DER(DERInd).Tau = Tau;
end

% Save all values in the DER struct
DER(DERInd).DER_Ind = DERInd;
DER(DERInd).DSSName = Name;
DER(DERInd).BusNum = BusNum;
DER(DERInd).BusName = myBus.Name;
DER(DERInd).Type = DERType;
DER(DERInd).Nodes = Nodes;
DER(DERInd).ConnType = ConnType;
DER(DERInd).NPhases = Nphase;
DER(DERInd).Mode = Mode;
DER(DERInd).Px = Px;
DER(DERInd).Qx = Qx;
DER(DERInd).Pmin = Pmin;
DER(DERInd).Pmax = Pmax;
DER(DERInd).Qmin = Qmin;
DER(DERInd).Qmax = Qmax;
DER(DERInd).ax = ax;
DER(DERInd).cx = cx;
DER(DERInd).MultiPhaseRatio = 1/Nphase; %1/3, 1/2 or 1 ratio for multi-phase elements (check QR definition in the notes).
DER(DERInd).W = [timespan', zeros(length(timespan),1)];
DER(DERInd).Var = [timespan', zeros(length(timespan),1)];
DER(DERInd).i_WVar = 1;
DER(DERInd).WSetpoint = [setpoint_timespan', zeros(length(setpoint_timespan),1)]; % [time, value]
DER(DERInd).VarSetpoint = [setpoint_timespan', zeros(length(setpoint_timespan),1)]; % [time, value]
DER(DERInd).WUnfilteredSetpoint = [zeros(length(setpoint_timespan),1)]; % [value]
DER(DERInd).VarUnvfilteredSetpoint = [zeros(length(setpoint_timespan),1)]; % [value]
DER(DERInd).i_setpoint = 1;
DER(DERInd).UpdateStatus = 1;
DER(DERInd).LLCMoveSetPoint = 0; % This is a bool to tell the higher level controller to move current set-point due to lower-level area constraints violations.


% update the number of DERs in the Sim struct.
MATDSS.Sim.nDER = MATDSS.Sim.nDER + 1; % Since the new element is added, let's update the counter




clearvars -except app MATDSS DER

end

%}