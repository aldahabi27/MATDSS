function MATDSS = MATDSS_YBusMatrix(MATDSS, Flag_ForceGenerateYBus, Exclude_Source_Load_Modeling_Flag)
% MATDSS_YBusMatrix(MATDSS, Flag_ForceGenerateYBus, Exclude_Source_Load_Modeling_Flag)
% This function obtains the traditional Y-bus matrix for the OpenDSS circuit/file.
% By default, the Y-matrix is generated while excluding load and source modeling,
% to align with the results obtained in voltages and powers from OpenDSS's internal functions.
% This exclusion flag is enabled by default.
%
% Parameters:
%   - MATDSS: The MATDSS application structure containing properties and methods 
%             related to the simulation and Y-bus matrix generation.
%   - Flag_ForceGenerateYBus: A flag to force the generation of the Y-bus matrix 
%                             even if the simulation data is already available.
%   - Exclude_Source_Load_Modeling_Flag: A flag to exclude load and source modeling 
%                                         when obtaining the Y-bus matrix.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Check the number of input arguments and set default values if necessary
if nargin < 3
    Exclude_Source_Load_Modeling_Flag = 1; % Default to excluding source/load modeling
    if nargin < 2
        Flag_ForceGenerateYBus = 0; % Default flag to not force Y-bus generation
    end
end

% Check if simulation data should be used or if Y-bus generation is forced
if ~MATDSS.UseSimDataFlag || Flag_ForceGenerateYBus
    % This flag can be used to get the OpenDSS Ybus matrix (not the format used
    % in textbooks). To do so, set the flag to 0.
    
    AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames; % Get all node names from the DSS circuit
    % To get Y-bus matrix (Traditional Y-bus matrix), execute these commands in
    % OpenDSS and solve the PF (power flow).
    if Exclude_Source_Load_Modeling_Flag
        MATDSS.Sim.DSSText.Command = 'vsource.source.enabled=no'; % Disable the voltage source
        MATDSS.Sim.DSSText.command = 'batchedit load..* enabled=no'; % Disable all loads
        MATDSS.Sim.DSSText.Command = 'solve'; % Solve the circuit to obtain Y-bus matrix
    end

    Nnodes = MATDSS.Sim.DSSCircuit.NumNodes; % This includes all phases considered.
    Ybus = MATDSS.Sim.DSSCircuit.SystemY; % Y-matrix contains only line and shunt admittances
    Ybus = Ybus(1:2:end) + 1i*Ybus(2:2:end); % Convert to complex Y matrix
    MATDSS.Sim.Ybus = reshape(Ybus, Nnodes, []); % Convert to square matrix

    YbusNodesNames = lower(MATDSS.Sim.DSSCircuit.YNodeOrder); % Get ordered node names in lowercase
    A = MATDSS_StrComp(YbusNodesNames, AllNodesNames); % Compare node names for ordering
    % TempAllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames; % this was
    % found to be the same as original AllNodesNames. It does not change,
    % but only Y-bus nodes order changes!
    NewYbus = MATDSS.Sim.Ybus(A, :); % Rearrange Y-bus matrix rows according to node order
    NewYbus = NewYbus(:, A); % Rearrange Y-bus matrix columns according to node order
    MATDSS.Sim.Ybus = NewYbus; % Update the Y-bus matrix with the new order
    % TempAllNodesNames = YbusNodesNames(A);
    % TempAllNodesNamesNew = TempAllNodesNames(A);
    % Once obtained the traditional Y-bus matrix, return the OpenDSS
    % environment to the default operating conditions.
    if Exclude_Source_Load_Modeling_Flag
        MATDSS.Sim.DSSText.Command = 'vsource.source.enabled=yes'; % Re-enable the voltage source
        MATDSS.Sim.DSSText.command = 'batchedit load..* enabled=yes'; % Re-enable all loads
        MATDSS.Sim.DSSText.Command = 'solve'; % Solve again to update circuit conditions
    end

    Nnodes = MATDSS.Sim.DSSCircuit.NumNodes; % Get the number of nodes again
    Ybus = MATDSS.Sim.DSSCircuit.SystemY; % Retrieve the Y-matrix again
    Ybus = Ybus(1:2:end) + 1i*Ybus(2:2:end); % Convert to complex Y matrix again
    MATDSS.Sim.DSSYbus = reshape(Ybus, Nnodes, []); % Convert to square matrix for DSSYbus

elseif MATDSS.UseSimDataFlag
    MATDSS.Sim.Ybus = MATDSS.LoadedSimData.MATDSS.Sim.Ybus; % Load pre-existing Y-bus from simulation data
    MATDSS.Sim.DSSYbus = MATDSS.LoadedSimData.MATDSS.Sim.DSSYbus; % Load pre-existing DSS Y-bus from simulation data
end

% MATDSS.Sim.DSSYbus = MATDSS.Sim.Ybus; % Uncomment this line if needed for further processing

end



%% Old function code

%{

function MATDSS = MATDSS_YBusMatrix(MATDSS, Flag_ForceGenerateYBus, Exclude_Source_Load_Modeling_Flag)
% Obtain traditional Y-bus matrix for OpenDSS circuit/File
%
% By default, we obtain the Y-matrix by excluding load and source modelling
% that is included in OpenDSS to get matching results in Voltages and
% Powers. This is what OpenDSS does even to obtain the powers from its
% internal functions. The flag is enabled by default.
%
%
%
%
if nargin < 3
     Exclude_Source_Load_Modeling_Flag = 1;
    if nargin < 2
        Flag_ForceGenerateYBus = 0;
    end
end

if ~MATDSS.UseSimDataFlag || Flag_ForceGenerateYBus
    % This flag can be used to get the OpenDSS Ybus matrix (not the format used
    % in textbooks). To do so, set the flag to 0.
    

    AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
    % To get Y-bus matrix (Traditional Y-bus matrix), excute these commands in
    % OpenDSS and solve the PF.
    if Exclude_Source_Load_Modeling_Flag
        MATDSS.Sim.DSSText.Command = 'vsource.source.enabled=no';
        MATDSS.Sim.DSSText.command = 'batchedit load..* enabled=no';
        MATDSS.Sim.DSSText.Command = 'solve';
    end

    Nnodes = MATDSS.Sim.DSSCircuit.NumNodes; % this includes all phases considered.
    Ybus = MATDSS.Sim.DSSCircuit.SystemY; % Ymatrix contain only line and shunt admittances
    Ybus = Ybus(1:2:end) + 1i*Ybus(2:2:end); % convert it complex Y matrix and then save it
    MATDSS.Sim.Ybus = reshape(Ybus,Nnodes,[]); % Convert it to square matrix
    
    
    YbusNodesNames = lower(MATDSS.Sim.DSSCircuit.YNodeOrder);
    A = MATDSS_StrComp(YbusNodesNames,AllNodesNames);
    % TempAllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames; % this was
    % found to be the same as original AllNodesNames. It does not change,
    % but only ybus nodes order changes!
    NewYbus = MATDSS.Sim.Ybus(A,:);
    NewYbus = NewYbus(:,A);
    MATDSS.Sim.Ybus = NewYbus;
    % TempAllNodesNames = YbusNodesNames(A);
    % TempAllNodesNamesNew = TempAllNodesNames(A);
    % Once obtained the traditional Y-bus matrix, return the OpenDSS
    % environment to the default operating conditions.
    if Exclude_Source_Load_Modeling_Flag
        MATDSS.Sim.DSSText.Command = 'vsource.source.enabled=yes';
        MATDSS.Sim.DSSText.command = 'batchedit load..* enabled=yes';
        MATDSS.Sim.DSSText.Command = 'solve';
    end


    Nnodes = MATDSS.Sim.DSSCircuit.NumNodes; % this includes all phases considered.
    Ybus = MATDSS.Sim.DSSCircuit.SystemY; % Ymatrix contain only line and shunt admittances
    Ybus = Ybus(1:2:end) + 1i*Ybus(2:2:end); % convert it complex Y matrix and then save it
    MATDSS.Sim.DSSYbus = reshape(Ybus,Nnodes,[]); % Convert it to square matrix


elseif MATDSS.UseSimDataFlag
    MATDSS.Sim.Ybus = MATDSS.LoadedSimData.MATDSS.Sim.Ybus;
    MATDSS.Sim.DSSYbus = MATDSS.LoadedSimData.MATDSS.Sim.DSSYbus;
end



% MATDSS.Sim.DSSYbus = MATDSS.Sim.Ybus;

end

%}