function [P0, Q0, P0Set, Q0Set] = MATDSS_CAP0Q0Calculation(MATDSS,DER,CA,Meast)
% MATDSS_CAP0Q0Calculation(MATDSS,DER,CA,Meast)
% This function computes active power (P0) and reactive power (Q0) as well
% as the respective setpoints (P0Set, Q0Set) for a given control area (CA)
% in a distributed energy resources (DER) system using MATDSS.
%
% The function distinguishes between the highest-level control area and 
% lower-level areas within the feeder and processes the corresponding bus 
% and line measurements to compute the power flow.
%
% Parameters:
%   - MATDSS: The main simulation structure containing system information,
%             control signals, measurements, and simulation settings.
%   - DER:    The structure array containing information on distributed 
%             energy resources in the system.
%   - CA:     The control area structure, representing an area being
%             controlled within the MATDSS simulation.
%   - Meast:  The index or identifier of the current measurement set 
%             being used for simulation purposes.
%
% The function handles both LLC and other control area types, retrieving 
% power measurements from buses and lines, and projects them to the active 
% control area. For interface buses, it computes the voltage and current 
% of the line to determine the power flow at the control area boundary.
%
% Last Update: MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com



% Extract current time and active time steps from the MATDSS simulation
t = MATDSS.Sim.t; % Current time step
at = MATDSS.Sim.at; % Active time

% Get the VDER (Voltage DER) index from the control area
VDERIndex = CA.VDERIndex;

% Check if this is the highest-level control area (no VDER interface)
if VDERIndex == -1 
    % Highest level control area, retrieve the active power setpoint P0Set
    P0Set = MATDSS.ControlSignals.P0Set(t); 
    
    % If the control area type is 'llc', check for reactive power setpoints (Q0Set)
    if strcmpi(CA.Type,'llc') 
        if ~isfield(MATDSS.ControlSignals,'Q0Set') 
            % If Q0Set is not available, compute the sum of reactive power from the profile
            Q0Set = sum(MATDSS.Meas.QProfile(CA.CABus0,Meast));
        else
            % Otherwise, retrieve the reactive power setpoint
            Q0Set = MATDSS.ControlSignals.Q0Set(t); 
        end
    else
        % For other control area types, sum reactive power from the profile
        Q0Set = sum(MATDSS.Meas.QProfile(CA.CABus0,Meast));
    end
    
    % Retrieve active and reactive power measurements at the bus
    P0 = MATDSS.Meas.P0(:,Meast); % Active power measurements
    Q0 = MATDSS.Meas.Q0(:,Meast); % Reactive power measurements

else
    % Lower-level control area: retrieve the interface bus details
    CAInterfaceBuses = MATDSS.Cont.CAInterfaceBuses;

    % Find the matching control area index for the interface buses
    j = find([CAInterfaceBuses{:,5}] == CA.Area);
    
    % Set the active bus in the simulation based on the found interface bus
    MATDSS.Sim.DSSCircuit.SetActiveBus(CAInterfaceBuses{j,2});
    
    % Retrieve the voltage information at the active bus
    MyBus = MATDSS.Sim.DSSCircuit.ActiveBus;
    MyBusV = MyBus.Voltages; % Get the voltage data
    MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end); % Complex representation
    MyBusV = MyBusV.'; % Transpose to row vector
    
    % Retrieve all the line elements in the simulation
    MyLines = MATDSS.Sim.DSSCircuit.Lines;
    MyLinesAllNames = MyLines.AllNames; % Get all line names
    
    % Loop through each line to find the interface line connected to the bus
    for k = 1:length(MyLinesAllNames)
        success = MATDSS.Sim.DSSCircuit.SetActiveElement(['Line.' MyLinesAllNames{k}]); % Set active element to line k
        if ~success
            error('Error in setting active element in Branch current measurements in MATDSS_IPQ.m');
        end

        % Retrieve the current active line and its bus names
        MyLine = MATDSS.Sim.DSSCircuit.ActiveElement;
        MyLineBuses = MyLine.BusNames;

        % Find the location of "dots" in the bus names to separate node info
        iDotBus1 = strfind(MyLineBuses{1},'.');
        if ~isempty(iDotBus1)
            MyBus1 = MyLineBuses{1};
            MyLineBuses(1) = {MyBus1(1:iDotBus1(1)-1)};
        end
        iDotBus2 = strfind(MyLineBuses{2},'.');
        if ~isempty(iDotBus2)
            MyBus2 = MyLineBuses{2};
            MyLineBuses(2) = {MyBus2(1:iDotBus2(1)-1)};
        end

        % Check if the current line is an interface line between the buses
        if sum(strcmpi(MyLineBuses,CAInterfaceBuses{j,2})) && sum(strcmpi(MyLineBuses,CAInterfaceBuses{j,4}))
            % Retrieve the line current and voltage values
            MyLineI = MyLine.Currents;
            MyLineI = MyLineI(1:2:end) + 1i*MyLineI(2:2:end);
            MyBusV = MyLine.Voltages;
            MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end);

            % Determine whether the interface is at bus 1 or bus 2
            if strcmpi(MyLineBuses{1},CAInterfaceBuses{j,2})
                MyLineI = MyLineI(1:length(MyLineI)/2).'; % Get currents from bus 1
                MyBusV = MyBusV(1:length(MyBusV)/2).'; % Get voltages from bus 1
                break;
            else
                MyLineI = MyLineI(length(MyLineI)/2 + 1:end).'; % Get currents from bus 2
                MyBusV = MyBusV(length(MyBusV)/2 + 1:end).'; % Get voltages from bus 2
                break;
            end
        end
    end
    
    % Compute complex power (S) at the control area interface
    MyCAInterfaceS0 = MyBusV.*conj(MyLineI);
    MyCAInterfaceP0 = real(MyCAInterfaceS0); % Active power component
    MyCAInterfaceQ0 = imag(MyCAInterfaceS0); % Reactive power component
    P0 = MyCAInterfaceP0; % Active power output
    Q0 = MyCAInterfaceQ0; % Reactive power output
    
    % For the first time step, set the DER setpoints for power
    if t == 1
        DER(CA.VDERIndex).WSetpoint(:,2) = ones(size(DER(CA.VDERIndex).WSetpoint,1),1).*(-sum(MyCAInterfaceP0));
        DER(CA.VDERIndex).VarSetpoint(:,2) = ones(size(DER(CA.VDERIndex).VarSetpoint,1),1).*(-sum(MyCAInterfaceQ0));
    end

    % Retrieve the power setpoints (P0Set, Q0Set) from the DER
    P0Set = -DER(CA.VDERIndex).WSetpoint(DER(CA.VDERIndex).i_setpoint,2);
    Q0Set = -DER(CA.VDERIndex).VarSetpoint(DER(CA.VDERIndex).i_setpoint,2);
end

end


%% Old function code

%{

function [P0, Q0, P0Set, Q0Set] = MATDSS_CAP0Q0Calculation(MATDSS,DER,CA,Meast)

t = MATDSS.Sim.t;
at = MATDSS.Sim.at;

VDERIndex = CA.VDERIndex;
if VDERIndex == -1 % this is the highest level control area in the feeder
    P0Set = MATDSS.ControlSignals.P0Set(t);
    if strcmpi(CA.Type,'llc')
        if ~isfield(MATDSS.ControlSignals,'Q0Set')
            Q0Set = sum(MATDSS.Meas.QProfile(CA.CABus0,Meast));
        else
            Q0Set = MATDSS.ControlSignals.Q0Set(t);
        end
    else
        Q0Set = sum(MATDSS.Meas.QProfile(CA.CABus0,Meast));
    end
    P0 = MATDSS.Meas.P0(:,Meast);
    Q0 = MATDSS.Meas.Q0(:,Meast);
else
    CAInterfaceBuses = MATDSS.Cont.CAInterfaceBuses;

    j = find([CAInterfaceBuses{:,5}] == CA.Area);
    %                 BusNum = find(strcmpi(AllBusNames,CAInterfaceBuses{j,2}));
    MATDSS.Sim.DSSCircuit.SetActiveBus(CAInterfaceBuses{j,2});
    MyBus = MATDSS.Sim.DSSCircuit.ActiveBus;
    MyBusV = MyBus.Voltages;
    MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end);
    MyBusV = MyBusV.';%.*1e3;

    MyLines = MATDSS.Sim.DSSCircuit.Lines;
    MyLinesAllNames = MyLines.AllNames;

    for k = 1:length(MyLinesAllNames)
        success = MATDSS.Sim.DSSCircuit.SetActiveElement(['Line.' MyLinesAllNames{k}]); % set active element to line k
        if ~success
            error('Error in setting active element in Branch current measurements in MATDSS_IPQ.m');
        end

        MyLine = MATDSS.Sim.DSSCircuit.ActiveElement; % handle of the current line i
        MyLineBuses = MyLine.BusNames;

        iDotBus1 = strfind(MyLineBuses{1},'.'); % find the location of "dots" in the node name to find the busname
        if ~isempty(iDotBus1) % If the bus name we got was only the name without nodes, add the nodes manually and then remove them. This is important so our code can work with all cases.
            MyBus1 = MyLineBuses{1};
            MyLineBuses(1) = {MyBus1(1:iDotBus1(1)-1)};
        end
        iDotBus2 = strfind(MyLineBuses{2},'.'); % find the location of "dots" in the node name to find the busname
        if ~isempty(iDotBus2) % If the bus name we got was only the name without nodes, add the nodes manually and then remove them. This is important so our code can work with all cases.
            MyBus2 = MyLineBuses{2};
            MyLineBuses(2) = {MyBus2(1:iDotBus2(1)-1)};
        end

        if sum(strcmpi(MyLineBuses,CAInterfaceBuses{j,2})) && sum(strcmpi(MyLineBuses,CAInterfaceBuses{j,4}))
            MyLineI = MyLine.Currents;
            MyLineI = MyLineI(1:2:end) + 1i*MyLineI(2:2:end);
            MyBusV = MyLine.Voltages;
            MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end);
            if strcmpi(MyLineBuses{1},CAInterfaceBuses{j,2})
                MyLineI = MyLineI(1:length(MyLineI)/2).';
                MyBusV = MyBusV(1:length(MyBusV)/2).';
                break;

            else
                MyLineI = MyLineI(length(MyLineI)/2 + 1:end).';
                MyBusV = MyBusV(length(MyBusV)/2 + 1:end).';
                break;
            end
        end
    end
    
    
    
    MyCAInterfaceS0 = MyBusV.*conj(MyLineI);
    MyCAInterfaceP0 = real(MyCAInterfaceS0);
    MyCAInterfaceQ0 = imag(MyCAInterfaceS0);
    P0 = MyCAInterfaceP0;
    Q0 = MyCAInterfaceQ0;
    % this has been moved to DER Update
    % DER(CA.VDERIndex).W = [DER(CA.VDERIndex).W, -sum(P0)];%-DER(CA.VDERIndex).W(1)];
    % DER(CA.VDERIndex).Var = [DER(CA.VDERIndex).Var, -sum(Q0)];%-DER(CA.VDERIndex).Var(1)];
    if t == 1
        DER(CA.VDERIndex).WSetpoint(:,2) = ones(size(DER(CA.VDERIndex).WSetpoint,1),1).*(-sum(MyCAInterfaceP0));
        DER(CA.VDERIndex).VarSetpoint(:,2) = ones(size(DER(CA.VDERIndex).VarSetpoint,1),1).*(-sum(MyCAInterfaceQ0));
        % DER(CA.VDERIndex).WSetpoint(DER(CA.VDERIndex).i_setpoint,2) = -sum(P0);% + DER(CA.VDERIndex).W(1);
        % DER(CA.VDERIndex).VarSetpoint(DER(CA.VDERIndex).i_setpoint,2) = -sum(Q0);% + DER(CA.VDERIndex).Var(1);
    end

    P0Set = -DER(CA.VDERIndex).WSetpoint(DER(CA.VDERIndex).i_setpoint,2);% + DER(CA.VDERIndex).W(1);
    Q0Set = -DER(CA.VDERIndex).VarSetpoint(DER(CA.VDERIndex).i_setpoint,2);% + DER(CA.VDERIndex).Var(1);


    

end

end


%}