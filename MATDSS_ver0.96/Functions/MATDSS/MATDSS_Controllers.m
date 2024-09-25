function [MATDSS, DER] = MATDSS_Controllers(MATDSS, DER)
% [MATDSS, DER] = MATDSS_Controllers(MATDSS, DER)
% This function is responsible for defining the control areas and assigning
% buses to controllers based on the settings provided in the configuration
% files within the MATDSS application.
%
% Each controller is created with specified parameters, and all relevant
% buses and nodes are assigned to the respective control areas. The function
% ensures that each control area (CA) is correctly defined and associated
% with all buses within it. Additionally, it provides the flexibility to
% handle different controller types and manage both L2C and non-L2C areas.
%
% Parameters:
%   - MATDSS: Structure containing the simulation and measurement data
%             relevant to the electrical distribution system simulation.
%   - DER:    Distributed Energy Resource (DER) data structure.
%
% Last Update: MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com



% Initialization of controller settings and control areas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ____                    _   _       _   _                  ____    _               ____            _
%  |  _ \    ___    _ __   ( ) | |_    | | | |  ___    ___    / ___|  (_)  _ __ ___   |  _ \    __ _  | |_    __ _
%  | | | |  / _ \  | '_ \  |/  | __|   | | | | / __|  / _ \   \___ \  | | | '_ ` _ \  | | | |  / _` | | __|  / _` |
%  | |_| | | (_) | | | | |     | |_    | |_| | \__ \ |  __/    ___) | | | | | | | | | | |_| | | (_| | | |_  | (_| |
%  |____/   \___/  |_| |_|      \__|    \___/  |___/  \___|   |____/  |_| |_| |_| |_| |____/   \__,_|  \__|  \__,_|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~MATDSS.UseSimDataFlag
    % If simulation data is not being used, initialize necessary control area parameters

    CA = [];  % Control Area (CA) structure to store controller settings
    AllBusNodesNames = {};  % List to store bus names extracted from node names
    AllNodesNames = MATDSS.Sim.Meas.AllNodesNames;  % Retrieve all node names from simulation data
    AllBusNames = MATDSS.Sim.Meas.AllBusNames;  % Retrieve all bus names from simulation data

    % Loop through all nodes to extract bus names
    for i = 1:length(AllNodesNames)
        iNodeName = AllNodesNames{i};
        iDot = strfind(iNodeName, '.');  % Find the location of the dot in the node name
        iBusname = iNodeName(1:iDot(1)-1);  % Extract the bus name from the node name
        AllBusNodesNames = [AllBusNodesNames; iBusname];  % Append bus name to list
    end

    % Assign source bus name as the first bus in the list of all bus names
    SourceBusName = AllBusNames{1};

    % Get all line and branch names from the simulation data
    AllLinesNames = MATDSS.Sim.Meas.AllBranchesPhasesNames;
    MyLines = MATDSS.Sim.DSSCircuit.Lines;  % Retrieve DSSCircuit line data
    MyLinesAllNames = MyLines.AllNames;  % Retrieve all line names from DSSCircuit

    % Determine the number of control areas based on unique settings in CASettings table
    NumberOfCA = length(unique(str2double(MATDSS.TableData.CASettings(:, 3))));

    % Define progress bar ranges (MyPer indicates starting point, MyPerEnd is end point)
    MyPer = 0.55;
    MyPerEnd = 0.75;
    MyPerSec = (MyPerEnd - MyPer) / 4;

    % Extracting Control Areas settings from the table (CASettings)
    for i = 1:NumberOfCA
        Controller = struct;  % Initialize a new controller structure

        % Assign controller settings from the table
        [...
            Controller.Area, ...
            Controller.Type, ...
            Controller.alpha, ...
            Controller.rp, ...
            Controller.rbard, ...
            Controller.E, ...
            Controller.vul, ...
            Controller.vll, ...
            Controller.iul, ...
            Controller.arho, ...
            Controller.asigma, ...
            Controller.alambda, ...
            Controller.amu, ...
            Controller.aeta, ...
            Controller.apsi, ...
            Controller.agamma, ...
            Controller.anu, ...
            Controller.azeta, ...
            Controller.crho, ...
            Controller.csigma, ...
            Controller.clambda, ...
            Controller.cmu, ...
            Controller.ceta, ...
            Controller.cpsi, ...
            Controller.cgamma, ...
            Controller.cnu, ...
            Controller.czeta] = deal(MATDSS.TableData.ControllersSettings{i, :});

        ControllerFields = fieldnames(Controller);  % Get field names of controller structure

        % Convert string entries to numeric values for specific fields
        for j = 1:length(ControllerFields)
            if j ~= [2, 3, 9]  % Skip Type, alpha, and iul fields (special cases)
                Controller.(ControllerFields{j}) = str2double(Controller.(ControllerFields{j}));  % Convert to double
            end

            if j == 3 && ~strcmp(Controller.(ControllerFields{j}), 'auto')
                Controller.(ControllerFields{j}) = str2double(Controller.(ControllerFields{j}));  % Convert alpha if not 'auto'
            end
        end

        % Initialize the bus list for this controller
        Controller.Buses = cell(0, 5);

        % Find buses that belong to this controller based on CASettings table
        for j = 1:size(MATDSS.TableData.CASettings, 1)
            if str2double(MATDSS.TableData.CASettings{j, 3}) == i
                Controller.Buses = [Controller.Buses; MATDSS.TableData.CASettings(j, 1), ...
                    str2double(MATDSS.TableData.CASettings(j, 2)), ...
                    str2double(MATDSS.TableData.CASettings(j, 3)), ...
                    MATDSS.TableData.CASettings(j, 4), ...
                    MATDSS.TableData.CASettings(j, 5)];
            end
        end

        % If the area is not L2C, find the interface bus and add it to the list of buses
        if MATDSS_StrComp(Controller.Buses(:, 1), SourceBusName) <= 0
            CABus0 = Controller.Buses(find(strcmp('t', Controller.Buses(:, 4))), 5);
            CASettingsIndex = find(strcmp(MATDSS.TableData.CASettings(:, 1), CABus0));
            Controller.Buses = [MATDSS.TableData.CASettings(CASettingsIndex, 1), ...
                str2double(MATDSS.TableData.CASettings(CASettingsIndex, 2)), ...
                str2double(MATDSS.TableData.CASettings(CASettingsIndex, 3)), ...
                MATDSS.TableData.CASettings(CASettingsIndex, 4), ...
                MATDSS.TableData.CASettings(CASettingsIndex, 5); ...
                Controller.Buses];
        end

        % Find the indices of all buses in this control area
        CABusesIndices = [];
        ControllerBusesIndices = [];
        for j = 2:size(Controller.Buses, 1)
            CABusesIndices = [CABusesIndices; find(strcmp(AllBusNodesNames, Controller.Buses(j, 1)))];
            if Controller.Buses{j, 3} == i
                ControllerBusesIndices = [ControllerBusesIndices; find(strcmp(AllBusNodesNames, Controller.Buses(j, 1)))];
            end
        end

        % Store nodes and bus indices in the controller structure
        Controller.NodesIndices = CABusesIndices;
        Controller.nBuses = size(Controller.Buses, 1);
        Controller.ControllerBusesIndices = ControllerBusesIndices;

        % Append the current controller to the control area list (CA)
        CA = [CA; Controller];

        % Update status (progress bar for simulation - commented out)
        % MATDSSApp_Status(app, MyPer + MyPerSec * i / NumberOfCA);

        clear Controller  % Clear the controller structure for the next iteration
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Get All Interface Buses and Their Corresponding Information %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This section processes the interface buses between control areas.
    % Interface buses are where control areas connect, and we need to get
    % details about these connections (e.g., bus name, control areas involved,
    % phases, etc.). These connections are critical for defining virtual DERs
    % in ABMH functions.

    CAInterfaceBuses = {};
    % Initialize an empty list to store interface buses and their attributes.
    % This variable will hold information such as the control area, buses,
    % number of phases, and connected nodes.

    % Format:
    % 1. Control area number (CA#) where the virtual DER is required.
    % 2. Bus1: The interface bus where the virtual DER will be located.
    % 3. Bus1CA: The control area of Bus1 (same as CA#).
    % 4. Bus2: The bus connected to Bus1.
    % 5. Bus2CA: The control area of Bus2, which is connected to a higher
    %    level control area. So it is the one required to provide the
    %    information about P,Q and their limits as a VDER.
    % 6. nPhases: The number of phases at the interface bus.
    % 7. Nodes: Nodes connected to Bus2 from Bus1.

    CABusesTIndex = find(strcmp(MATDSS.TableData.CASettings(:,4), 't'));
    % Find all interface buses marked with 't' in CASettings.
    CABus0 = {};

    % Loop through all interface buses.
    for i = 1:length(CABusesTIndex)
        % Get the bus (Bus1) from the CASettings table.
        Bus1 = MATDSS.TableData.CASettings{CABusesTIndex(i), 5};

        % Check which control area (CA) Bus1 belongs to.
        for j = 1:NumberOfCA
            if MATDSS_StrComp(CA(j).Buses(:, 1), Bus1) > 0
                Bus1CA = j; % Assign the control area number for Bus1.
                CANumber = j; % Assign the current control area number (CANumber).
                break;
            end
        end

        % Get Bus2 from the table, and find its control area (Bus2CA).
        Bus2 = MATDSS.TableData.CASettings{CABusesTIndex(i), 1};
        Bus2CA = MATDSS.TableData.CASettings{CABusesTIndex(i), 3};

        % Check if Bus2 is the source bus (if L2C control area).
        if strcmpi(Bus2, SourceBusName)
            nPhases = 3; % Assume 3-phase connection with transmission network (can be updated later).
            CABus0Phases = {};
            % Loop through all phases of Bus2.
            for j = 1:nPhases
                CABus0Phases = [CABus0Phases; [Bus2, '.', num2str(j)]];
            end
            % Save the control area and nodes associated with Bus2 in CABus0.
            CABus0 = [CABus0; {Bus2CA, MATDSS_StrComp(AllNodesNames, CABus0Phases)}];
        else
            % If not an L2C control area, search for the line connecting the two interface buses.
            for k = 1:length(MyLinesAllNames)
                % Loop through all lines to find the one connecting Bus1 and Bus2.
                MyLines.Name = MyLinesAllNames{k}; % Assign the line name.
                LineBus1 = MyLines.Bus1; % Get the name of Bus1.
                iDotBus1 = strfind(LineBus1, '.'); % Locate the dot ('.') in the name to separate bus and node.
                if isempty(iDotBus1)
                    LineBus1 = [LineBus1 '.1.2.3']; % Append phases (1.2.3) if missing.
                    iDotBus1 = strfind(LineBus1, '.'); % Locate the new dot positions.
                end
                LineBus1FullName = LineBus1; % Save the full name of Bus1 (with phases).
                LineBus1 = LineBus1(1:iDotBus1(1)-1); % Extract the bus name without phases.

                % Repeat for Bus2.
                LineBus2 = MyLines.Bus2;
                iDotBus2 = strfind(LineBus2, '.');
                if isempty(iDotBus2)
                    LineBus2 = [LineBus2 '.1.2.3'];
                    iDotBus2 = strfind(LineBus2, '.');
                end
                LineBus2FullName = LineBus2;
                LineBus2 = LineBus2(1:iDotBus2(1)-1);

                % Check if this line connects the interface buses we are processing.
                if (strcmp(Bus1, LineBus1) && strcmp(Bus2, LineBus2)) || ...
                        ((strcmp(Bus2, LineBus1) && strcmp(Bus1, LineBus2)))

                    nPhases = MyLines.Phases; % Get the number of phases for this line.

                    % If Bus1 is not the interface bus, switch Bus1 and Bus2.
                    if ~strcmp(Bus1, LineBus1)
                        LineBus1FullName = LineBus2FullName;
                        iDotBus1 = iDotBus2;
                    end

                    CABus0Phases = {};
                    % Retrieve the phases associated with the interface bus.
                    for j = 1:length(iDotBus1)
                        CABus0Phases = [CABus0Phases; [LineBus1FullName(1:iDotBus1(1)), LineBus1FullName(iDotBus1(j)+1)]];
                    end

                    % Save the interface bus information in CABus0.
                    CABus0 = [CABus0; {Bus2CA, MATDSS_StrComp(AllNodesNames, CABus0Phases)}];
                    break; % Stop the loop once the line is found.
                end
            end

            % Generate a list of nodes connected at Bus1.
            bus1nodes = '[';
            for j = 1:length(iDotBus1)
                bus1nodes = [bus1nodes, LineBus1FullName(iDotBus1(j)+1), ','];
            end
            bus1nodes = [bus1nodes(1:end-1), ']']; % Remove the last comma.

            % Add the control area, buses, and node information to CAInterfaceBuses.
            CAInterfaceBuses = [CAInterfaceBuses; {CANumber, Bus1, Bus1CA, Bus2, str2num(Bus2CA), nPhases, bus1nodes}];
        end

        % Update the status bar (commented out in this version).
        % MATDSSApp_Status(app, MyPer + MyPerSec*(1 + i / length(CABusesTIndex)));
    end % End of loop over all interface buses.

    % Save the control areas (CA) and interface buses (CAInterfaceBuses) in the MATDSS structure.
    MATDSS.Cont.CA = CA;
    MATDSS.Cont.CAInterfaceBuses = CAInterfaceBuses;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                Defining Virtual DERs                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop over all the interface buses and define Virtual DERs (VDERs)
    for i = 1:size(CAInterfaceBuses, 1)

        % Get the table of VDER settings
        VDERsTable = MATDSS.TableData.VDERSettings;

        % Find the corresponding VDER settings based on the control area (CA)
        VDERIndexTable = find(CAInterfaceBuses{i,5} == str2double(VDERsTable(:,1)));

        % Extract specific VDER information for the current control area interface bus
        VDERInfo = VDERsTable(VDERIndexTable,[3,4,5,6,7,8,9]); % Tau, Connection type, DER Mode, P, Q, ax, cx

        % Create a new Virtual DER (VDER) in the MATDSS structure
        [MATDSS, DER] = MATDSS_DERNew(MATDSS,DER,...
            num2str(size(DER,2)+1),...                                    % DER Index (incrementing number)
            ['VDER_CA' num2str(CAInterfaceBuses{i,1}) '_For_CA' num2str(CAInterfaceBuses{i,5})],...  % DER Name with control area info
            strcat("'", CAInterfaceBuses{i,2}, "'"),...                   % DER Bus Location (Bus1 in CA)
            VDERInfo(1),...                                               % Tau (time constant)
            "VDER",...                                                    % DER Type --> Virtual DER (VDER)
            CAInterfaceBuses{i,7},...                                     % Nodes (nodes connected to the DER)
            VDERInfo(2),...                                               % Type of connection (Wye or Delta)
            num2str(CAInterfaceBuses{i,6}),...                            % Number of phases
            VDERInfo(3),...                                               % DER Mode (e.g., Power factor control)
            VDERInfo(4),...                                               % Active power (P)
            VDERInfo(5),...                                               % Reactive power (Q)
            '-1e99','1e99','-1e99','1e99',...                             % Limits for Pmin, Pmax, Qmin, Qmax
            VDERInfo(6),VDERInfo(7));                                     % Coefficients ax and cx

        % MATDSSApp_Status(app, MyPer + MyPerSec*(2 + i / size(CAInterfaceBuses,1)));
    end

    % Now, we are ready to re-define our control areas (CA) as isolated networks.
    % We map the location of DERs and VDERs to the bus numbers inside each control area.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%         GET CA Properties and Information         %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop through all control areas to define CA properties
    for i = 1:NumberOfCA

        % Initialize variables for each control area
        CAnDER = 0;            % Number of DERs in the CA
        CADERIndex = [];       % Indexes of DERs in this CA
        DERBuses = [DER(:).BusNum]';  % List of buses with DERs
        Pmax = 0; Pmin = 0;    % Initialize P limits for the control area
        Qmax = 0; Qmin = 0;    % Initialize Q limits for the control area

        % Loop through all DERs and check if they belong to this CA
        for j = 1:length(DERBuses)
            IsDERInCA = find(DERBuses(j) == [MATDSS.Cont.CA(i).Buses{:,2}]');
            if ~isempty(IsDERInCA)  % If DER is in this CA
                if MATDSS.Cont.CA(i).Buses{IsDERInCA,3} == i  % Check if it matches the CA number
                    CAnDER = CAnDER + 1;  % Increment the number of DERs
                    CADERIndex = [CADERIndex; j];  % Add DER index to the list
                    DER(j).CAIndex = i;  % Assign CA index to the DER

                    % Update the Pmin, Pmax, Qmin, Qmax values for the control area
                    Pmin = Pmin - DER(j).Pmax;
                    Pmax = Pmax - DER(j).Pmin;
                    Qmin = Qmin - DER(j).Qmax;
                    Qmax = Qmax - DER(j).Qmin;
                end
            end
        end

        % Store the information back into MATDSS for the control area
        MATDSS.Cont.CA(i).nDER = CAnDER;  % Number of DERs in the control area
        MATDSS.Cont.CA(i).DERIndex = CADERIndex;  % Store DER indexes for this CA
        area_index = find(i == str2double([CABus0(:,1)]));  % Find the corresponding CA index in CABus0
        MATDSS.Cont.CA(i).CABus0 = CABus0{area_index,2};  % Get the corresponding Bus0 information
        MATDSS.Cont.CA(i).NodesIndices = [MATDSS.Cont.CA(i).CABus0; MATDSS.Cont.CA(i).NodesIndices];  % Update node indices including Bus0
        MATDSS.Cont.CA(i).nNodes = length(MATDSS.Cont.CA(i).NodesIndices);  % Total number of nodes in the CA
        MATDSS.Cont.CA(i).nBuses = size(MATDSS.Cont.CA(i).Buses,1);  % Total number of buses in the CA


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%                    YBus of CA                     %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use Y-reduced method to shrink the Y-bus matrix of the entire network to
        % the Y-bus of this particular control area (CA)

        FullYbus = MATDSS.Sim.Ybus;  % Full Y-bus matrix of the entire network

        % Extract relevant Y-bus portions for the control area
        ControllerYbus = [FullYbus(:,MATDSS.Cont.CA(i).NodesIndices),FullYbus(:,setdiff(1:length(AllBusNodesNames),MATDSS.Cont.CA(i).NodesIndices))];
        ControllerYbus = [ControllerYbus(MATDSS.Cont.CA(i).NodesIndices,:);ControllerYbus(setdiff(1:length(AllBusNodesNames),MATDSS.Cont.CA(i).NodesIndices),:)];

        % Split Y-bus into blocks: Ynn, Yrr, Ynr, Yrn
        nphase = length(AllBusNodesNames);  % Number of phases in the entire network
        nCAphase = length(MATDSS.Cont.CA(i).NodesIndices);  % Number of phases in the control area
        Ynn = ControllerYbus(1:nCAphase,1:nCAphase);  % Ynn block (CA nodes)
        Yrr = ControllerYbus(nCAphase+1:nphase,nCAphase+1:nphase);  % Yrr block (other network nodes)
        Ynr = ControllerYbus(1:nCAphase,nCAphase+1:nphase);  % Ynr block (CA to network)
        Yrn = ControllerYbus(nCAphase+1:nphase,1:nCAphase);  % Yrn block (network to CA)

        % Compute the reduced Y-bus for the control area
        ControllerYred = Ynr/(Yrr);  % Inversion step for Y-bus reduction
        ControllerYred = Ynn - ControllerYred * Yrn;  % Reduced Y-bus for the control area
        MATDSS.Cont.CA(i).Ybus = ControllerYred;  % Store the reduced Y-bus in the MATDSS structure



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%      Find P0, Q0, P0Set, Q0Set of CA and VDER     %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % If this is the highest control area (L2C), no VDER is needed.
        if strcmp(MATDSS.Cont.CA(i).Buses(1,1), SourceBusName)
            % Set active power (P0Set) and reactive power (Q0Set) to zeros for the highest CA.
            MATDSS.Cont.CA(i).P0Set = zeros(1, length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).Q0Set = zeros(1, length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).TimeIndex = 1; % Initialize time index

            % Remove bus indices corresponding to CA Bus0
            temp = MATDSS.Cont.CA(i).CABus0;
            while ~isempty(temp)
                MATDSS.Cont.CA(i).ControllerBusesIndices(MATDSS.Cont.CA(i).ControllerBusesIndices == temp(1)) = [];
                temp(1) = [];
            end

            MATDSS.Cont.CA(i).VDERIndex = -1; % No VDER, get P0 and Q0 from control signals

        else
            % For non-highest CA, find the corresponding VDER index to retrieve setpoints
            for j = 1:size(CAInterfaceBuses, 1)
                if CAInterfaceBuses{j,5} == MATDSS.Cont.CA(i).Area
                    nNotVDER = MATDSS.Sim.nDER - size(CAInterfaceBuses, 1);
                    MATDSS.Cont.CA(i).VDERIndex = nNotVDER + j; % Assign VDER index

                    if MATDSS.Cont.CA(i).VDERIndex > 0
                        % Activate the bus for the VDER interface
                        MATDSS.Sim.DSSCircuit.SetActiveBus(CAInterfaceBuses{j,2});
                        MyBus = MATDSS.Sim.DSSCircuit.ActiveBus;
                        MyBusV = MyBus.Voltages;
                        MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end); % Voltage as complex number
                        MyBusV = MyBusV.';

                        % Get line information to retrieve current measurements
                        MyLines = MATDSS.Sim.DSSCircuit.Lines;
                        MyLinesAllNames = MyLines.AllNames;

                        % Loop through lines to find the corresponding one
                        for k = 1:length(MyLinesAllNames)
                            success = MATDSS.Sim.DSSCircuit.SetActiveElement(['Line.' MyLinesAllNames{k}]);
                            if ~success
                                error('Error in setting active element in Branch current measurements in MATDSS_Controllers.m');
                            end

                            MyLine = MATDSS.Sim.DSSCircuit.ActiveElement;
                            MyLineBuses = MyLine.BusNames;

                            % Ensure bus names are processed correctly by removing node identifiers
                            iDotBus1 = strfind(MyLineBuses{1}, '.');
                            if ~isempty(iDotBus1)
                                MyBus1 = MyLineBuses{1};
                                MyLineBuses(1) = {MyBus1(1:iDotBus1(1)-1)};
                            end
                            iDotBus2 = strfind(MyLineBuses{2}, '.');
                            if ~isempty(iDotBus2)
                                MyBus2 = MyLineBuses{2};
                                MyLineBuses(2) = {MyBus2(1:iDotBus2(1)-1)};
                            end

                            % Check if this line is between the VDER and CA interface
                            if sum(strcmpi(MyLineBuses, CAInterfaceBuses{j,2})) && sum(strcmpi(MyLineBuses, CAInterfaceBuses{j,4}))
                                % Get the currents and voltages from the line element
                                MyLineI = MyLine.Currents;
                                MyLineI = MyLineI(1:2:end) + 1i*MyLineI(2:2:end); % Current as complex numbers
                                MyBusV = MyLine.Voltages;
                                MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end); % Voltage as complex numbers

                                % Select the appropriate half of the data based on bus alignment
                                if strcmpi(MyLineBuses{1}, CAInterfaceBuses{j,2})
                                    MyLineI = MyLineI(1:length(MyLineI)/2).';
                                    MyBusV = MyBusV(1:length(MyBusV)/2).';
                                else
                                    MyLineI = MyLineI(length(MyLineI)/2 + 1:end).';
                                    MyBusV = MyBusV(length(MyBusV)/2 + 1:end).';
                                end
                                break;
                            end
                        end

                        % Calculate interface powers (S = VI*)
                        MyCAInterfaceS0 = MyBusV .* conj(MyLineI);
                        MyCAInterfaceP0 = real(MyCAInterfaceS0); % Active power
                        MyCAInterfaceQ0 = imag(MyCAInterfaceS0); % Reactive power

                        % Assign CA power values
                        MATDSS.Cont.CA(i).MyCAInterfaceP0 = MyCAInterfaceP0;
                        MATDSS.Cont.CA(i).MyCAInterfaceQ0 = MyCAInterfaceQ0;

                        % Set DER power setpoints
                        i_WVar = DER(MATDSS.Cont.CA(i).VDERIndex).i_WVar;
                        DER(MATDSS.Cont.CA(i).VDERIndex).W(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).W,1),1) .* (-sum(MyCAInterfaceP0));
                        DER(MATDSS.Cont.CA(i).VDERIndex).Var(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).Var,1),1) .* (-sum(MyCAInterfaceQ0));
                        DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint,1),1) .* (-sum(MyCAInterfaceP0));
                        DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint,1),1) .* (-sum(MyCAInterfaceQ0));

                        % Set CA power setpoints
                        MATDSS.Cont.CA(i).P0Set = ones(1, length(MATDSS.Cont.at)) .* sum(MyCAInterfaceP0);
                        MATDSS.Cont.CA(i).Q0Set = ones(1, length(MATDSS.Cont.at)) .* sum(MyCAInterfaceQ0);

                        % Adjust DER power limits based on the interface power
                        DER(MATDSS.Cont.CA(i).VDERIndex).Pmax = -sum(MyCAInterfaceP0) + Pmax;
                        DER(MATDSS.Cont.CA(i).VDERIndex).Pmin = -sum(MyCAInterfaceP0) + Pmin;
                        DER(MATDSS.Cont.CA(i).VDERIndex).Qmax = -sum(MyCAInterfaceQ0) + Qmax;
                        DER(MATDSS.Cont.CA(i).VDERIndex).Qmin = -sum(MyCAInterfaceQ0) + Qmin;
                        MATDSS.Cont.CA(i).TimeIndex = 1;
                    end
                end
            end
        end

        % Find the bus indices in the measurement set that correspond to the CA
        [MATDSS.Cont.CA(i).k_v, MATDSS.Cont.CA(i).k_vIndex, ~] = intersect(MATDSS.Meas.k_v, MATDSS.Cont.CA(i).ControllerBusesIndices);


        % Part 5: Processing controller variables, dual variables, and error handles for Control Areas (CA)

        % Find k_i: This part extracts measurement indices (k_i) related to currents in the system
        k_i = []; % Initialize empty array for k_i values
        if MATDSS.Meas.k_i > 0 % Check if there are current measurements available
            k_iLines = AllLinesNames(MATDSS.Meas.k_i); % Get the names of the lines corresponding to current measurements
            for j = 1:length(k_iLines) % Loop through each line with a measurement
                iDotline = strfind(k_iLines{j},'.'); % Find the dot separating bus/phase info in line name
                LineName = k_iLines{j}; % Get the full line name
                MyLines.Name = LineName(1:iDotline(1)-1); % Extract the line name (without phase info)
                LineBus1 = MyLines.Bus1; % Retrieve the name of Bus1 (connected to the line)
                LineBus2 = MyLines.Bus2; % Retrieve the name of Bus2 (other end of the line)

                % Remove phase information from bus names (keep only base name)
                iDotBus1 = strfind(LineBus1,'.');
                iDotBus2 = strfind(LineBus2,'.');
                if ~isempty(iDotBus1)
                    LineBus1 = LineBus1(1:iDotBus1(1)-1); % Remove phase info from Bus1
                end
                if ~isempty(iDotBus2)
                    LineBus2 = LineBus2(1:iDotBus2(1)-1); % Remove phase info from Bus2
                end

                % Check if both buses are part of the current control area (CA)
                iBus1 = find(strcmp(MATDSS.Cont.CA(i).Buses(:,1),LineBus1)); % Find if Bus1 belongs to the CA
                iBus2 = find(strcmp(MATDSS.Cont.CA(i).Buses(:,1),LineBus2)); % Find if Bus2 belongs to the CA

                % If both buses belong to the CA, record the current measurement index
                if ~isempty(iBus1) && ~isempty(iBus2)
                    k_i = [k_i; MATDSS.Meas.k_i(j)]; % Append current measurement index
                end
            end
        end
        MATDSS.Cont.CA(i).k_i = k_i; % Store k_i values in the CA
        [~,MATDSS.Cont.CA(i).k_iIndex,~] = intersect(MATDSS.Meas.k_i,MATDSS.Cont.CA(i).k_i); % Find index of k_i in overall measurement list

        % Compute system matrices (A, B, M, H) using the MATDSS_ABMH function, tailored for the current CA
        [A, B, M, H, Mv, MvIndex, Mi, MiBuses] = MATDSS_ABMH(MATDSS, DER, MATDSS.Cont.CA(i));
        MATDSS.Cont.CA(i).ABMH.A = A; % Store matrix A in the control area
        MATDSS.Cont.CA(i).ABMH.B = B; % Store matrix B in the control area
        MATDSS.Cont.CA(i).ABMH.M = M; % Store matrix M in the control area
        MATDSS.Cont.CA(i).ABMH.H = H; % Store matrix H in the control area

        % Initialize dual variables for the controller for the current CA
        MATDSS.Cont.CA(i).Duals.rho = zeros(size(M,1), length(MATDSS.Cont.at)); % Initialize rho
        MATDSS.Cont.CA(i).Duals.sigma = zeros(size(H,1), length(MATDSS.Cont.at)); % Initialize sigma
        MATDSS.Cont.CA(i).Duals.lambda = zeros(1, length(MATDSS.Cont.at)); % Initialize lambda
        MATDSS.Cont.CA(i).Duals.mu = zeros(1, length(MATDSS.Cont.at)); % Initialize mu
        MATDSS.Cont.CA(i).Duals.eta = zeros(1, length(MATDSS.Cont.at)); % Initialize eta
        MATDSS.Cont.CA(i).Duals.psi = zeros(1, length(MATDSS.Cont.at)); % Initialize psi
        MATDSS.Cont.CA(i).Duals.gamma = zeros(max(size(Mv)), length(MATDSS.Cont.at)); % Initialize gamma (voltage control)
        MATDSS.Cont.CA(i).Duals.nu = zeros(max(size(Mv)), length(MATDSS.Cont.at)); % Initialize nu (voltage control)
        MATDSS.Cont.CA(i).Duals.zeta = zeros(max(size(Mi)), length(MATDSS.Cont.at)); % Initialize zeta (current control)

        % Initialize error handles for the dual variables in the controller
        MATDSS.Cont.CA(i).Duals.Errors.rho = zeros(size(M,1), length(MATDSS.Cont.at)); % Error handle for rho
        MATDSS.Cont.CA(i).Duals.Errors.sigma = zeros(size(H,1), length(MATDSS.Cont.at)); % Error handle for sigma
        MATDSS.Cont.CA(i).Duals.Errors.lambda = zeros(1, length(MATDSS.Cont.at)); % Error handle for lambda
        MATDSS.Cont.CA(i).Duals.Errors.mu = zeros(1, length(MATDSS.Cont.at)); % Error handle for mu
        MATDSS.Cont.CA(i).Duals.Errors.eta = zeros(1, length(MATDSS.Cont.at)); % Error handle for eta
        MATDSS.Cont.CA(i).Duals.Errors.psi = zeros(1, length(MATDSS.Cont.at)); % Error handle for psi
        MATDSS.Cont.CA(i).Duals.Errors.gamma = zeros(max(size(Mv)), length(MATDSS.Cont.at)); % Error handle for gamma
        MATDSS.Cont.CA(i).Duals.Errors.nu = zeros(max(size(Mv)), length(MATDSS.Cont.at)); % Error handle for nu
        MATDSS.Cont.CA(i).Duals.Errors.zeta = zeros(max(size(Mi)), length(MATDSS.Cont.at)); % Error handle for zeta

        % Initialize error handles for the measurement-related variables in the controller
        MATDSS.Cont.CA(i).Duals.ErrorsY.rho = zeros(size(M,1), length(MATDSS.Cont.at)); % Measurement error handle for rho
        MATDSS.Cont.CA(i).Duals.ErrorsY.sigma = zeros(size(H,1), length(MATDSS.Cont.at)); % Measurement error handle for sigma
        MATDSS.Cont.CA(i).Duals.ErrorsY.lambda = zeros(1, length(MATDSS.Cont.at)); % Measurement error handle for lambda
        MATDSS.Cont.CA(i).Duals.ErrorsY.mu = zeros(1, length(MATDSS.Cont.at)); % Measurement error handle for mu
        MATDSS.Cont.CA(i).Duals.ErrorsY.eta = zeros(1, length(MATDSS.Cont.at)); % Measurement error handle for eta
        MATDSS.Cont.CA(i).Duals.ErrorsY.psi = zeros(1, length(MATDSS.Cont.at)); % Measurement error handle for psi
        MATDSS.Cont.CA(i).Duals.ErrorsY.gamma = zeros(max(size(Mv)), length(MATDSS.Cont.at)); % Measurement error handle for gamma
        MATDSS.Cont.CA(i).Duals.ErrorsY.nu = zeros(max(size(Mv)), length(MATDSS.Cont.at)); % Measurement error handle for nu
        MATDSS.Cont.CA(i).Duals.ErrorsY.zeta = zeros(max(size(Mi)), length(MATDSS.Cont.at)); % Measurement error handle for zeta

        % Set alpha gain values for the control area, either 'auto' or user-defined
        if strcmpi(MATDSS.Cont.CA(i).alpha,'auto')
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.Gain .* MATDSS.Time.Cont.TimeStep; % Auto alpha based on controller gain and time step
        else
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.CA(i).alpha .* MATDSS.Cont.Gain .* MATDSS.Time.Cont.TimeStep; % User-defined alpha
        end

        % Initialize the P0 and Q0 setpoints (real and reactive power) for the control area
        MATDSS.Cont.CA(i).P0 = zeros(length(MATDSS.Cont.CA(i).CABus0), length(MATDSS.Cont.at)); % Initialize P0
        MATDSS.Cont.CA(i).Q0 = zeros(length(MATDSS.Cont.CA(i).CABus0), length(MATDSS.Cont.at)); % Initialize Q0

        % Status update in the controller (optional for logging purposes)
        % MATDSSApp_Status(app, MyPer + MyPerSec*(3 + i / NumberOfCA)); % Update the status based on CA index

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  _   _                  ____    _               ____            _            %
    % | | | |  ___    ___    / ___|  (_)  _ __ ___   |  _ \    __ _  | |_    __ _  %
    % | | | | / __|  / _ \   \___ \  | | | '_ ` _ \  | | | |  / _` | | __|  / _` | %
    % | |_| | \__ \ |  __/    ___) | | | | | | | | | | |_| | | (_| | | |_  | (_| | %
    %  \___/  |___/  \___|   |____/  |_| |_| |_| |_| |____/   \__,_|  \__|  \__,_| %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    % If not in simulation mode, copy the controller data from previously loaded simulation data
    MATDSS.Cont.CA = MATDSS.LoadedSimData.MATDSS.Cont.CA; % Load controller data
    MATDSS.Cont.CAInterfaceBuses = MATDSS.LoadedSimData.MATDSS.Cont.CAInterfaceBuses; % Load interface bus data

    % Get all bus names from the measurements
    AllBusNames = MATDSS.Sim.Meas.AllBusNames;
    SourceBusName = AllBusNames{1}; % Set the first bus name as the source bus
    CAInterfaceBuses = MATDSS.Cont.CAInterfaceBuses; % Get the interface buses

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                Defining Virtual DERs                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now we are ready to define our Virtual Distributed Energy Resources (VDERs)!
    for i = 1:size(CAInterfaceBuses,1) % Loop through each interface bus

        % Get the settings for VDERs from the table
        VDERsTable = MATDSS.TableData.VDERSettings;
        VDERIndexTable = find(CAInterfaceBuses{i,5} == str2double(VDERsTable(:,1))); % Find matching VDER index
        VDERInfo = VDERsTable(VDERIndexTable, [3, 4, 5, 6, 7, 8, 9]); % Extract relevant VDER information

        % Create new VDERs based on extracted settings
        [MATDSS, DER] = MATDSS_DERNew(MATDSS, DER, ...
            num2str(size(DER, 2) + 1), ... % DER Index (new DER number)
            ['VDER_CA' num2str(CAInterfaceBuses{i,1}) '_For_CA' num2str(CAInterfaceBuses{i,5})], ... % DER Name
            strcat("'", CAInterfaceBuses{i,2}, "'"), ... % DER Bus Location
            VDERInfo(1), ... % Tau (time constant)
            "VDER", ... % DER Type (specifically VDER)
            CAInterfaceBuses{i,7}, ... % Nodes for the VDER
            VDERInfo(2), ... % Type of connection (e.g., Wye)
            num2str(CAInterfaceBuses{i,6}), ... % Number of phases
            VDERInfo(3), ... % DER Mode
            VDERInfo(4), ... % P(x)
            VDERInfo(5), ... % Q(x)
            '-1e99', '1e99', '-1e99', '1e99', ... % Pmin, Pmax, Qmin, Qmax (set very wide limits)
            VDERInfo(6), VDERInfo(7)); % Additional parameters for ax and cx
    end

    % Get the number of different CA interfaces
    NumberOfCA = length(unique(str2double(MATDSS.TableData.CASettings(:,3))));
    DERBuses = [DER(:).BusNum]'; % Collect all DER bus numbers
    for i = 1:NumberOfCA % Loop through each CA interface
        CAnDER = 0; % Initialize count of DERs in the CA
        CADERIndex = []; % Initialize index array for corresponding DERs
        % Initialize limits for DERs in this CA
        Pmax = 0;
        Pmin = 0;
        Qmax = 0;
        Qmin = 0;

        % Loop through all DERs to check if they are within the current CA interface
        for j = 1:length(DERBuses)
            IsDERInCA = find(DERBuses(j) == [MATDSS.Cont.CA(i).Buses{:,2}]'); % Check if the DER bus matches CA bus
            if ~isempty(IsDERInCA) % If DER is found in CA
                if MATDSS.Cont.CA(i).Buses{IsDERInCA,3} == i % Ensure the bus matches the current CA index
                    CAnDER = CAnDER + 1; % Increment DER count
                    CADERIndex = [CADERIndex; j]; % Store the index of this DER
                    DER(j).CAIndex = i; % Assign the CA index to the DER
                    % Update the limits based on the DER's specifications
                    Pmin = Pmin - DER(j).Pmax;
                    Pmax = Pmax - DER(j).Pmin;
                    Qmin = Qmin - DER(j).Qmax;
                    Qmax = Qmax - DER(j).Qmin;
                end
            end
        end

        % Check if the first bus matches the source bus
        if strcmp(MATDSS.Cont.CA(i).Buses(1,1), SourceBusName)
            % Initialize power set points for the source bus
            MATDSS.Cont.CA(i).P0Set = zeros(1,length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).Q0Set = zeros(1,length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).TimeIndex = 1; % Initialize time index
        else
            % If not the source bus, get current interface power and reactive power
            MyCAInterfaceP0 = MATDSS.Cont.CA(i).MyCAInterfaceP0;
            MyCAInterfaceQ0 = MATDSS.Cont.CA(i).MyCAInterfaceQ0;

            % Initialize the power set points for the interface bus
            MATDSS.Cont.CA(i).P0Set = ones(1,length(MATDSS.Cont.at)) .* (sum(MyCAInterfaceP0));
            MATDSS.Cont.CA(i).Q0Set = ones(1,length(MATDSS.Cont.at)) .* (sum(MyCAInterfaceQ0));

            % Update the DER settings based on the current interface power and reactive power
            DER(MATDSS.Cont.CA(i).VDERIndex).W(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).W, 1), 1) .* (-sum(MyCAInterfaceP0));
            DER(MATDSS.Cont.CA(i).VDERIndex).Var(:,2)  = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).Var, 1), 1) .* (-sum(MyCAInterfaceQ0));
            DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint, 1), 1) .* (-sum(MyCAInterfaceP0));
            DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint, 1), 1) .* (-sum(MyCAInterfaceQ0));
            % Update the maximum and minimum limits of the DER based on the current CA settings
            DER(MATDSS.Cont.CA(i).VDERIndex).Pmax = -sum(MyCAInterfaceP0) + Pmax;
            DER(MATDSS.Cont.CA(i).VDERIndex).Pmin = -sum(MyCAInterfaceP0) + Pmin;
            DER(MATDSS.Cont.CA(i).VDERIndex).Qmax = -sum(MyCAInterfaceQ0) + Qmax;
            DER(MATDSS.Cont.CA(i).VDERIndex).Qmin = -sum(MyCAInterfaceQ0) + Qmin;
        end

        % Initialize dual variables for control area i
        MATDSS.Cont.CA(i).Duals.rho = zeros(size(MATDSS.Cont.CA(i).Duals.rho,1), length(MATDSS.Cont.at)); % Initialize rho with zeros
        MATDSS.Cont.CA(i).Duals.sigma = zeros(size(MATDSS.Cont.CA(i).Duals.sigma,1), length(MATDSS.Cont.at)); % Initialize sigma with zeros
        MATDSS.Cont.CA(i).Duals.lambda = zeros(1, length(MATDSS.Cont.at)); % Initialize lambda with zeros
        MATDSS.Cont.CA(i).Duals.mu = zeros(1, length(MATDSS.Cont.at)); % Initialize mu with zeros
        MATDSS.Cont.CA(i).Duals.eta = zeros(1, length(MATDSS.Cont.at)); % Initialize eta with zeros
        MATDSS.Cont.CA(i).Duals.psi = zeros(1, length(MATDSS.Cont.at)); % Initialize psi with zeros

        % Uncommented block for dual variables gamma and nu initialization
        MATDSS.Cont.CA(i).Duals.gamma = zeros(max(size(MATDSS.Cont.CA(i).Duals.gamma,1)), length(MATDSS.Cont.at)); % Initialize gamma with zeros
        MATDSS.Cont.CA(i).Duals.nu = zeros(max(size(MATDSS.Cont.CA(i).Duals.nu,1)), length(MATDSS.Cont.at)); % Initialize nu with zeros

        % Initialize zeta dual variable with zeros
        MATDSS.Cont.CA(i).Duals.zeta = zeros(max(size(MATDSS.Cont.CA(i).Duals.zeta,1)), length(MATDSS.Cont.at)); % Initialize zeta with zeros

        % Initialize primary control variables P0 and Q0 for control area i
        MATDSS.Cont.CA(i).P0 = zeros(length(MATDSS.Cont.CA(i).CABus0), length(MATDSS.Cont.at)); % Initialize P0 with zeros
        MATDSS.Cont.CA(i).Q0 = zeros(length(MATDSS.Cont.CA(i).CABus0), length(MATDSS.Cont.at)); % Initialize Q0 with zeros

        % Define error handles for dual variables
        MATDSS.Cont.CA(i).Duals.Errors.rho = zeros(size(MATDSS.Cont.CA(i).Duals.rho,1), length(MATDSS.Cont.at)); % Initialize error handle for rho
        MATDSS.Cont.CA(i).Duals.Errors.sigma = zeros(size(MATDSS.Cont.CA(i).Duals.sigma,1), length(MATDSS.Cont.at)); % Initialize error handle for sigma
        MATDSS.Cont.CA(i).Duals.Errors.lambda = zeros(1, length(MATDSS.Cont.at)); % Initialize error handle for lambda
        MATDSS.Cont.CA(i).Duals.Errors.mu = zeros(1, length(MATDSS.Cont.at)); % Initialize error handle for mu
        MATDSS.Cont.CA(i).Duals.Errors.eta = zeros(1, length(MATDSS.Cont.at)); % Initialize error handle for eta
        MATDSS.Cont.CA(i).Duals.Errors.psi = zeros(1, length(MATDSS.Cont.at)); % Initialize error handle for psi
        MATDSS.Cont.CA(i).Duals.Errors.gamma = zeros(max(size(MATDSS.Cont.CA(i).Duals.gamma,1)), length(MATDSS.Cont.at)); % Initialize error handle for gamma
        MATDSS.Cont.CA(i).Duals.Errors.nu = zeros(max(size(MATDSS.Cont.CA(i).Duals.nu,1)), length(MATDSS.Cont.at)); % Initialize error handle for nu
        MATDSS.Cont.CA(i).Duals.Errors.zeta = zeros(max(size(MATDSS.Cont.CA(i).Duals.zeta,1)), length(MATDSS.Cont.at)); % Initialize error handle for zeta

        % Define error handles for Y measurements in the D controller
        MATDSS.Cont.CA(i).Duals.ErrorsY.rho = zeros(size(MATDSS.Cont.CA(i).Duals.rho,1), length(MATDSS.Cont.at)); % Initialize Y measurement error handle for rho
        MATDSS.Cont.CA(i).Duals.ErrorsY.sigma = zeros(size(MATDSS.Cont.CA(i).Duals.sigma,1), length(MATDSS.Cont.at)); % Initialize Y measurement error handle for sigma
        MATDSS.Cont.CA(i).Duals.ErrorsY.lambda = zeros(1, length(MATDSS.Cont.at)); % Initialize Y measurement error handle for lambda
        MATDSS.Cont.CA(i).Duals.ErrorsY.mu = zeros(1, length(MATDSS.Cont.at)); % Initialize Y measurement error handle for mu
        MATDSS.Cont.CA(i).Duals.ErrorsY.eta = zeros(1, length(MATDSS.Cont.at)); % Initialize Y measurement error handle for eta
        MATDSS.Cont.CA(i).Duals.ErrorsY.psi = zeros(1, length(MATDSS.Cont.at)); % Initialize Y measurement error handle for psi
        MATDSS.Cont.CA(i).Duals.ErrorsY.gamma = zeros(max(size(MATDSS.Cont.CA(i).Duals.gamma,1)), length(MATDSS.Cont.at)); % Initialize Y measurement error handle for gamma
        MATDSS.Cont.CA(i).Duals.ErrorsY.nu = zeros(max(size(MATDSS.Cont.CA(i).Duals.nu,1)), length(MATDSS.Cont.at)); % Initialize Y measurement error handle for nu
        MATDSS.Cont.CA(i).Duals.ErrorsY.zeta = zeros(max(size(MATDSS.Cont.CA(i).Duals.zeta,1)), length(MATDSS.Cont.at)); % Initialize Y measurement error handle for zeta

        % Assigning control parameters from the ControllersSettings table
        [... % Using deal to assign multiple values at once
            MATDSS.Cont.CA(i).Area,... % Area of control area i
            MATDSS.Cont.CA(i).Type,... % Type of controller
            MATDSS.Cont.CA(i).alpha,... % Alpha parameter for control
            MATDSS.Cont.CA(i).rp,... % Control parameter rp
            MATDSS.Cont.CA(i).rbard,... % Control parameter rbar_d
            MATDSS.Cont.CA(i).E,... % Control parameter E
            MATDSS.Cont.CA(i).vul,... % Control parameter vul
            MATDSS.Cont.CA(i).vll,... % Control parameter vll
            MATDSS.Cont.CA(i).iul,... % Control parameter iul
            MATDSS.Cont.CA(i).arho,... % Dual variable arho
            MATDSS.Cont.CA(i).asigma,... % Dual variable asigma
            MATDSS.Cont.CA(i).alambda,... % Dual variable alambda
            MATDSS.Cont.CA(i).amu,... % Dual variable amu
            MATDSS.Cont.CA(i).aeta,... % Dual variable aeta
            MATDSS.Cont.CA(i).apsi,... % Dual variable apsi
            MATDSS.Cont.CA(i).agamma,... % Dual variable agamma
            MATDSS.Cont.CA(i).anu,... % Dual variable anu
            MATDSS.Cont.CA(i).azeta,... % Dual variable azeta
            MATDSS.Cont.CA(i).crho,... % Dual variable crho
            MATDSS.Cont.CA(i).csigma,... % Dual variable csigma
            MATDSS.Cont.CA(i).clambda,... % Dual variable clambda
            MATDSS.Cont.CA(i).cmu,... % Dual variable cmu
            MATDSS.Cont.CA(i).ceta,... % Dual variable ceta
            MATDSS.Cont.CA(i).cpsi,... % Dual variable cpsi
            MATDSS.Cont.CA(i).cgamma,... % Dual variable cgamma
            MATDSS.Cont.CA(i).cnu,... % Dual variable cnu
            MATDSS.Cont.CA(i).czeta] = deal(MATDSS.TableData.ControllersSettings{i,:}); % Assigning control parameters from table

        % Convert string parameters to numeric values where applicable
        ControllerFields = fieldnames(MATDSS.Cont.CA(i)); % Get field names of control area i
        for j = 1:27 % Loop through all fields
            if j ~= [2,3,9] % Skip certain indices
                MATDSS.Cont.CA(i).(ControllerFields{j}) = str2double(MATDSS.Cont.CA(i).(ControllerFields{j})); % Convert string to double
            end

            % Special case for alpha field
            if j == 3 && ~strcmp(MATDSS.Cont.CA(i).(ControllerFields{j}), 'auto') % Check if alpha is not set to 'auto'
                MATDSS.Cont.CA(i).(ControllerFields{j}) = str2double(MATDSS.Cont.CA(i).(ControllerFields{j})); % Convert string to double
            end
        end

        % Set alpha value based on conditions
        if strcmpi(MATDSS.Cont.CA(i).alpha, 'auto') % Check if alpha is set to 'auto'
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep; % Compute alpha automatically
        else
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.CA(i).alpha.*MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep; % Adjust alpha based on current value
        end

    end
end
end



%% Old function Code

%{

function [MATDSS, DER] = MATDSS_Controllers(MATDSS,DER)
% function MATDSS = MATDSS_Controllers(app,MATDSS,DER)
% This function will define the controllers needed following the provided
% information in the Configuration files. Use the configuration option in
% MATDSS Application to define the control areas and assign the buses to
% them. Also specify main controller parameters from within the application
% too.
%
% Make sure you specify the area to include all buses within it. Otherwise,
% the controller will be confused when defining the branches. This
% functionality will be autmated later to ensure the areas are being
% defined as they should be.
%
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ____                    _   _       _   _                  ____    _               ____            _            %
% |  _ \    ___    _ __   ( ) | |_    | | | |  ___    ___    / ___|  (_)  _ __ ___   |  _ \    __ _  | |_    __ _  %
% | | | |  / _ \  | '_ \  |/  | __|   | | | | / __|  / _ \   \___ \  | | | '_ ` _ \  | | | |  / _` | | __|  / _` | %
% | |_| | | (_) | | | | |     | |_    | |_| | \__ \ |  __/    ___) | | | | | | | | | | |_| | | (_| | | |_  | (_| | %
% |____/   \___/  |_| |_|      \__|    \___/  |___/  \___|   |____/  |_| |_| |_| |_| |____/   \__,_|  \__|  \__,_| %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~MATDSS.UseSimDataFlag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%     Parameters needed to define the conrollers    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    CA = [];
    AllBusNodesNames = {};
    AllNodesNames = MATDSS.Sim.Meas.AllNodesNames;
    AllBusNames = MATDSS.Sim.Meas.AllBusNames;
    for i = 1:length(AllNodesNames)
        iNodeName = AllNodesNames{i};
        iDot = strfind(iNodeName,'.'); % find the location of "dots" in the node name to find the busname
        iBusname = iNodeName(1:iDot(1)-1);
        AllBusNodesNames = [AllBusNodesNames;iBusname];
    end
    SourceBusName = AllBusNames{1};
    % Get all lines/branches names
    AllLinesNames = MATDSS.Sim.Meas.AllBranchesPhasesNames;
    % Get all lines/branches names
    MyLines = MATDSS.Sim.DSSCircuit.Lines;
    MyLinesAllNames = MyLines.AllNames;


    NumberOfCA = length(unique(str2double(MATDSS.TableData.CASettings(:,3))));


    MyPer = 0.55; MyPerEnd = 0.75;
    MyPerSec = (MyPerEnd - MyPer)/4;
    % MATDSSApp_Status(app,MyPer);

    % Extracting Control Areas settings from the tables.
    for i = 1:NumberOfCA
        Controller = struct;

        [...
            Controller.Area,...
            Controller.Type,...
            Controller.alpha,...
            Controller.rp,...
            Controller.rbard,...
            Controller.E,...
            Controller.vul,...
            Controller.vll,...
            Controller.iul,...
            Controller.arho,...
            Controller.asigma,...
            Controller.alambda,...
            Controller.amu,...
            Controller.aeta,...
            Controller.apsi,...
            Controller.agamma,...
            Controller.anu,...
            Controller.azeta,...
            Controller.crho,...
            Controller.csigma,...
            Controller.clambda,...
            Controller.cmu,...
            Controller.ceta,...
            Controller.cpsi,...
            Controller.cgamma,...
            Controller.cnu,...
            Controller.czeta] = deal(MATDSS.TableData.ControllersSettings{i,:});

        ControllerFields = fieldnames(Controller);
        % Changing strings to numbers in the corresponding locations.
        for j = 1:length(ControllerFields)
            if j ~= [2,3,9]
                Controller.(ControllerFields{j}) = str2double(Controller.(ControllerFields{j}));
            end

            if j == 3 && ~strcmp(Controller.(ControllerFields{j}),'auto')
                Controller.(ControllerFields{j}) = str2double(Controller.(ControllerFields{j}));
            end

        end

        Controller.Buses = cell(0,5); % Check which buses belongs to this controller
        % Adding list of buses that belong to the controller
        for j = 1:size(MATDSS.TableData.CASettings,1)
            if str2double(MATDSS.TableData.CASettings{j,3}) == i
                Controller.Buses = [Controller.Buses; MATDSS.TableData.CASettings(j,1),str2double(MATDSS.TableData.CASettings(j,2)),...
                    str2double(MATDSS.TableData.CASettings(j,3)),MATDSS.TableData.CASettings(j,4),...
                    MATDSS.TableData.CASettings(j,5)];
            end
        end


        % Check if this area is not L2C, then find the interface bus, and add
        % it to the list of buses.
        if MATDSS_StrComp(Controller.Buses(:,1),SourceBusName) <= 0
            CABus0 = Controller.Buses(find(strcmp('t',Controller.Buses(:,4))),5);
            CASettingsIndex = find(strcmp(MATDSS.TableData.CASettings(:,1),CABus0));
            Controller.Buses = [MATDSS.TableData.CASettings(CASettingsIndex,1),str2double(MATDSS.TableData.CASettings(CASettingsIndex,2)),...
                str2double(MATDSS.TableData.CASettings(CASettingsIndex,3)),MATDSS.TableData.CASettings(CASettingsIndex,4),...
                MATDSS.TableData.CASettings(CASettingsIndex,5);...
                Controller.Buses];
        end


        % Find all indecies of all phases in this area
        CABusesIndices = [];
        ControllerBusesIndices = [];
        for j = 2:size(Controller.Buses,1)
            CABusesIndices = [CABusesIndices; find(strcmp(AllBusNodesNames,Controller.Buses(j,1)))];
            if Controller.Buses{j,3} == i
                ControllerBusesIndices = [ControllerBusesIndices;find(strcmp(AllBusNodesNames,Controller.Buses(j,1)))];
            end
        end


        % Save these variables in the controller struct, then save it in CA
        % variable.
        Controller.NodesIndices = CABusesIndices;
        Controller.nBuses = size(Controller.Buses,1);
        Controller.ControllerBusesIndices = ControllerBusesIndices;


        CA = [CA; Controller];


        % MATDSSApp_Status(app, MyPer + MyPerSec * i / NumberOfCA);

        clear Controller

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Get All Interface Buses and Their Corresponding Information %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate informationm about virtual DERs and their location (to be used
    % mainly in ABMH function). We do that by looking at the branches, where
    % the buses belongs to different control areas.

    CAInterfaceBuses = {};
    % CA#, Bus 1, Bus1CA, Bus 2, Bus2CA, nPhases, Nodes
    % 1 CA#:is the Control area where the virtual DER is needed
    % 2 Bus1: is where the virtual DER is located.
    % 3 Bus1CA: is the control area where Bus1 is located (this is same as CA#)
    % 4 Bus2: is the bus connected to Bus1
    % 5 Bus2CA: is the Control area connected to the higher level control area. So
    % it is the one required to provide the information about P,Q and their
    % limits as a VDER.
    % 6 nPhases: number of phases at the interface bus. In some cases, this can
    % be less than 3 phases interface.
    % 7 Nodes: The nodes connected at Bus2 (to Bus1).
    %


    CABusesTIndex = find(strcmp(MATDSS.TableData.CASettings(:,4),'t'));
    CABus0 = {};

    % Loop over all interface buses/lines/connections
    for i = 1:length(CABusesTIndex)
        % Take bus 1 from teh table
        Bus1 = MATDSS.TableData.CASettings{CABusesTIndex(i),5};
        % Check to which area it belongs
        for j = 1:NumberOfCA
            if MATDSS_StrComp(CA(j).Buses(:,1),Bus1) > 0
                Bus1CA = j; % Get CA number of Bus1
                CANumber = j; % This is the area for which we are finding the connection line.
                break;
            end
        end

        % Get bus 2 from the table, also get its area too
        Bus2 = MATDSS.TableData.CASettings{CABusesTIndex(i),1};
        Bus2CA = MATDSS.TableData.CASettings{CABusesTIndex(i),3};
        if strcmpi(Bus2,SourceBusName) % If here we have it with the source bus, this means we are in L2C, simple!
            nPhases = 3; % Currently we are considering 3phase interface with transmission network. We can later add code to automate this with all possible cases.
            CABus0Phases = {};
            for j = 1:nPhases
                CABus0Phases = [CABus0Phases; [Bus2, '.', num2str(j)]];
            end
            CABus0 = [CABus0; {Bus2CA, MATDSS_StrComp(AllNodesNames,CABus0Phases)}]; % save it then in CABus0
        else % if this is not L2C (we are in lower levels)

            % In this case, we need to do more work. Let's find the line
            % connecting the bus0 and then get all phases information/order.
            for k = 1:length(MyLinesAllNames)
                % Loop over all lines in the circuit, to find the branch that
                % connects our interface bus of LLC to L2C or other LLC.
                MyLines.Name = MyLinesAllNames{k};
                LineBus1 = MyLines.Bus1; % Get bus 1 name
                iDotBus1 = strfind(LineBus1,'.'); % find the location of "dots" in the node name to find the busname
                if isempty(iDotBus1) % If the bus name we got was only the name without nodes, add the nodes manually and then remove them. This is important so our code can work with all cases.
                    LineBus1 = [LineBus1 '.1.2.3'];
                    iDotBus1 = strfind(LineBus1,'.'); % find the location of "dots" in the node name to find the busname
                end
                LineBus1FullName = LineBus1; % Keep a copy of the full name for later
                LineBus1 = LineBus1(1:iDotBus1(1)-1); % here is the short name without phases


                % Repeat what we did to Bus1 for Bus2.
                LineBus2 = MyLines.Bus2;
                iDotBus2 = strfind(LineBus2,'.'); % find the location of "dots" in the node name to find the busname
                if isempty(iDotBus2)
                    LineBus2 = [LineBus2 '.1.2.3'];
                    iDotBus2 = strfind(LineBus2,'.'); % find the location of "dots" in the node name to find the busname
                end
                LineBus2FullName = LineBus2;
                LineBus2 = LineBus2(1:iDotBus2(1)-1);


                % Now check, if one of those buses that we have on this branch
                % is the interface bus that we are interested in.
                if (strcmp(Bus1,LineBus1) && strcmp(Bus2,LineBus2)) || ((strcmp(Bus2,LineBus1) && strcmp(Bus1,LineBus2)))
                    nPhases = MyLines.Phases; % Get the number of phases of this CA

                    % Check if Bus1 is the interface bus we want, else flip
                    % them
                    if ~strcmp(Bus1,LineBus1)
                        LineBus1FullName = LineBus2FullName;
                        iDotBus1 = iDotBus2;
                    end

                    CABus0Phases = {};
                    % Get phases order
                    for j = 1:length(iDotBus1)
                        CABus0Phases = [CABus0Phases; [LineBus1FullName(1:iDotBus1(1)) LineBus1FullName(iDotBus1(j)+1)]];
                    end

                    % Save Interface bus information in CABus0
                    CABus0 = [CABus0; {Bus2CA, MATDSS_StrComp(AllNodesNames,CABus0Phases)}];
                    break;
                end
            end % Done with detecting buses and lines that connects interface bus of CA we are considering here

            % generating list of nodes to build our CAInterfaceBuses variable.
            bus1nodes = '[';
            for j = 1:length(iDotBus1)
                bus1nodes = [bus1nodes, LineBus1FullName(iDotBus1(j)+1),','];
            end
            bus1nodes = [bus1nodes(1:end-1),']'];

            CAInterfaceBuses = [CAInterfaceBuses; {CANumber, Bus1, Bus1CA, Bus2, str2num(Bus2CA), nPhases, bus1nodes}];
        end
        % MATDSSApp_Status(app, MyPer + MyPerSec*(1 + i / length(CABusesTIndex)));

    end % End of our for loop over all interface buses


    % Save CA and CAInterfaceBuses variable in MATDSS.Cont
    MATDSS.Cont.CA = CA;
    MATDSS.Cont.CAInterfaceBuses = CAInterfaceBuses;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                Defing Virtual DERs                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % So now, we are ready to define our VDERs!
    for i = 1:size(CAInterfaceBuses,1)

        VDERsTable = MATDSS.TableData.VDERSettings;
        VDERIndexTable = find(CAInterfaceBuses{i,5}==str2double(VDERsTable(:,1)));
        VDERInfo = VDERsTable(VDERIndexTable,[3,4,5,6,7,8,9]);


        [MATDSS, DER] = MATDSS_DERNew(MATDSS,DER,...
            num2str(size(DER,2)+1),...                                                              DER Index
            ['VDER_CA' num2str(CAInterfaceBuses{i,1}) '_For_CA' num2str(CAInterfaceBuses{i,5})],... DER Name
            strcat("'", CAInterfaceBuses{i,2}, "'"),...                                             DER Bus Location
            VDERInfo(1),...                                                                              Tau
            "VDER",...                                                                              DER Type --> VDER
            CAInterfaceBuses{i,7},...                                                               Nodes for the VDER
            VDERInfo(2),...                                                                               Type of connection (for now will go with Wye, but we might need to have a smarter way of selecting why or delta)
            num2str(CAInterfaceBuses{i,6}),...                                                      Number of phases
            VDERInfo(3),...                                                                          DER Mode
            VDERInfo(4),...                                                                          P(x)
            VDERInfo(5),...                                                                          Q(x)
            '-1e99','1e99','-1e99','1e99',...                                                       Pmin,Pmax,Qmin,Qmax
            VDERInfo(6),VDERInfo(7));                                                                              %ax and cx

        % MATDSSApp_Status(app, MyPer + MyPerSec*(2 + i / size(CAInterfaceBuses,1)));
    end

    % Now, we are ready to re-define our control areas (as isolated networks).
    % We define our CA-DER list for each CA. We map the location of the DERs
    % (DERs and VDERs) to the bus numbers inside the redefined control areas.
    %
    % Then we will be able to call our ABMH matrix. We define a new measurement
    % function that would map the voltages, currents and powers to the new
    % isolated networks.


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%         GET CA Properties and Information         %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:NumberOfCA


        % Initialize some properties
        CAnDER = 0; %nDER of CA
        CADERIndex = []; % Corresponding indecies
        DERBuses = [DER(:).BusNum]'; % DER Buses
        %Pmax, Pmin, Qmax, Qmin of CA VDER. Currenlty we sum up all available
        %resources within CA to get an estimate of available energy.
        %
        % This section needs to be rewritten to start from lower areas up to
        % L2C
        %
        % If DER is in the area, add up its limits to our VDER
        Pmax = 0;
        Pmin = 0;
        Qmax = 0;
        Qmin = 0;
        for j = 1:length(DERBuses)
            IsDERInCA = find(DERBuses(j) == [MATDSS.Cont.CA(i).Buses{:,2}]');
            if ~isempty(IsDERInCA)
                if MATDSS.Cont.CA(i).Buses{IsDERInCA,3} == i
                    CAnDER = CAnDER + 1;
                    CADERIndex = [CADERIndex; j];
                    DER(j).CAIndex = i;
                    Pmin = Pmin - DER(j).Pmax;
                    Pmax = Pmax - DER(j).Pmin;
                    Qmin = Qmin - DER(j).Qmax;
                    Qmax = Qmax - DER(j).Qmin;
                end
            end
        end

        MATDSS.Cont.CA(i).nDER = CAnDER; %Save the number of DERs in CA
        MATDSS.Cont.CA(i).DERIndex = CADERIndex; % Save DERs indecies in CA
        area_index = find (i==str2double([CABus0(:,1)])); % CA index in the list of CABus0
        MATDSS.Cont.CA(i).CABus0 = CABus0{area_index,2}; % Get corresponding Bus0 information
        MATDSS.Cont.CA(i).NodesIndices = [MATDSS.Cont.CA(i).CABus0;MATDSS.Cont.CA(i).NodesIndices]; % All nodes indicies including Bus0
        MATDSS.Cont.CA(i).nNodes = length(MATDSS.Cont.CA(i).NodesIndices); %nNodes = length of all nodes in CA
        MATDSS.Cont.CA(i).nBuses = size(MATDSS.Cont.CA(i).Buses,1); % nBuses = length of Buses in CA



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%                    YBus of CA                     %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % We use Yreduced method to shrink Ybus of the whole network to
        % Yreduced of CA
        FullYbus = MATDSS.Sim.Ybus;
        ControllerYbus = [];
        ControllerYbus = [FullYbus(:,MATDSS.Cont.CA(i).NodesIndices),FullYbus(:,setdiff(1:length(AllBusNodesNames),MATDSS.Cont.CA(i).NodesIndices))];
        ControllerYbus = [ControllerYbus(MATDSS.Cont.CA(i).NodesIndices,:);ControllerYbus(setdiff(1:length(AllBusNodesNames),MATDSS.Cont.CA(i).NodesIndices),:)];

        nphase = length(AllBusNodesNames);
        nCAphase = length(MATDSS.Cont.CA(i).NodesIndices);
        Ynn = ControllerYbus(1:nCAphase,1:nCAphase);
        Yrr = ControllerYbus(nCAphase+1:nphase,nCAphase+1:nphase);
        Ynr = ControllerYbus(1:nCAphase,nCAphase+1:nphase);
        Yrn = ControllerYbus(nCAphase+1:nphase,1:nCAphase);
        ControllerYred = Ynr/(Yrr);
        ControllerYred = Ynn-ControllerYred*Yrn;
        % ControllerYred = Ynn-Ynr*(Yrr^-1)*Yrn;
        MATDSS.Cont.CA(i).Ybus = ControllerYred;




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%      Find P0, Q0, P0Set, Q0Set of CA and VDER     %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % If this is L2C (highest control area), then no VDER is needed.
        if strcmp(MATDSS.Cont.CA(i).Buses(1,1),SourceBusName)
            MATDSS.Cont.CA(i).P0Set = zeros(1,length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).Q0Set = zeros(1,length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).TimeIndex = 1;
            temp = MATDSS.Cont.CA(i).CABus0;
            while ~isempty(temp)
                MATDSS.Cont.CA(i).ControllerBusesIndices(MATDSS.Cont.CA(i).ControllerBusesIndices == temp(1)) = [];
                temp(1) = [];
            end
            MATDSS.Cont.CA(i).VDERIndex = -1; % Indicating that we should use The P0 and Q0 from MATDSS.ControlSignals
        else % Get the VDER index for the CA
            for j = 1:size(CAInterfaceBuses,1) % Record the corresponding VDER index for the CA to retrieve the P and Q Setpoints
                if CAInterfaceBuses{j,5} == MATDSS.Cont.CA(i).Area
                    nNotVDER = MATDSS.Sim.nDER - size(CAInterfaceBuses,1);
                    MATDSS.Cont.CA(i).VDERIndex = nNotVDER + j;

                    if MATDSS.Cont.CA(i).VDERIndex > 0 % This if statement is not needed
                        % Foimd
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
                                %                             if strcmpi(MyLineBuses{1},CAInterfaceBuses{j,2})
                                %                                 MyLineI = MyLineI(1:length(MyLineI)/2).';
                                %                                 break;
                                %                             else
                                %                                 MyLineI = MyLineI(length(MyLineI)/2 + 1:end).';
                                %                                 break;
                                %                             end
                            end
                        end

                        MyCAInterfaceS0 = MyBusV.*conj(MyLineI);
                        MyCAInterfaceP0 = real(MyCAInterfaceS0);
                        MyCAInterfaceQ0 = imag(MyCAInterfaceS0);
                        MATDSS.Cont.CA(i).MyCAInterfaceP0 = MyCAInterfaceP0;
                        MATDSS.Cont.CA(i).MyCAInterfaceQ0 = MyCAInterfaceQ0;
                        i_WVar = DER(MATDSS.Cont.CA(i).VDERIndex).i_WVar;
                        DER(MATDSS.Cont.CA(i).VDERIndex).W(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).W,1),1).*(-sum(MyCAInterfaceP0));
                        DER(MATDSS.Cont.CA(i).VDERIndex).Var(:,2)  = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).Var,1),1).*(-sum(MyCAInterfaceQ0));
                        DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint,1),1).*(-sum(MyCAInterfaceP0));
                        DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint,1),1).*(-sum(MyCAInterfaceQ0));
                        MATDSS.Cont.CA(i).P0Set = ones(1,length(MATDSS.Cont.at)).*(sum(MyCAInterfaceP0));
                        MATDSS.Cont.CA(i).Q0Set = ones(1,length(MATDSS.Cont.at)).*(sum(MyCAInterfaceQ0));


                        DER(MATDSS.Cont.CA(i).VDERIndex).Pmax = -sum(MyCAInterfaceP0) + Pmax;
                        DER(MATDSS.Cont.CA(i).VDERIndex).Pmin = -sum(MyCAInterfaceP0) + Pmin;
                        DER(MATDSS.Cont.CA(i).VDERIndex).Qmax = -sum(MyCAInterfaceQ0) + Qmax;
                        DER(MATDSS.Cont.CA(i).VDERIndex).Qmin = -sum(MyCAInterfaceQ0) + Qmin;
                        MATDSS.Cont.CA(i).TimeIndex = 1;

                        %                     DER(MATDSS.Cont.CA(i).VDERIndex).Pmax = DER(MATDSS.Cont.CA(i).VDERIndex).W + 10e3;
                        %                     DER(MATDSS.Cont.CA(i).VDERIndex).Pmin = DER(MATDSS.Cont.CA(i).VDERIndex).W + -10e3;
                        %                     DER(MATDSS.Cont.CA(i).VDERIndex).Qmax = DER(MATDSS.Cont.CA(i).VDERIndex).Var + 10e3;
                        %                     DER(MATDSS.Cont.CA(i).VDERIndex).Qmin = DER(MATDSS.Cont.CA(i).VDERIndex).Var + -10e3;
                    end

                end
            end

            %         MATDSS.Cont.CA(i).CABus0 = setdiff(MATDSS.Cont.CA(i).NodesIndices,MATDSS.Cont.CA(i).ControllerBusesIndices);
        end


        %     MATDSS.Cont.CA(i).CABus0 = CABus0;
        [MATDSS.Cont.CA(i).k_v, MATDSS.Cont.CA(i).k_vIndex,~] = intersect(MATDSS.Meas.k_v,MATDSS.Cont.CA(i).ControllerBusesIndices);


        % Find k_i
        k_i = [];
        if MATDSS.Meas.k_i > 0
            k_iLines = AllLinesNames(MATDSS.Meas.k_i);
            for j = 1:length(k_iLines)
                iDotline = strfind(k_iLines{j},'.');
                LineName = k_iLines{j};
                MyLines.Name = LineName(1:iDotline(1)-1); % Get the line name first and handle that line in mylines variable
                LineBus1 = MyLines.Bus1; % Get bus 1 name with all phases
                LineBus2 = MyLines.Bus2; % Get bus 2 name with all phases

                iDotBus1 = strfind(LineBus1,'.');
                iDotBus2 = strfind(LineBus2,'.');
                if ~isempty(iDotBus1)
                    LineBus1 = LineBus1(1:iDotBus1(1)-1);
                end
                if ~isempty(iDotBus2)
                    LineBus2 = LineBus2(1:iDotBus2(1)-1);
                end
                iBus1 = find(strcmp(MATDSS.Cont.CA(i).Buses(:,1),LineBus1));
                iBus2 = find(strcmp(MATDSS.Cont.CA(i).Buses(:,1),LineBus2));

                if ~isempty(iBus1) && ~isempty(iBus2)
                    k_i = [k_i; MATDSS.Meas.k_i(j)];
                end

            end
        end
        MATDSS.Cont.CA(i).k_i = k_i;
        [~,MATDSS.Cont.CA(i).k_iIndex,~] = intersect(MATDSS.Meas.k_i,MATDSS.Cont.CA(i).k_i);


        [A, B, M, H, Mv, MvIndex, Mi,MiBuses] = MATDSS_ABMH(MATDSS,DER,MATDSS.Cont.CA(i));
        MATDSS.Cont.CA(i).ABMH.A = A;
        MATDSS.Cont.CA(i).ABMH.B = B;
        MATDSS.Cont.CA(i).ABMH.M = M;
        MATDSS.Cont.CA(i).ABMH.H = H;
        % Initialize the controller dual variables
        MATDSS.Cont.CA(i).Duals.rho = zeros(size(M,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.sigma = zeros(size(H,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.lambda = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.mu = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.eta = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.psi = zeros(1,length(MATDSS.Cont.at));
        %     if MvIndex == 0
        %         MATDSS.Cont.CA(i).Duals.gamma(:,1) = 0;
        %         MATDSS.Cont.CA(i).Duals.nu(:,1) = 0;
        %     else
        MATDSS.Cont.CA(i).Duals.gamma = zeros(max(size(Mv)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.nu = zeros(max(size(Mv)),length(MATDSS.Cont.at));
        %     end
        % if isempty(k_i)
        % MATDSS.Cont.CA(i).Duals.zeta(:,1) = 0;
        % else
        MATDSS.Cont.CA(i).Duals.zeta = zeros(max(size(Mi)),length(MATDSS.Cont.at));
        % end



        % Defining error handles
        MATDSS.Cont.CA(i).Duals.Errors.rho = zeros(size(M,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.sigma = zeros(size(H,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.lambda = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.mu = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.eta = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.psi = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.gamma = zeros(max(size(Mv)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.nu = zeros(max(size(Mv)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.zeta = zeros(max(size(Mi)),length(MATDSS.Cont.at));


        % Defining error handles for measurements in D controller
        MATDSS.Cont.CA(i).Duals.ErrorsY.rho = zeros(size(M,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.sigma = zeros(size(H,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.lambda = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.mu = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.eta = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.psi = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.gamma = zeros(max(size(Mv)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.nu = zeros(max(size(Mv)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.zeta = zeros(max(size(Mi)),length(MATDSS.Cont.at));

        if strcmpi(MATDSS.Cont.CA(i).alpha,'auto')
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
        else
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.CA(i).alpha.*MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
        end
        % MATDSS.Cont.CA(i).ABMH = []; % Check why!
        %     MATDSS.Cont.CA(i).P0Set = [];
        %     MATDSS.Cont.CA(i).Q0Set = [];

        MATDSS.Cont.CA(i).P0 = zeros(length(MATDSS.Cont.CA(i).CABus0),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Q0 = zeros(length(MATDSS.Cont.CA(i).CABus0),length(MATDSS.Cont.at));
        % MATDSSApp_Status(app, MyPer + MyPerSec*(3 + i / NumberOfCA));
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  _   _                  ____    _               ____            _            %
    % | | | |  ___    ___    / ___|  (_)  _ __ ___   |  _ \    __ _  | |_    __ _  %
    % | | | | / __|  / _ \   \___ \  | | | '_ ` _ \  | | | |  / _` | | __|  / _` | %
    % | |_| | \__ \ |  __/    ___) | | | | | | | | | | |_| | | (_| | | |_  | (_| | %
    %  \___/  |___/  \___|   |____/  |_| |_| |_| |_| |____/   \__,_|  \__|  \__,_| %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    % Copying Controllers from SimData
    MATDSS.Cont.CA = MATDSS.LoadedSimData.MATDSS.Cont.CA;
    MATDSS.Cont.CAInterfaceBuses = MATDSS.LoadedSimData.MATDSS.Cont.CAInterfaceBuses;


    AllBusNames = MATDSS.Sim.Meas.AllBusNames;
    SourceBusName = AllBusNames{1};
    CAInterfaceBuses = MATDSS.Cont.CAInterfaceBuses;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%                Defing Virtual DERs                %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % So now, we are ready to define our VDERs!
    for i = 1:size(CAInterfaceBuses,1)

        VDERsTable = MATDSS.TableData.VDERSettings;
        VDERIndexTable = find(CAInterfaceBuses{i,5}==str2double(VDERsTable(:,1)));
        VDERInfo = VDERsTable(VDERIndexTable,[3,4,5,6,7,8,9]);


        [MATDSS, DER] = MATDSS_DERNew(MATDSS,DER,...
            num2str(size(DER,2)+1),...                                                                       DER Index
            ['VDER_CA' num2str(CAInterfaceBuses{i,1}) '_For_CA' num2str(CAInterfaceBuses{i,5})],... DER Name
            strcat("'", CAInterfaceBuses{i,2}, "'"),...                                             DER Bus Location
            VDERInfo(1),...                                                                              Tau
            "VDER",...                                                                              DER Type --> VDER
            CAInterfaceBuses{i,7},...                                                               Nodes for the VDER
            VDERInfo(2),...                                                                               Type of connection (for now will go with Wye, but we might need to have a smarter way of selecting why or delta)
            num2str(CAInterfaceBuses{i,6}),...                                                      Number of phases
            VDERInfo(3),...                                                                          DER Mode
            VDERInfo(4),...                                                                          P(x)
            VDERInfo(5),...                                                                          Q(x)
            '-1e99','1e99','-1e99','1e99',...                                                       Pmin,Pmax,Qmin,Qmax
            VDERInfo(6),VDERInfo(7));                                                                              %ax and cx
    end


    NumberOfCA = length(unique(str2double(MATDSS.TableData.CASettings(:,3))));
    DERBuses = [DER(:).BusNum]';
    for i = 1:NumberOfCA
        CAnDER = 0; %nDER of CA
        CADERIndex = []; % Corresponding indecies
        % If DER is in the area, add up its limits to our VDER
        Pmax = 0;
        Pmin = 0;
        Qmax = 0;
        Qmin = 0;
        for j = 1:length(DERBuses)
            IsDERInCA = find(DERBuses(j) == [MATDSS.Cont.CA(i).Buses{:,2}]');
            if ~isempty(IsDERInCA)
                if MATDSS.Cont.CA(i).Buses{IsDERInCA,3} == i
                    CAnDER = CAnDER + 1;
                    CADERIndex = [CADERIndex; j];
                    DER(j).CAIndex = i;
                    Pmin = Pmin - DER(j).Pmax;
                    Pmax = Pmax - DER(j).Pmin;
                    Qmin = Qmin - DER(j).Qmax;
                    Qmax = Qmax - DER(j).Qmin;
                end
            end
        end
        if strcmp(MATDSS.Cont.CA(i).Buses(1,1),SourceBusName)
            MATDSS.Cont.CA(i).P0Set = zeros(1,length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).Q0Set = zeros(1,length(MATDSS.Cont.at));
            MATDSS.Cont.CA(i).TimeIndex = 1;
        else
            MyCAInterfaceP0 = MATDSS.Cont.CA(i).MyCAInterfaceP0;
            MyCAInterfaceQ0 = MATDSS.Cont.CA(i).MyCAInterfaceQ0;

            MATDSS.Cont.CA(i).P0Set = ones(1,length(MATDSS.Cont.at)).*(sum(MyCAInterfaceP0));
            MATDSS.Cont.CA(i).Q0Set = ones(1,length(MATDSS.Cont.at)).*(sum(MyCAInterfaceQ0));

            DER(MATDSS.Cont.CA(i).VDERIndex).W(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).W,1),1).*(-sum(MyCAInterfaceP0));
            DER(MATDSS.Cont.CA(i).VDERIndex).Var(:,2)  = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).Var,1),1).*(-sum(MyCAInterfaceQ0));
            DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint,1),1).*(-sum(MyCAInterfaceP0));
            DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint(:,2) = ones(size(DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint,1),1).*(-sum(MyCAInterfaceQ0));
            DER(MATDSS.Cont.CA(i).VDERIndex).Pmax = -sum(MyCAInterfaceP0) + Pmax;
            DER(MATDSS.Cont.CA(i).VDERIndex).Pmin = -sum(MyCAInterfaceP0) + Pmin;
            DER(MATDSS.Cont.CA(i).VDERIndex).Qmax = -sum(MyCAInterfaceQ0) + Qmax;
            DER(MATDSS.Cont.CA(i).VDERIndex).Qmin = -sum(MyCAInterfaceQ0) + Qmin;
        end

        MATDSS.Cont.CA(i).Duals.rho = zeros(size(MATDSS.Cont.CA(i).Duals.rho,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.sigma = zeros(size(MATDSS.Cont.CA(i).Duals.sigma,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.lambda = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.mu = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.eta = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.psi = zeros(1,length(MATDSS.Cont.at));
        %     if MvIndex == 0
        %         MATDSS.Cont.CA(i).Duals.gamma(:,1) = 0;
        %         MATDSS.Cont.CA(i).Duals.nu(:,1) = 0;
        %     else
        MATDSS.Cont.CA(i).Duals.gamma = zeros(max(size(MATDSS.Cont.CA(i).Duals.gamma,1)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.nu = zeros(max(size(MATDSS.Cont.CA(i).Duals.nu,1)),length(MATDSS.Cont.at));
        %     end
        % if isempty(k_i)
        % MATDSS.Cont.CA(i).Duals.zeta(:,1) = 0;
        % else
        MATDSS.Cont.CA(i).Duals.zeta = zeros(max(size(MATDSS.Cont.CA(i).Duals.zeta,1)),length(MATDSS.Cont.at));
        % end


        % if strcmpi(MATDSS.Cont.CA(i).alpha,'auto')
        %     MATDSS.Cont.CA(i).alpha = MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
        % else
        %     MATDSS.Cont.CA(i).alpha = MATDSS.Cont.CA(i).alpha.*MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
        % end
        % MATDSS.Cont.CA(i).ABMH = []; % Check why!
        %     MATDSS.Cont.CA(i).P0Set = [];
        %     MATDSS.Cont.CA(i).Q0Set = [];

        MATDSS.Cont.CA(i).P0 = zeros(length(MATDSS.Cont.CA(i).CABus0),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Q0 = zeros(length(MATDSS.Cont.CA(i).CABus0),length(MATDSS.Cont.at));

        % Defining error handles
        MATDSS.Cont.CA(i).Duals.Errors.rho = zeros(size(MATDSS.Cont.CA(i).Duals.rho,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.sigma = zeros(size(MATDSS.Cont.CA(i).Duals.sigma,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.lambda = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.mu = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.eta = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.psi = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.gamma = zeros(max(size(MATDSS.Cont.CA(i).Duals.gamma,1)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.nu = zeros(max(size(MATDSS.Cont.CA(i).Duals.nu,1)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.Errors.zeta = zeros(max(size(MATDSS.Cont.CA(i).Duals.zeta,1)),length(MATDSS.Cont.at));

        % Defining error handles for Y measurements in D controller
        MATDSS.Cont.CA(i).Duals.ErrorsY.rho = zeros(size(MATDSS.Cont.CA(i).Duals.rho,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.sigma = zeros(size(MATDSS.Cont.CA(i).Duals.sigma,1),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.lambda = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.mu = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.eta = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.psi = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.gamma = zeros(max(size(MATDSS.Cont.CA(i).Duals.gamma,1)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.nu = zeros(max(size(MATDSS.Cont.CA(i).Duals.nu,1)),length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Duals.ErrorsY.zeta = zeros(max(size(MATDSS.Cont.CA(i).Duals.zeta,1)),length(MATDSS.Cont.at));

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
            MATDSS.Cont.CA(i).czeta] = deal(MATDSS.TableData.ControllersSettings{i,:});

        ControllerFields = fieldnames(MATDSS.Cont.CA(i));
        % Changing strings to numbers in the corresponding locations.
        for j = 1:27
            if j ~= [2,3,9]
                MATDSS.Cont.CA(i).(ControllerFields{j}) = str2double(MATDSS.Cont.CA(i).(ControllerFields{j}));
            end

            if j == 3 && ~strcmp(MATDSS.Cont.CA(i).(ControllerFields{j}),'auto')
                MATDSS.Cont.CA(i).(ControllerFields{j}) = str2double(MATDSS.Cont.CA(i).(ControllerFields{j}));
            end

        end


        if strcmpi(MATDSS.Cont.CA(i).alpha,'auto')
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
        else
            MATDSS.Cont.CA(i).alpha = MATDSS.Cont.CA(i).alpha.*MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
        end
    end

end

end




%}