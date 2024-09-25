function MATDSS = MATDSS_Disturbance(MATDSS)
% MATDSS_Disturbance(MATDSS)
% This function generates disturbances in OpenDSS based on the
% configurations defined in the MATDSS application. It assigns
% various parameters such as load type (wye or delta) and
% initializes the disturbances for the simulation.
%
% Parameters:
%   - MATDSS: The MATDSS application structure containing properties
%             and methods related to the simulation and disturbance
%             configurations.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com

DisturbanceTableData = MATDSS.TableData.Disturbance; % Retrieve disturbance data

nDis = size(DisturbanceTableData, 1); % Number of Disturbances

MATDSS.Sim.nDis = nDis; % Update the number of disturbances in simulation
for i = 1:nDis
    ID = str2double(DisturbanceTableData(i, 1)); % Extract ID for disturbance
    Bus = DisturbanceTableData(i, 2); % Extract bus name/number
    Phases = str2num(DisturbanceTableData(i, 3)); % Extract phases as numeric array
    NodesText = []; % Initialize nodes text for OpenDSS command
    myNodes = Phases; % Temporary variable to manipulate phases

    % Constructing the NodesText for OpenDSS syntax
    while ~isempty(myNodes)
        NodesText = [NodesText, '.', num2str(myNodes(1))]; % Append node number
        myNodes = myNodes(2:end); % Remove the first element
    end

    Nphase = str2num(DisturbanceTableData(i, 4)); % Extract number of phases

    % Determine connection type (Wye or Delta)
    if strcmpi(DisturbanceTableData(i, 5), 'y')
        ConnType = 'Wye'; % Set connection type to Wye
    else
        ConnType = 'Delta'; % Set connection type to Delta
    end

    P = str2num(DisturbanceTableData(i, 6)); % Extract active power
    pf = str2num(DisturbanceTableData(i, 7)); % Extract power factor

    % Saving time information for disturbance activation
    t = str2num(DisturbanceTableData(i, 8));

    % Determine if the disturbance is enabled
    if strcmpi(DisturbanceTableData(i, 9), 't')
        Enabled = true; % Set enabled flag to true
    else
        Enabled = false; % Set enabled flag to false
    end

    % Construct DSS name for the disturbance
    BusChar = convertStringsToChars(Bus); % Convert bus name to character array
    DSSName = strcat('Disturb_Load_', num2str(ID), '_bus_', BusChar(2:end-1)); % Create unique name
    % Store disturbance parameters in MATDSS structure
    MATDSS.Disturbance(i).DSSName = DSSName;
    MATDSS.Disturbance(i).ID = ID;
    MATDSS.Disturbance(i).Bus = Bus;
    MATDSS.Disturbance(i).Phases = Phases;
    MATDSS.Disturbance(i).Nphase = Nphase;
    MATDSS.Disturbance(i).ConnType = ConnType;
    MATDSS.Disturbance(i).P = P;
    MATDSS.Disturbance(i).pf = pf;
    MATDSS.Disturbance(i).t = t;
    MATDSS.Disturbance(i).Enabled = Enabled;

    % Determine the bus number for OpenDSS
    charBus = convertStringsToChars(Bus); % Convert bus to character array
    if charBus(1) == char(39) % Check if bus name is given (enclosed in single quotes)
        BusNum = MATDSS_StrComp(MATDSS.Sim.DSSCircuit.AllBusNames, charBus(2:end-1)); % Get bus number by name
    else % If a bus number is provided
        BusNum = str2double(Bus); % Convert bus number to double
    end

    % Select the corresponding bus in OpenDSS
    myBus = MATDSS.Sim.DSSCircuit.Buses(int16(BusNum - 1)); % Access bus object by index

    % Create the OpenDSS command to define the disturbance
    DSSNewDisturbanceCommand = ['New Load.' DSSName '  Bus1=' BusChar(2:end-1) NodesText '  Phases=' num2str(Nphase) ...
        '  Conn=' ConnType '  Model=1 kV=' num2str(myBus.kVBase * sqrt(3)) '  kW=0  pf=' num2str(pf) ' Vminpu=0.1']; % Command string

    MATDSS.Sim.DSSText.Command = DSSNewDisturbanceCommand; % Assign the command for execution
end

end

%% Old function code

%{

function MATDSS = MATDSS_Disturbance(MATDSS)
% MATDSS = MATDSS_Disturbance(app, MATDSS)
%
% This function will generate loads in OpenDSS following the
% details/configurations in MATDSS Application. The function will assign
% load type (wye-delta)
%




DisturbanceTableData = MATDSS.TableData.Disturbance;

nDis = size(DisturbanceTableData,1); % Number of Disturbances

MATDSS.Sim.nDis = nDis;
for i = 1:nDis
    ID = str2double(DisturbanceTableData(i,1));
    Bus = DisturbanceTableData(i,2);
    Phases = str2num(DisturbanceTableData(i,3));
    NodesText = [];
    myNodes = Phases;
    while ~isempty(myNodes)
        NodesText = [NodesText,'.',num2str(myNodes(1))];
        myNodes = myNodes(2:end);
    end

    Nphase = str2num(DisturbanceTableData(i,4));
    

    if strcmpi(DisturbanceTableData(i,5),'y')
        ConnType = 'Wye';
    else
        ConnType = 'Delta';
    end

    P = str2num(DisturbanceTableData(i,6));
    pf = str2num(DisturbanceTableData(i,7));

    % Saving time information

    t = str2num(DisturbanceTableData(i,8));

    if strcmpi(DisturbanceTableData(i,9),'t')
        Enabled = true;
    else
        Enabled = false;
    end


    
    BusChar = convertStringsToChars(Bus);
    DSSName = strcat('Disturb_Load_', num2str(ID), '_bus_', BusChar(2:end-1));
    MATDSS.Disturbance(i).DSSName = DSSName;
    MATDSS.Disturbance(i).ID = ID;
    MATDSS.Disturbance(i).Bus = Bus;
    MATDSS.Disturbance(i).Phases = Phases;
    MATDSS.Disturbance(i).Nphase = Nphase;
    MATDSS.Disturbance(i).ConnType = ConnType;
    MATDSS.Disturbance(i).P = P;
    MATDSS.Disturbance(i).pf = pf;
    MATDSS.Disturbance(i).t = t;
    MATDSS.Disturbance(i).Enabled = Enabled;

    charBus = convertStringsToChars(Bus);
    % Getting bus number (check if bus name is given)
    if charBus(1) == char(39) %string of bus name
        BusNum = MATDSS_StrComp(MATDSS.Sim.DSSCircuit.AllBusNames,charBus(2:end-1));
    else %string of bus number
        BusNum = str2double(Bus);
    end
    %Select the bus in OpenDSS to define the DER OpenDSS load
    myBus = MATDSS.Sim.DSSCircuit.Buses(int16(BusNum-1));

    DSSNewDisturbanceCommand = ['New Load.' DSSName '  Bus1=' BusChar(2:end-1) NodesText '  Phases=' num2str(Nphase) ...
    '  Conn=' ConnType  '  Model=1 kV=' num2str(myBus.kVBase*sqrt(3)) '  kW=0  pf=' num2str(pf) ' Vminpu=0.1'];% pf=1'];

    MATDSS.Sim.DSSText.Command = DSSNewDisturbanceCommand;


    
end


end
%}