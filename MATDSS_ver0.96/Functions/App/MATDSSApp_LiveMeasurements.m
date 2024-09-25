function MATDSSApp_LiveMeasurements(app, MATDSS, DER)
% MATDSSApp_LiveMeasurements(app, MATDSS, DER)
% This function is responsible for displaying live measurements of various
% parameters such as voltage (V), current (I), active power (P), and reactive
% power (Q) for all buses and DERs (Distributed Energy Resources) in the
% MATDSS system. It focuses on the DERs currently under control.
%
% Measurements include variables such as voltages, currents, power flow,
% and dual variables (used for optimization/control purposes). The function
% updates the respective tables within the app's GUI to reflect real-time
% data during the simulation process.
%
% Parameters:
%   - app: The MATDSS application instance containing GUI properties and methods.
%   - MATDSS: The MATDSS system object that contains simulation variables and
%             measurement data.
%   - DER: The specific Distributed Energy Resource being measured.
%
% This function is executed when the "Live Measurements" button is activated
% in the app interface.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Get the current simulation time
t = MATDSS.Sim.t;   % Current simulation time step
at = MATDSS.Sim.at; % Absolute time in the simulation

% Check if live measurements are enabled by the user
if app.LiveMeasButton.Value
    % Voltage Measurements Update
    % Initialize voltage table headers and variable types
    VoltageTableHeaders = {'Vname', 'CA', 'Vll (kV)', 'V (kV)', 'Vul (kV)', 'gamma (Vul)', 'nu (Vll)'};
    VoltageTableVarTypes = repmat({'string'}, 1, length(VoltageTableHeaders)); % Define the data types for table columns
    VoltageTableData = {}; % Empty table for voltage data to be filled later

    % Determine the delay index based on measurement delay settings
    DelayIndex = ceil((MATDSS.Time.Meas.Delay / MATDSS.Time.Meas.TimeStep));

    % Calculate the number of valid measurements based on simulation progress
    pos = MATDSS.Meas.at - at;
    pos(pos > 0) = []; % Remove future times
    Meast = length(pos) - DelayIndex; % Number of valid measurement steps

    % Loop through each Control Area (CA) to update voltage measurements
    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i); % Get the control area

        % Loop through all voltage nodes within the current control area
        for j = 1:length(CA.k_v)
            Vname = MATDSS.Sim.Meas.AllNodesNames{CA.k_v(j)}; % Node name
            CAArea = CA.Area; % Control area identifier
            VBase = MATDSS.Meas.VBases(CA.k_v(j)); % Base voltage at this node
            Vll = CA.vll * VBase; % Line-to-line voltage
            V = MATDSS.Meas.VMagProfile(CA.k_vIndex(j), Meast) / 1e3; % Magnitude of voltage in kV
            Vul = CA.vul * VBase; % Upper limit voltage
            gamma = CA.Duals.gamma(j, CA.TimeIndex); % Dual variable gamma
            nu = CA.Duals.nu(j, CA.TimeIndex); % Dual variable nu

            % Store the voltage data in the table
            VoltageTableData = [VoltageTableData; {Vname, CAArea, Vll, V, Vul, gamma, nu}];
        end
    end

    % Loop through all nodes to update general voltage measurements
    for i = 1:length(MATDSS.Sim.Meas.AllNodesNames)
        Vname = MATDSS.Sim.Meas.AllNodesNames{i}; % Node name
        CAArea = -1; % General nodes not associated with a control area
        VBase = MATDSS.Meas.VBases(i); % Base voltage
        Vll = CA.vll * VBase; % Line-to-line voltage
        V = MATDSS.Sim.Meas.VMagProfile(i, t) / 1e3; % Voltage magnitude in kV
        Vul = CA.vul * VBase; % Upper limit voltage
        gamma = -1; % Default value for general nodes
        nu = -1; % Default value for general nodes

        % Store the voltage data in the table
        VoltageTableData = [VoltageTableData; {Vname, CAArea, Vll, V, Vul, gamma, nu}];
    end

    % Create the voltage table and populate the GUI table component
    VoltageTable = table('Size', size(VoltageTableData), 'VariableTypes', VoltageTableVarTypes, 'VariableNames', VoltageTableHeaders);
    VoltageTable.Variables = VoltageTableData;
    app.VTable.Data = VoltageTable; % Update the app's voltage table
    app.VTable.ColumnName = VoltageTableHeaders; % Set column headers

    % Current Measurements Update
    % Initialize current table headers and variable types
    CurrentTableHeaders = {'Iname', 'CA', 'I (A)', 'Iul (A)', 'zeta (Vul)'};
    CurrentTableVarTypes = repmat({'string'}, 1, length(CurrentTableHeaders));
    CurrentTableData = {}; % Empty table for current data

    % Loop through each Control Area to update current measurements
    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i);

        % Loop through all current branches within the current control area
        for j = 1:length(CA.k_i)
            Iname = MATDSS.Sim.Meas.AllBranchesPhasesNames{CA.k_i(j)}; % Current branch name
            CAArea = CA.Area; % Control area identifier
            I = MATDSS.Meas.IProfile(CA.k_iIndex(j), Meast); % Current magnitude in A
            Iul = MATDSS.Meas.Imax(CA.k_iIndex(j)); % Upper limit for current
            zeta = CA.Duals.zeta(j, CA.TimeIndex); % Dual variable zeta

            % Store the current data in the table
            CurrentTableData = [CurrentTableData; {Iname, CAArea, I, Iul, zeta}];
        end
    end

    % Loop through all current branches to update general current measurements
    Imax = []; % Initialize Imax for general nodes
    for i = 1:length(MATDSS.Sim.Meas.Imax)
        Imax = [Imax; MATDSS.Sim.Meas.Imax{i}]; % Collect maximum current values
    end

    for i = 1:length(MATDSS.Sim.Meas.AllBranchesPhasesNames)
        Iname = MATDSS.Sim.Meas.AllBranchesPhasesNames{i}; % Current branch name
        CAArea = -1; % General branches not associated with a control area
        I = MATDSS.Sim.Meas.IProfile(i, t); % Current magnitude in A
        Iul = Imax(i); % Maximum current value
        zeta = -1; % Default value for general nodes

        % Store the current data in the table
        CurrentTableData = [CurrentTableData; {Iname, CAArea, I, Iul, zeta}];
    end

    % Create the current table and populate the GUI table component
    CurrentTable = table('Size', size(CurrentTableData), 'VariableTypes', CurrentTableVarTypes, 'VariableNames', CurrentTableHeaders);
    CurrentTable.Variables = CurrentTableData;
    app.ITable.Data = CurrentTable; % Update the app's current table
    app.ITable.ColumnName = CurrentTableHeaders; % Set column headers


    % Updating P0 measurements
    PTableHeaders = {'Pname', 'CA', 'P0 (kW)', 'Sum(P0) (kW)', 'P0_set (kW)', 'rho', 'lambda', 'mu'};
    PTableVarTypes = repmat({'string'}, 1, length(PTableHeaders)); % Define variable types for the table
    PTableData = {}; % Initialize empty cell array for storing data

    % Loop through each Control Area (CA) to update P0 measurements
    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i); % Get the current Control Area
        % Loop through all buses in the current Control Area
        for j = 1:length(CA.CABus0)
            Pname = MATDSS.Sim.Meas.AllNodesNames{CA.CABus0(j)}; % Get node name
            CAArea = CA.Area; % Get Control Area identifier
            P0 = (CA.P0(j, CA.TimeIndex) - CA.P0(j, 2)) / 1e3; % Active power in kW
            sumP0 = (sum(CA.P0(:, CA.TimeIndex)) - sum(CA.P0(:, 2))) / 1e3; % Total active power in kW
            P0Set = (CA.P0Set(CA.TimeIndex) - CA.P0Set(2)) / 1e3; % Set point for active power in kW
            rho = CA.Duals.rho(j, CA.TimeIndex); % Dual variable rho
            lambda = CA.Duals.lambda(CA.TimeIndex); % Dual variable lambda
            mu = CA.Duals.mu(CA.TimeIndex); % Dual variable mu

            % Append data to the PTableData cell array
            PTableData = [PTableData; {Pname, CAArea, P0, sumP0, P0Set, rho, lambda, mu}];
        end
    end

    % Create a table for P0 measurements and assign data
    PTable = table('Size', size(PTableData), 'VariableTypes', PTableVarTypes, 'VariableNames', PTableHeaders);
    PTable.Variables = PTableData; % Populate the table with collected data

    % Update the app's PTable with the newly created PTable
    app.PTable.Data = PTable;
    app.PTable.ColumnName = PTableHeaders; % Set column names for the app's table

    % Updating Q0 measurements
    QTableHeaders = {'Qname', 'CA', 'Q0 (kVar)', 'Sum(Q0) (kVar)', 'Q0_set (kVar)', 'sigma', 'eta', 'psi'};
    QTableVarTypes = repmat({'string'}, 1, length(QTableHeaders)); % Define variable types for the table
    QTableData = {}; % Initialize empty cell array for storing data

    % Loop through each Control Area to update Q0 measurements
    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i); % Get the current Control Area
        % Loop through all buses in the current Control Area
        for j = 1:length(CA.CABus0)
            Qname = MATDSS.Sim.Meas.AllNodesNames{CA.CABus0(j)}; % Get node name
            CAArea = CA.Area; % Get Control Area identifier
            Q0 = (CA.Q0(j, CA.TimeIndex) - CA.Q0(j, 2)) / 1e3; % Reactive power in kVar
            sumQ0 = (sum(CA.Q0(:, CA.TimeIndex)) - sum(CA.Q0(:, 2))) / 1e3; % Total reactive power in kVar
            Q0Set = (CA.Q0Set(CA.TimeIndex) - CA.Q0Set(2)) / 1e3; % Set point for reactive power in kVar
            sigma = CA.Duals.sigma(j, CA.TimeIndex); % Dual variable sigma
            eta = CA.Duals.eta(CA.TimeIndex); % Dual variable eta
            psi = CA.Duals.psi(CA.TimeIndex); % Dual variable psi

            % Append data to the QTableData cell array
            QTableData = [QTableData; {Qname, CAArea, Q0, sumQ0, Q0Set, sigma, eta, psi}];
        end
    end

    % Create a table for Q0 measurements and assign data
    QTable = table('Size', size(QTableData), 'VariableTypes', QTableVarTypes, 'VariableNames', QTableHeaders);
    QTable.Variables = QTableData; % Populate the table with collected data

    % Update the app's QTable with the newly created QTable
    app.QTable.Data = QTable;
    app.QTable.ColumnName = QTableHeaders; % Set column names for the app's table

end % End of live measurements function

pause(0.001); % Pause briefly to allow GUI updates

end

%% Old Function Code

%{

function MATDSSApp_LiveMeasurements(app,MATDSS,DER)
% this function will show live measurements and values for different
% parameters related to V, I, P and Q of all buses and DERs (focusing on
% the one under control).

t = MATDSS.Sim.t;
at = MATDSS.Sim.at;
if app.LiveMeasButton.Value
    % Updating voltage measurements
    VoltageTableHeaders = {'Vname','CA','Vll (kV)','V (kV)','Vul (kV)', 'gamma (Vul)','nu (Vll)'};
    VoltageTableVarTypes = repmat({'string'},1,length(VoltageTableHeaders));
    VoltageTableData = {};
    DelayIndex = ceil((MATDSS.Time.Meas.Delay/MATDSS.Time.Meas.TimeStep));
    pos = MATDSS.Meas.at - at;
    pos(pos > 0) = [];
    Meast = length(pos);
    Meast = Meast - DelayIndex;
    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i);
        for j = 1:length(CA.k_v)
            Vname = MATDSS.Sim.Meas.AllNodesNames{CA.k_v(j)};
            CAArea = CA.Area;
            VBase = MATDSS.Meas.VBases(CA.k_v(j));
            Vll = CA.vll.*VBase;
            V = MATDSS.Meas.VMagProfile(CA.k_vIndex(j),Meast)./1e3;
            Vul = CA.vul.*VBase;
            gamma = CA.Duals.gamma(j,CA.TimeIndex);
            nu = CA.Duals.nu(j,CA.TimeIndex);
            VoltageTableData = [VoltageTableData;...
                {Vname,CAArea,Vll,V,Vul,gamma,nu}];
        end
    end

    for i = 1:length(MATDSS.Sim.Meas.AllNodesNames)
        Vname = MATDSS.Sim.Meas.AllNodesNames{i};
        CAArea = -1;
        VBase = MATDSS.Meas.VBases(i);
        Vll = CA.vll.*VBase;
        V = MATDSS.Sim.Meas.VMagProfile(i,t)./1e3;
        Vul = CA.vul.*VBase;
        gamma = -1;
        nu = -1;
        VoltageTableData = [VoltageTableData;...
            {Vname,CAArea,Vll,V,Vul,gamma,nu}];

    end

    VoltageTable = table('Size',size(VoltageTableData),'VariableTypes',VoltageTableVarTypes,'VariableNames',VoltageTableHeaders);
    VoltageTable.Variables = VoltageTableData;

    app.VTable.Data = VoltageTable;
    app.VTable.ColumnName = VoltageTableHeaders;




    % Updating current measurements
    CurrentTableHeaders = {'Iname','CA','I (A)', 'Iul (A)', 'zeta (Vul)'};
    CurrentTableVarTypes = repmat({'string'},1,length(CurrentTableHeaders));
    CurrentTableData = {};

    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i);
        for j = 1:length(CA.k_i)
            Iname = MATDSS.Sim.Meas.AllBranchesPhasesNames{CA.k_i(j)};
            CAArea = CA.Area;
            I = MATDSS.Meas.IProfile(CA.k_iIndex(j),Meast);
            Iul = MATDSS.Meas.Imax(CA.k_iIndex(j));
            zeta = CA.Duals.zeta(j,CA.TimeIndex);
            CurrentTableData = [CurrentTableData;...
                {Iname,CAArea,I,Iul,zeta}];
        end
    end
    Imax = [];
    for i = 1:length(MATDSS.Sim.Meas.Imax)
        Imax = [Imax; MATDSS.Sim.Meas.Imax{i}];
    end
    for i = 1:length(MATDSS.Sim.Meas.AllBranchesPhasesNames)
        Iname = MATDSS.Sim.Meas.AllBranchesPhasesNames{i};
        CAArea = -1;
        I = MATDSS.Sim.Meas.IProfile(i,t);
        Iul = Imax(i);
        zeta = -1;
        CurrentTableData = [CurrentTableData;...
            {Iname,CAArea,I,Iul,zeta}];

    end

    CurrentTable = table('Size',size(CurrentTableData),'VariableTypes',CurrentTableVarTypes,'VariableNames',CurrentTableHeaders);
    CurrentTable.Variables = CurrentTableData;

    app.ITable.Data = CurrentTable;
    app.ITable.ColumnName = CurrentTableHeaders;



    % Updating P0 measurements
    PTableHeaders = {'Pname','CA','P0 (kW)', 'Sum(P0) (kW)', 'P0_set (kW)','rho','lambda','mu'};
    PTableVarTypes = repmat({'string'},1,length(PTableHeaders));
    PTableData = {};

    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i);
        for j = 1:length(CA.CABus0)
            Pname = MATDSS.Sim.Meas.AllNodesNames{CA.CABus0(j)};
            CAArea = CA.Area;
            P0 = (CA.P0(j,CA.TimeIndex)-CA.P0(j,2))./1e3;
            sumP0 = (sum(CA.P0(:,CA.TimeIndex))-sum(CA.P0(:,2)))./1e3;
            P0Set = (CA.P0Set(CA.TimeIndex)-CA.P0Set(2))./1e3;
            rho = CA.Duals.rho(j,CA.TimeIndex);
            lambda = CA.Duals.lambda(CA.TimeIndex);
            mu = CA.Duals.mu(CA.TimeIndex);
            PTableData = [PTableData;...
                {Pname,CAArea,P0,sumP0,P0Set,rho,lambda,mu}];
        end
    end

    PTable = table('Size',size(PTableData),'VariableTypes',PTableVarTypes,'VariableNames',PTableHeaders);
    PTable.Variables = PTableData;

    app.PTable.Data = PTable;
    app.PTable.ColumnName = PTableHeaders;



    % Updating Q0 measurements
    QTableHeaders = {'Qname','CA','Q0 (kVar)', 'Sum(Q0) (kVar)', 'Q0_set (kVar)','sigma','eta','psi'};
    QTableVarTypes = repmat({'string'},1,length(QTableHeaders));
    QTableData = {};

    for i = 1:length(MATDSS.Cont.CA)
        CA = MATDSS.Cont.CA(i);
        for j = 1:length(CA.CABus0)
            Qname = MATDSS.Sim.Meas.AllNodesNames{CA.CABus0(j)};
            CAArea = CA.Area;
            Q0 = (CA.Q0(j,CA.TimeIndex) - CA.Q0(j,2))./1e3;
            sumQ0 = (sum(CA.Q0(:,CA.TimeIndex)) - sum(CA.Q0(:,2)))./1e3;
            Q0Set = (CA.Q0Set(CA.TimeIndex)-CA.Q0Set(2))./1e3;
            sigma = CA.Duals.sigma(j,CA.TimeIndex);
            eta = CA.Duals.eta(CA.TimeIndex);
            psi = CA.Duals.psi(CA.TimeIndex);
            QTableData = [QTableData;...
                {Qname,CAArea,Q0,sumQ0,Q0Set,sigma,eta,psi}];
        end
    end

    QTable = table('Size',size(QTableData),'VariableTypes',QTableVarTypes,'VariableNames',QTableHeaders);
    QTable.Variables = QTableData;

    app.QTable.Data = QTable;
    app.QTable.ColumnName = QTableHeaders;


end


pause(0.001);

end



%}