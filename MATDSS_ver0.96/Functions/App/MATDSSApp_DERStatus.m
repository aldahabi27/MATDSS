function MATDSSApp_DERStatus(app, MATDSS, DER)
% MATDSSApp_DERStatus updates the status of Distributed Energy Resources (DER)
% in the MATDSS application.
%
% This function retrieves the current status of each DER and displays it 
% in a table format within the application interface. The data includes 
% real power (P), reactive power (Q), and their set points, as well as 
% changes over the current simulation step.
%
% The function only executes updates if the Live DER button is active and 
% if the DER's UpdateStatus flag is set to true. It populates the DER 
% status table with the appropriate data and updates the column names.
%
% Parameters:
%   - app: The MATDSS application instance containing properties and 
%          methods related to the simulation and GUI management.
%   - MATDSS: A structure containing all necessary data and settings 
%             for the simulation, including time parameters and control 
%             signals.
%   - DER: An array of structures containing information about each 
%          Distributed Energy Resource, including their statuses and 
%          measurements.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

if DER(1).UpdateStatus && app.LiveDERButton.Value
    nDER = MATDSS.Sim.nDER; % Number of DERs in the simulation
    t = MATDSS.Sim.t; % Current simulation time index

    DERTableData = {}; % Initialize cell array for DER table data
    % Define the headers for the DER status table
    DERTableHeaders = {'#', 'DER', 'P (kW)', 'P_set (kW)', 'Q (kVar)', 'Q_set (kVar)', [char(916) 'P (kW)'], [char(916) 'Q (kVar)']};
    DERTableVarTypes = repmat({'string'}, 1, length(DERTableHeaders)); % Variable types for table columns

    % Loop through each DER to populate the table data
    for i = 1:size(DER, 2)
        i_WVar = DER(i).i_WVar; % Index for W variable
        i_setpoint = DER(i).i_setpoint; % Index for setpoint variable
        
        if t == 1 % At the first time step, initialize table data
            DERTableData = [DERTableData; {DER(i).DER_Ind, DER(i).DSSName, ...
                DER(i).W(i_WVar, 2)./1e3, DER(i).WSetpoint(i_setpoint, 2)./1e3, ...
                DER(i).Var(i_WVar, 2)./1e3, DER(i).VarSetpoint(i_setpoint, 2)./1e3, 0, 0}];
            app.DERStatusTable.ColumnName = DERTableHeaders; % Set column names for the table
        else
            % Handle the specific case for 'vder' type DER
            if strcmpi(DER(i).Type, 'vder')
                Contt_window = MATDSS.Time.Cont.TimeStep / MATDSS.Time.Sim.TimeStep; % Control window for vder
                DERTableData = [DERTableData; {DER(i).DER_Ind, DER(i).DSSName, ...
                    (DER(i).W(i_WVar, 2) - DER(i).W(1, 2))./1e3, ...
                    (DER(i).WSetpoint(i_setpoint, 2) - DER(i).WSetpoint(1, 2))./1e3, ...
                    (DER(i).Var(i_WVar, 2) - DER(i).Var(1, 2))./1e3, ...
                    (DER(i).VarSetpoint(i_setpoint, 2) - DER(i).VarSetpoint(1, 2))./1e3, ...
                    (DER(i).W(i_WVar, 2) - DER(i).W(i_WVar - Contt_window, 2))./1e3, ...
                    (DER(i).Var(i_WVar, 2) - DER(i).Var(i_WVar - Contt_window, 2))./1e3}];
            else
                % For all other DER types, calculate the differences for the table
                DERTableData = [DERTableData; {DER(i).DER_Ind, DER(i).DSSName, ...
                    (DER(i).W(i_WVar, 2) - DER(i).W(1, 2))./1e3, ...
                    (DER(i).WSetpoint(i_setpoint, 2) - DER(i).WSetpoint(1, 2))./1e3, ...
                    (DER(i).Var(i_WVar, 2) - DER(i).Var(1, 2))./1e3, ...
                    (DER(i).VarSetpoint(i_setpoint, 2) - DER(i).VarSetpoint(1, 2))./1e3, ...
                    (DER(i).W(i_WVar, 2) - DER(i).W(i_WVar - 1, 2))./1e3, ...
                    (DER(i).Var(i_WVar, 2) - DER(i).Var(i_WVar - 1, 2))./1e3}];
            end            
        end
    end

    % Create a table with the updated data for the DER status
    DERTable = table('Size', [nDER, length(DERTableHeaders)], 'VariableTypes', DERTableVarTypes, ...
        'VariableNames', DERTableHeaders);
    DERTable.Variables = DERTableData; % Assign the data to the table
    app.DERStatusTable.Data = DERTable; % Update the application table with new data
    app.DERStatusTable.ColumnName = DERTableHeaders; % Set the column names for the table
end


%% Old function code

%{
function MATDSSApp_DERStatus(app,MATDSS,DER)


if DER(1).UpdateStatus && app.LiveDERButton.Value
    nDER = MATDSS.Sim.nDER;
    t = MATDSS.Sim.t;

    DERTableData = {};
    DERTableHeaders = {'#', 'DER', 'P (kW)', 'P_set (kW)', 'Q (kVar)', 'Q_set (kVar)', [char(916) 'P (kW)'], [char(916) 'Q (kVar)']};
%     DERTableHeaders = {'#', 'DER', 'P (kW)', 'Q (kVar)', [char(916) 'P (kW)'], [char(916) 'Q (kVar)']};
    DERTableVarTypes = repmat({'string'},1,length(DERTableHeaders));

    for i = 1:size(DER,2)
        i_WVar = DER(i).i_WVar;
        i_setpoint = DER(i).i_setpoint;
        if t == 1
            DERTableData = [DERTableData; {DER(i).DER_Ind, DER(i).DSSName, DER(i).W(i_WVar,2)./1e3, DER(i).WSetpoint(i_setpoint,2)./1e3, DER(i).Var(i_WVar,2)./1e3, DER(i).VarSetpoint(i_setpoint,2)./1e3, 0,0}];
            app.DERStatusTable.ColumnName = DERTableHeaders;
        else
            if strcmpi(DER(i).Type, 'vder')
                Contt_window = MATDSS.Time.Cont.TimeStep/MATDSS.Time.Sim.TimeStep;
                DERTableData = [DERTableData; {DER(i).DER_Ind, DER(i).DSSName, (DER(i).W(i_WVar,2)-DER(i).W(1,2))./1e3, (DER(i).WSetpoint(i_setpoint,2)-DER(i).WSetpoint(1,2))./1e3, (DER(i).Var(i_WVar,2)-DER(i).Var(1,2))./1e3, (DER(i).VarSetpoint(i_setpoint,2)-DER(i).VarSetpoint(1,2))./1e3, (DER(i).W(i_WVar,2)-DER(i).W(i_WVar-Contt_window,2))./1e3,(DER(i).Var(i_WVar,2)-DER(i).Var(i_WVar-Contt_window,2))./1e3}];
            else
                DERTableData = [DERTableData; {DER(i).DER_Ind, DER(i).DSSName, (DER(i).W(i_WVar,2)-DER(i).W(1,2))./1e3, (DER(i).WSetpoint(i_setpoint,2)-DER(i).WSetpoint(1,2))./1e3, (DER(i).Var(i_WVar,2)-DER(i).Var(1,2))./1e3, (DER(i).VarSetpoint(i_setpoint,2)-DER(i).VarSetpoint(1,2))./1e3, (DER(i).W(i_WVar,2)-DER(i).W(i_WVar-1,2))./1e3,(DER(i).Var(i_WVar,2)-DER(i).Var(i_WVar-1,2))./1e3}];
            end            
        end
    end

    DERTable = table('Size',[nDER,length(DERTableHeaders)],'VariableTypes',DERTableVarTypes,'VariableNames',DERTableHeaders);
    DERTable.Variables = DERTableData;
    app.DERStatusTable.Data = DERTable;
    app.DERStatusTable.ColumnName = DERTableHeaders;
end

%}