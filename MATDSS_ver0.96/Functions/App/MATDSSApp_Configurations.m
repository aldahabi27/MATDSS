% MATDSSApp_Configurations(app)
% This function loads the configurations saved in Excel files for the 
% MATDSS application. These configurations are used for simulation settings 
% such as DERs but are not directly used to generate Plot properties by default.
%
% If the necessary configuration files are not found, it creates them and 
% initializes the corresponding settings with default values.
%
% Last Update for this function was on MATDSS App Ver 0.96 (19 Sept. 2024)
%
%   MATDSS Application
%   Copyright (c) 2024, Ilyas Farhat
%   by Ilyas Farhat
%
%   This file is part of MATDSS Application
%   Contact the developer at ilyas.farhat@outlook.com

function MATDSSApp_Configurations(app)
    % This function loads the configurations saved in Excel files to MATDSS
    % Application. Those configurations are not used to generate Plot
    % properties panel by default!

    MATDSSApp_Status(app, 'Loading Configuration Files', ...
        ['Loading configuration files for "' app.OpenDSSFilesListBox.Value '"']);

    % Tag all tabs as 'Not Loaded'
    for i = 1:numel(app.CircuitConfigurationsTabGroup.Children)
        app.CircuitConfigurationsTabGroup.Children(i).Tag = 'Not Loaded';
    end

    % Set the target sheet in the Excel file
    app.MATDSSExcelSheet = app.OpenDSSFilesListBox.Value;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DER Tab Setup                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MATDSSApp_Details(app, 'Loading DER table');
    DERTableHeaders = {'Index', 'Name', 'Bus (# or ''Bus Name'')', 'Tau', 'DER Type', 'Nodes', ...
                       'Connection Type', 'Nphase', 'Mode', 'P(x)', 'Q(x)', 'Pmin', ...
                       'Pmax', 'Qmin', 'Qmax', 'ax', 'cx'};
    DERTableVarTypes = repmat({'string'}, 1, length(DERTableHeaders));

    % Create the MATDSS_DER.xlsx file if it doesn't exist
    derFilePath = fullfile(pwd, 'Config', 'MATDSS_DER.xlsx');
    if ~exist(derFilePath, 'file')
        DERTable = table('Size', [0, length(DERTableHeaders)], ...
                         'VariableTypes', DERTableVarTypes, ...
                         'VariableNames', DERTableHeaders);
        writetable(DERTable, derFilePath, 'Sheet', app.MATDSSExcelSheet, ...
                   'WriteMode', 'overwritesheet', 'autofitwidth', 1);
    end

    % Check if the target sheet exists in the Excel file
    MATDSSDERSheetNames = cellstr(sheetnames(derFilePath)); % List all sheets in the Excel file
    MATDSSDERFlag = MATDSS_StrComp(MATDSSDERSheetNames, app.MATDSSExcelSheet); % Check for the target sheet

    % Create a new empty table if the target sheet is not found
    if MATDSSDERFlag <= 0
        DERTable = table('Size', [0, length(DERTableHeaders)], ...
                         'VariableTypes', DERTableVarTypes, ...
                         'VariableNames', DERTableHeaders);
        writetable(DERTable, derFilePath, 'Sheet', app.MATDSSExcelSheet, ...
                   'WriteMode', 'overwritesheet', 'autofitwidth', 1);
    end

    % Set import options for reading the table from the Excel sheet
    DERReadOptions = detectImportOptions(derFilePath, 'Sheet', app.MATDSSExcelSheet, ...
                                         'VariableNamingRule', 'preserve');
    DERReadOptions = setvartype(DERReadOptions, DERReadOptions.VariableNames(:), 'string');

    % Read DER table from Excel
    DERTable = readtable(derFilePath, DERReadOptions);
    DERTableContent = DERTable.Variables;

    % Default values for missing entries in the DER table
    DERTableDefaultValues = ["";"DER_";"";"0.2";"PV";"[1,2,3]";"Wye";"3";"DSS_load";"[2,0,0]"; ...
                             "[2,0,0]";"-1e6";"1e6";"-1e6";"1e6";"1";"1"];

    % Fill missing entries in the table
    for i = 1:size(DERTableContent, 2)
        if i == 1
            DERTableContent(find(ismissing(DERTableContent(:, i))), i) = ...
                num2str(find(ismissing(DERTableContent(:, i))));
        elseif i == 2
            DERTableContent(find(ismissing(DERTableContent(:, i))), i) = ...
                strcat(DERTableDefaultValues(i), num2str(find(ismissing(DERTableContent(:, i)))));
        else
            DERTableContent(find(ismissing(DERTableContent(:, i))), i) = DERTableDefaultValues(i);
        end
    end

    % Update table data and properties
    DERTable.Variables = DERTableContent;
    app.DERsTable.Data = DERTable;
    app.DERsTable.ColumnName = DERTable.Properties.VariableNames;
    app.DERsTable.ColumnWidth = 'auto';
    app.DERsTab.Tag = 'Loaded';

    % Run MATDSS quick simulation to generate list of Buses, Phases, and Branches
    [MATDSS] = MATDSS_Initialize(app, 1);
    MATDSS.Table.DER = DERTable;
    MATDSS.TableData.DER = DERTableContent;
    if MATDSS.Sim.DSSStartOk % If DSS is loaded and linked successfully
        % Get all node names and bus names
        AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
        AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames; 
        AllBusNames(1) = [];

        % Get all lines/branches names
        MyLines = MATDSS.Sim.DSSCircuit.Lines;
        MyLinesAllNames = MyLines.AllNames;
        MyLinesAllPhasesNames = {};
        for i = 1:numel(MyLinesAllNames)
            MyLines.Name = MyLinesAllNames{i};
            nphases = MyLines.Phases;
            for j = 1:nphases
                MyLinesAllPhasesNames = [MyLinesAllPhasesNames; [MyLinesAllNames{i} '.' num2str(j)]];
            end
        end
    else
        AllNodesNames = {'Error! Could not initiate the link with OpenDSS'};
        MyLinesAllPhasesNames = {'Error! Could not initiate the link with OpenDSS'};
    end



% This code was deprecated with the introduction of multi-control area 
% structure. It is kept here for reference only. In future releases, the
% code will be refined and this will be removed.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 L2C Tab Setup                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MATDSSL2CExcelSheetMissingFlag = true; % Flag indicating if L2C Excel sheet is missing
if exist([pwd '\Config\MATDSS_L2C.xlsx'], 'file') % Check if the L2C Excel file exists
    MATDSSL2CSheetnames = cellstr(sheetnames('MATDSS_L2C.xlsx')); % List all sheets in the Excel file
    MATDSSL2CFlag = MATDSS_StrComp(MATDSSL2CSheetnames, app.MATDSSExcelSheet); % Check if the target sheet exists
    if MATDSSL2CFlag > 0 % If sheet exists
        MATDSSL2CExcelSheetMissingFlag = false; % Mark the Excel sheet as found
        % Configure import options
        L2CReadOptions = detectImportOptions('MATDSS_L2C.xlsx', 'Sheet', app.MATDSSExcelSheet, 'VariableNamingRule', 'preserve');
        L2CReadOptions = setvartype(L2CReadOptions, L2CReadOptions.VariableNames(:), 'string'); % Set all variable types to string
        L2CTable = readtable([pwd '\Config\MATDSS_L2C.xlsx'], L2CReadOptions); % Read table data from Excel sheet
        L2CTableData = L2CTable.Variables; % Get table data
        
        % Assign values to respective edit fields
        app.alphaEditField.Value = L2CTableData(1);
        app.rpEditField.Value = L2CTableData(2);
        app.rbardEditField.Value = L2CTableData(3);
        app.EEditField.Value = L2CTableData(4);
        app.vulEditField.Value = L2CTableData(5);
        app.vllEditField.Value = L2CTableData(6);
        app.iulEditField.Value = L2CTableData(7);
        app.a_rhoEditField.Value = L2CTableData(8);
        app.a_lambdaEditField.Value = L2CTableData(9);
        app.a_muEditField.Value = L2CTableData(10);
        app.a_gammaEditField.Value = L2CTableData(11);
        app.a_nuEditField.Value = L2CTableData(12);
        app.a_zetaEditField.Value = L2CTableData(13);
        app.c_rhoEditField.Value = L2CTableData(14);
        app.c_lambdaEditField.Value = L2CTableData(15);
        app.c_muEditField.Value = L2CTableData(16);
        app.c_gammaEditField.Value = L2CTableData(17);
        app.c_nuEditField.Value = L2CTableData(18);
        app.c_zetaEditField.Value = L2CTableData(19);
    end
end

% Handle missing Excel sheet by creating a new one with default values
if MATDSSL2CExcelSheetMissingFlag
    L2CTableHeaders = {'alpha', 'rp', 'rbard', 'E (W)', 'vul (p.u.)', 'vll (p.u.)', 'iul (p.u.)', 'a_rho', 'a_lambda', 'a_mu', 'a_gamma (W^2/V^2)', 'a_nu (W^2/V^2)', 'a_zeta (W^2/A^2)', 'c_rho', 'c_lambda', 'c_mu', 'c_gamma (V^2/W^2)', 'c_nu (V^2/W^2)', 'c_zeta (A^2/W^2)'};
    L2CTableVarTypes = repmat({'string'}, 1, length(L2CTableHeaders)); % Set variable types as string
    L2CTableSize = [1, length(L2CTableHeaders)];
    L2CTableData = repmat("", L2CTableSize); % Initialize empty table data
    % Assign current app values to table data
    L2CTableData(1,1) = app.alphaEditField.Value;
    L2CTableData(1,2) = app.rpEditField.Value;
    L2CTableData(1,3) = app.rbardEditField.Value;
    L2CTableData(1,4) = app.EEditField.Value;
    L2CTableData(1,5) = app.vulEditField.Value;
    L2CTableData(1,6) = app.vllEditField.Value;
    L2CTableData(1,7) = app.iulEditField.Value;
    L2CTableData(1,8) = app.a_rhoEditField.Value;
    L2CTableData(1,9) = app.a_lambdaEditField.Value;
    L2CTableData(1,10) = app.a_muEditField.Value;
    L2CTableData(1,11) = app.a_gammaEditField.Value;
    L2CTableData(1,12) = app.a_nuEditField.Value;
    L2CTableData(1,13) = app.a_zetaEditField.Value;
    L2CTableData(1,14) = app.c_rhoEditField.Value;
    L2CTableData(1,15) = app.c_lambdaEditField.Value;
    L2CTableData(1,16) = app.c_muEditField.Value;
    L2CTableData(1,17) = app.c_gammaEditField.Value;
    L2CTableData(1,18) = app.c_nuEditField.Value;
    L2CTableData(1,19) = app.c_zetaEditField.Value;
    
    L2CTable = cell2table(cellstr(L2CTableData)); % Convert cell array to table
    L2CTable.Properties.VariableNames = L2CTableHeaders; % Set table headers
    writetable(L2CTable, [pwd '\Config\MATDSS_L2C.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1); % Write table to Excel
end
app.L2CTab.Tag = 'Loaded'; % Mark the tab as loaded


%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               V & I Tab Setup                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSSApp_Details(app, 'Loading V & I tables'); % Display loading message

% Populate the listboxes with data from OpenDSS
if size(AllNodesNames, 1) == 1 % If error in OpenDSS link
    app.VTrackingListBox.Items = AllNodesNames; % Set listbox items to error message
    app.ITrackingListBox.Items = MyLinesAllPhasesNames;
else
    app.VTrackingListBox.Items = AllNodesNames(3+1:end); % Set listbox items for voltage
    app.ITrackingListBox.Items = MyLinesAllPhasesNames; % Set listbox items for current
    
    if exist([pwd '\Config\MATDSS_VI.xlsx'], 'file') % Check if VI Excel file exists
        MATDSSVISheetnames = cellstr(sheetnames('MATDSS_VI.xlsx')); % List all sheets in the Excel file
        MATDSSVIFlag = MATDSS_StrComp(MATDSSVISheetnames, app.MATDSSExcelSheet); % Check if the target sheet exists
        
        if MATDSSVIFlag > 0 % If the sheet exists
            % Configure import options
            VIReadOptions = detectImportOptions('MATDSS_VI.xlsx', 'Sheet', app.MATDSSExcelSheet, 'VariableNamingRule', 'preserve');
            VIReadOptions = setvartype(VIReadOptions, VIReadOptions.VariableNames(:), 'string'); % Set all variable types to string
            VITable = readtable([pwd '\Config\MATDSS_VI.xlsx'], VIReadOptions); % Read table data from Excel sheet
            VITableData = VITable.Variables; % Get table data

            % Populate VTrackingListBox with saved values
            Mv = VITableData(:, 1); 
            Mv(ismissing(Mv)) = []; % Remove missing values
            v_size = length(Mv);
            A = app.VTrackingListBox.Items'; % Get current listbox items
            B = convertStringsToChars(VITableData(1:v_size, 1)); % Convert table data to characters
            [C, iA, iB] = intersect(A, B, 'stable'); % Find intersection between saved and current items
            k = iA;
            if k >= 1 % If matching item found, set the listbox value
                app.VTrackingListBox.Value = app.VTrackingListBox.Items(k);
            end

            % Populate ITrackingListBox with saved values
            Mi = VITableData(:, 3);
            Mi(ismissing(Mi)) = [];
            i_size = length(Mi);
            A = app.ITrackingListBox.Items';
            B = convertStringsToChars(VITableData(1:i_size, 3));
            [C, iA, iB] = intersect(A, B, 'stable');
            k = iA;
            if k >= 1 % If matching item found, set the listbox value
                app.ITrackingListBox.Value = app.ITrackingListBox.Items(k);
            end
        end
    end
end
app.VITab.Tag = 'Loaded'; % Mark the tab as loaded


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Control Areas Tab Setup             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSSApp_Details(app, 'Loading controllers tables'); % Status message for loading control areas
MATDSSControlAreasExcelSheetMissingFlag = true; % Initialize the flag for missing control areas Excel sheet
if exist([pwd '\Config\MATDSS_ControlAreas.xlsx'], 'file') % Check if the Control Areas Excel file exists
    MATDSSControlAreasSheetnames = cellstr(sheetnames('MATDSS_ControlAreas.xlsx')); % List all sheets in the Excel file
    MATDSSControlAreasFlag = MATDSS_StrComp(MATDSSControlAreasSheetnames, app.MATDSSExcelSheet); % Check if the target sheet exists
    if MATDSSControlAreasFlag > 0 % If the target sheet is available
        MATDSSControlAreasExcelSheetMissingFlag = false; % Update the missing sheet flag
        ControlAreasReadOptions = detectImportOptions('MATDSS_ControlAreas.xlsx', 'Sheet', app.MATDSSExcelSheet, 'VariableNamingRule', 'preserve'); % Set import options
        ControlAreasReadOptions = setvartype(ControlAreasReadOptions, ControlAreasReadOptions.VariableNames(:), 'string'); % Force variable types to string
        clear ControlAreasTable
        ControlAreasTable = readtable('MATDSS_ControlAreas.xlsx', ControlAreasReadOptions); % Read the table from the Excel sheet
        ControlAreasTableContent = ControlAreasTable.Variables; % Extract table content
        ControlAreasTableDefaultValues = ["";"";"1";"F";""]; % Define default values for missing entries
        for i = 1:size(ControlAreasTable.Variables, 2)
            ControlAreasTableContent(find(ismissing(ControlAreasTableContent(:, i))), i) = ControlAreasTableDefaultValues(i); % Replace missing values with defaults
        end

        ControlAreasTable.Variables = ControlAreasTableContent; % Update table with corrected content
        app.ControlAreasTable.Data = ControlAreasTable; % Assign data to the app's UI table
        app.ControlAreasTable.ColumnName = ControlAreasTable.Properties.VariableNames; % Set column names in the UI table
        app.ControlAreasTable.ColumnWidth = 'auto'; % Auto-adjust column width
    end
end

if MATDSSControlAreasExcelSheetMissingFlag % If the Control Areas Excel sheet is missing
    AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames; % Get all bus names from the simulation

    % Control Areas Settings Section
    ControlAreasTableHeaders = {'Bus Name', 'Bus Number', 'Assigned Control Area (#)', 'Interface Bus (T/F)', 'Connected to Bus name'}; % Define table headers
    ControlAreasTableVarTypes = repmat({'string'}, 1, length(ControlAreasTableHeaders)); % Set variable types as strings for all columns
    ControlAreasTableSize = [length(AllBusNames), length(ControlAreasTableHeaders)]; % Define the table size based on the number of buses and headers
    ControlAreasTableData = repmat("", ControlAreasTableSize); % Initialize table data with empty strings
    ControlAreasTableData(:, 1) = AllBusNames; % Set bus names in the first column
    ControlAreasTableData(:, 2) = num2str([1:length(AllBusNames)]'); % Set bus numbers in the second column
    ControlAreasTableData(:, 3) = num2str(1); % Default all control areas to 1
    ControlAreasTableData(:, 4) = 'F'; % Default interface bus column to 'F' (false)
    ControlAreasTableData(:, 5) = ''; % Empty connected bus names by default
    ControlAreasTable = cell2table(cellstr(ControlAreasTableData)); % Convert to a table structure
    ControlAreasTable.Properties.VariableNames = ControlAreasTableHeaders; % Set table headers
    writetable(ControlAreasTable, [pwd '\Config\MATDSS_ControlAreas.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1); % Write table to Excel file
    app.ControlAreasTable.Data = ControlAreasTable; % Assign data to the app's UI table
    app.ControlAreasTable.ColumnName = ControlAreasTable.Properties.VariableNames; % Set column names in the UI table
    app.ControlAreasTable.ColumnWidth = 'auto'; % Auto-adjust column width
    app.ControlAreasTable.Tag = 'Loaded'; % Mark the table as loaded
end

app.ControlAreasTableData = ControlAreasTable; % Store the control areas table data in the app's properties

% Below are commented-out sections that were previously used for setting
% up the Control Areas table, but replaced with the current implementation
% using the loaded or default Excel sheet.

%             if MATDSSControlAreaFlag <= 0 % if Sheet is not available
%                 ControlAreasTable = table('Size', [0, length(ControlAreasTableHeaders)], 'VariableTypes', ControlAreasTableVarTypes, 'VariableNames', ControlAreasTableHeaders);
%                 writetable(ControlAreasTable,[pwd '\Config\MATDSS_ControlAreas.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
%             end
%
%             ControlAreasReadOptions = detectImportOptions('MATDSS_ControlAreas.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
%             ControlAreasReadOptions = setvartype(ControlAreasReadOptions,ControlAreasReadOptions.VariableNames(:),'string');
%             clear ControlAreasTable
%             ControlAreasTable = readtable('MATDSS_ControlAreas.xlsx',ControlAreasReadOptions);
%             ControlAreasTableContent = ControlAreasTable.Variables;
%             ControlAreasTableDefaultValues = ["";"1"];
%             for i = 1:size(ControlAreasTable.Variables,2)
%                 ControlAreasTableContent(find(ismissing(ControlAreasTableContent(:,i))),i) = ControlAreasTableDefaultValues(i);
%             end
%
%             ControlAreasTable.Variables = ControlAreasTableContent;
%             app.ControlAreasTable.Data = ControlAreasTable;
%             app.ControlAreasTable.ColumnName = ControlAreasTable.Properties.VariableNames;
%             app.ControlAreasTable.ColumnWidth = 'auto';
%             app.ControlAreasTable.Tag = 'Loaded';

% Setting default values for Controllers Settings
ControllersSettingsTableDefaultValues = ["1"; "LLC"; "auto"; "1e-4"; "1e-3"; "1e2"; "1.05"; "0.95"; "Specified in DSS File";...
    "1"; "1"; "1e3"; "1e3"; "1e3"; "1e3"; "1e12"; "1e12"; "1e8";...
    "1"; "1"; "1e-3"; "1e-3"; "1e-3"; "1e-3"; "1e-12"; "1e-12"; "1e-8"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Controllers Settings Section        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSSControllersSettingsExcelSheetMissingFlag = true; % Initialize flag for missing controllers settings Excel sheet
if exist([pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'file') % Check if the Controllers Settings Excel file exists
    MATDSSControllersSettingsSheetnames = cellstr(sheetnames('MATDSS_ControllersSettings.xlsx')); % List all sheets in the Excel file
    MATDSSControllersSettingsFlag = MATDSS_StrComp(MATDSSControllersSettingsSheetnames, app.MATDSSExcelSheet); % Check if the target sheet exists
    if MATDSSControllersSettingsFlag > 0 % If the target sheet is available
        MATDSSControllersSettingsExcelSheetMissingFlag = false; % Update the missing sheet flag
        ControllersSettingsReadOptions = detectImportOptions('MATDSS_ControllersSettings.xlsx', 'Sheet', app.MATDSSExcelSheet, 'VariableNamingRule', 'preserve'); % Set import options
        ControllersSettingsReadOptions = setvartype(ControllersSettingsReadOptions, ControllersSettingsReadOptions.VariableNames(:), 'string'); % Force variable types to string
        clear ControllersSettingsTable
        ControllersSettingsTable = readtable('MATDSS_ControllersSettings.xlsx', ControllersSettingsReadOptions); % Read the table from the Excel sheet
        ControllersSettingsTableContent = ControllersSettingsTable.Variables; % Extract table content
        for i = 1:size(ControllersSettingsTable.Variables, 2)
            ControllersSettingsTableContent(find(ismissing(ControllersSettingsTableContent(:, i))), i) = ControllersSettingsTableDefaultValues(i); % Replace missing values with defaults
        end

        ControllersSettingsTable.Variables = ControllersSettingsTableContent; % Update table with corrected content
        app.ControllersSettingsTable.Data = ControllersSettingsTable; % Assign data to the app's UI table
        app.ControllersSettingsTable.ColumnName = ControllersSettingsTable.Properties.VariableNames; % Set column names in the UI table
        app.ControllersSettingsTable.ColumnWidth = 'auto'; % Auto-adjust column width
    end
end



if MATDSSControllersSettingsExcelSheetMissingFlag
    % If the controllers settings sheet is missing, create and populate a new one
    ControllersSettingsTableHeaders = {'# Control Area','Controller Type','alpha','r_p','rbar_d','E (W)','v_ul (p.u.)','v_ll (p.u.)','i_ul (p.u.)',...
        'a_rho','a_sigma','a_lambda','a_mu','a_eta','a_psi','a_gamma (W^2/V^2)','a_nu (W^2/V^2)','a_zeta (W^2/A^2)',...
        'c_rho','c_sigma','c_lambda','c_mu','c_eta','c_psi','c_gamma (V^2/W^2)','c_nu (V^2/W^2)','c_zeta (A^2/W^2)'};
    ControllersSettingsTableSize = [1,length(ControllersSettingsTableHeaders)];  % Define size of the table (1 row for defaults)
    ControllersSettingsTableData = repmat("",ControllersSettingsTableSize);  % Initialize with empty strings
    ControllersSettingsTableData(1,:) = ControllersSettingsTableDefaultValues;  % Populate with default values
    ControllersSettingsTable = cell2table(cellstr(ControllersSettingsTableData));  % Convert to cell table format
    ControllersSettingsTable.Properties.VariableNames = ControllersSettingsTableHeaders;  % Assign headers to the table
    writetable(ControllersSettingsTable,[pwd '\Config\MATDSS_ControllersSettings.xlsx'],... % Write table to Excel file
        'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
    app.ControllersSettingsTable.Data = ControllersSettingsTable;  % Update app table with new data
    app.ControllersSettingsTable.ColumnName = ControllersSettingsTable.Properties.VariableNames;  % Set table column names
    app.ControllersSettingsTable.ColumnWidth = 'auto';  % Auto-adjust column widths
    app.ControllersSettingsTable.Tag = 'Loaded';  % Mark the table as loaded
end

% VDER Configs (Voltage-Derived Emission Reduction) Configuration
MATDSSVDERsExcelSheetMissingFlag = true;
if exist([pwd '\Config\MATDSS_VDERs.xlsx'], 'file')  % Check if VDER Excel sheet exists
    MATDSSVDERsSheetnames =  cellstr(sheetnames('MATDSS_VDERs.xlsx'));  % Get sheet names from the file
    MATDSSVDERsFlag = MATDSS_StrComp(MATDSSVDERsSheetnames,app.MATDSSExcelSheet);  % Compare to ensure the correct sheet exists
    if MATDSSVDERsFlag > 0
        MATDSSVDERsExcelSheetMissingFlag = false;  % If sheet exists, set missing flag to false
        VDERsReadOptions = detectImportOptions('MATDSS_VDERs.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % Configure import settings
        VDERsReadOptions = setvartype(VDERsReadOptions,VDERsReadOptions.VariableNames(:),'string');  % Set variable types to string
        clear VDERsTable
        VDERsTable = readtable('MATDSS_VDERs.xlsx',VDERsReadOptions);  % Read the table from the Excel sheet
        VDERsTableContent = VDERsTable.Variables;  % Extract table data
        VDERsTableDefaultValues = ["";"";"1e-6";"Wye";"DSS_load"; '[10,0,0]'; '[2,0,0]'; '1';'1'];  % Default values for VDERs
        if size(VDERsTable,1)  % If the table is not empty
            for i = 1:size(VDERsTable.Variables,2)
                VDERsTableContent(find(ismissing(VDERsTableContent(:,i))),i) = VDERsTableDefaultValues(i);  % Fill missing values
            end
        end
        VDERsTable.Variables = VDERsTableContent;  % Update table data
        app.VDERsTable.Data = VDERsTable;  % Update app table with new data
        app.VDERsTable.ColumnName = VDERsTable.Properties.VariableNames;  % Set column names in the app table
        app.VDERsTable.ColumnWidth = 'auto';  % Auto-adjust column widths
    end
end

if MATDSSVDERsExcelSheetMissingFlag
    % If the VDER sheet is missing, create a new one with default headers and empty data
    AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames;  % Get all bus names for control areas

    VDERsTableHeaders = {'CA#','VDER Name', 'Tau', ['Y/' char(916)], 'DER Mode', 'P(x)', 'Q(x)', 'ax', 'cx'};  % Define headers
    VDERsTableVarTypes = repmat({'string'},1,length(VDERsTableHeaders));  % Set variable types as strings
    VDERsTableSize = [0,length(VDERsTableHeaders)];  % Initialize table with zero rows
    VDERsTable = table('Size', [0, length(VDERsTableHeaders)], 'VariableTypes', VDERsTableVarTypes, 'VariableNames', VDERsTableHeaders);  % Create table
    writetable(VDERsTable,[pwd '\Config\MATDSS_VDERs.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);  % Write table to Excel
    VDERsTable.Properties.VariableNames = VDERsTableHeaders;  % Assign headers to the table
    app.VDERsTable.Data = VDERsTable;  % Update app table with new data
    app.VDERsTable.ColumnName = VDERsTable.Properties.VariableNames;  % Set column names in the app table
    app.VDERsTable.ColumnWidth = 'auto';  % Auto-adjust column widths
    app.VDERsTable.Tag = 'Loaded';  % Mark the table as loaded
end

% Assign the final VDERsTable data to app.VDERsTableData for use in other parts of the app
app.VDERsTableData = VDERsTable;

% Update the app status to indicate successful loading of the configuration for the selected OpenDSS file
MATDSSApp_Status(app,'Ready', ['Configurations for "' app.OpenDSSFilesListBox.Value '" loaded successfully'])

% clearvars -except app  % Clear all variables except 'app' if necessary (commented out)



%% Old function Code

%{

function MATDSSApp_Configurations(app)
% This function loads the configurations saved in Excel files to MATDSS
% Application. Those configurations are not used to generate Plot
% properties panel by default!

MATDSSApp_Status(app,'Loading Configruation Files', ['Loading configuration files for "' app.OpenDSSFilesListBox.Value '"']);


% Tag all tabs as unloaded
for i = 1:size(app.CircuitConfigurationsTabGroup.Children,1)
    app.CircuitConfigurationsTabGroup.Children(i).Tag = 'Not Loaded';
end

% Read DER Information from MATDSS_DER.xlsx (if corresponding
% sheet is not available, the program will make an empty one)
app.MATDSSExcelSheet = app.OpenDSSFilesListBox.Value; % Target sheet in the Excel file




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DER Tab Setup                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSSApp_Details(app, 'Loading DER table');
DERTableHeaders = {'Index', 'Name', 'Bus (# or ''Bus Name'')', 'Tau', 'DER Type', 'Nodes', 'Connection Type', 'Nphase', 'Mode', 'P(x)', 'Q(x)','Pmin','Pmax','Qmin','Qmax', 'ax', 'cx'};
DERTableVarTypes = repmat({'string'},1,length(DERTableHeaders));

if ~exist([pwd '\Config\MATDSS_DER.xlsx'], 'file')
    DERTable = table('Size', [0, length(DERTableHeaders)], 'VariableTypes', DERTableVarTypes, 'VariableNames', DERTableHeaders);
    writetable(DERTable,[pwd '\Config\MATDSS_DER.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
end
MATDSSDERSheetnames =  cellstr(sheetnames('MATDSS_DER.xlsx')); % List of all sheets in the Excel file
MATDSSDERFlag = MATDSS_StrComp(MATDSSDERSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file

if MATDSSDERFlag <= 0 % if Sheet is not available
    DERTable = table('Size', [0, length(DERTableHeaders)], 'VariableTypes', DERTableVarTypes, 'VariableNames', DERTableHeaders);
    writetable(DERTable,[pwd '\Config\MATDSS_DER.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
end

DERReadOptions = detectImportOptions('MATDSS_DER.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
DERReadOptions = setvartype(DERReadOptions,DERReadOptions.VariableNames(:),'string');
clear DERTable
DERTable = readtable('MATDSS_DER.xlsx',DERReadOptions);
DERTableContent = DERTable.Variables;
DERTableDefaultValues = ["";"DER_";"";"0.2";"PV";"[1,2,3]";"Wye";"3";"DSS_load";"[2,0,0]";"[2,0,0]";"-1e6";"1e6";"-1e6";"1e6";"1";"1"];
for i = 1:size(DERTable.Variables,2)
    if i == 1
        DERTableContent(find(ismissing(DERTableContent(:,i))),i) = num2str(find(ismissing(DERTableContent(:,i))));
    elseif i == 2
        DERTableContent(find(ismissing(DERTableContent(:,i))),i) = strcat(DERTableDefaultValues(i), num2str(find(ismissing(DERTableContent(:,i)))));
    else
        DERTableContent(find(ismissing(DERTableContent(:,i))),i) = DERTableDefaultValues(i);
    end
end

DERTable.Variables = DERTableContent;
app.DERsTable.Data = DERTable;
app.DERsTable.ColumnName = DERTable.Properties.VariableNames;
app.DERsTable.ColumnWidth = 'auto';
app.DERsTab.Tag = 'Loaded';




% Run MATDSS quick simulation to generate list of Buses,
% Phases and Branches

[MATDSS] = MATDSS_Initialize(app,1);
MATDSS.Table.DER = DERTable;
MATDSS.TableData.DER =DERTableContent;
if MATDSS.Sim.DSSStartOk % if DSS is loaded and linked successfully
    % Defining a DER device with default values at bus #
    %     DER_Bus = 52; % IEEE123Master
    [MATDSS, DER] = MATDSS_DER(MATDSS);
    AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
    AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames; AllBusNames(1) = [];

    % Get all lines/branches names
    MyLines = MATDSS.Sim.DSSCircuit.Lines;
    MyLinesAllNames = MyLines.AllNames;
    MyLinesAllPhasesNames = {};
    for i = 1:size(MyLinesAllNames,1)
        MyLines.Name = MyLinesAllNames{i};
        nphases = MyLines.Phases;
        for j = 1:nphases
            MyLinesAllPhasesNames = [MyLinesAllPhasesNames;[MyLinesAllNames{i} '.' num2str(j)]];
        end
    end
else
    AllNodesNames = {'Error! Could not initiate the link with OpenDSS'};
    MyLinesAllPhasesNames = {'Error! Could not initiate the link with OpenDSS'};
end

% This code was depricated with the introduction of multi-control area
% structure. it is kept here for reference only. In future releases, the
% code will be refined and this will be removed.


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                 L2C Tab Setup                 %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% MATDSSL2CExcelSheetMissingFlag = true;
% if exist([pwd '\Config\MATDSS_L2C.xlsx'], 'file')
%     MATDSSL2CSheetnames =  cellstr(sheetnames('MATDSS_L2C.xlsx')); % List of all sheets in the Excel file
%     MATDSSL2CFlag = MATDSS_StrComp(MATDSSL2CSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file
%     if MATDSSL2CFlag > 0
%         MATDSSL2CExcelSheetMissingFlag = false;
%         L2CReadOptions = detectImportOptions('MATDSS_L2C.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
%         L2CReadOptions = setvartype(L2CReadOptions,L2CReadOptions.VariableNames(:),'string');
%         L2CTable = readtable([pwd '\Config\MATDSS_L2C.xlsx'],L2CReadOptions);%'Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve',);
%         L2CTableData = L2CTable.Variables;
%         app.alphaEditField.Value = L2CTableData(1);
%         app.rpEditField.Value = L2CTableData(2);
%         app.rbardEditField.Value = L2CTableData(3);
%         app.EEditField.Value = L2CTableData(4);
%         app.vulEditField.Value = L2CTableData(5);
%         app.vllEditField.Value = L2CTableData(6);
%         app.iulEditField.Value = L2CTableData(7);
%         app.a_rhoEditField.Value = L2CTableData(8);
%         app.a_lambdaEditField.Value = L2CTableData(9);
%         app.a_muEditField.Value = L2CTableData(10);
%         app.a_gammaEditField.Value = L2CTableData(11);
%         app.a_nuEditField.Value = L2CTableData(12);
%         app.a_zetaEditField.Value = L2CTableData(13);
%         app.c_rhoEditField.Value = L2CTableData(14);
%         app.c_lambdaEditField.Value = L2CTableData(15);
%         app.c_muEditField.Value = L2CTableData(16);
%         app.c_gammaEditField.Value = L2CTableData(17);
%         app.c_nuEditField.Value = L2CTableData(18);
%         app.c_zetaEditField.Value = L2CTableData(19);
%     end
% end
% if MATDSSL2CExcelSheetMissingFlag
%     L2CTableHeaders = {'alpha','rp','rbard','E (W)','vul (p.u.)','vll (p.u.)','iul (p.u.)','a_rho','a_lambda','a_mu','a_gamma (W^2/V^2)','a_nu (W^2/V^2)','a_zeta (W^2/A^2)','c_rho','c_lambda','c_mu','c_gamma (V^2/W^2)','c_nu (V^2/W^2)', 'c_zeta (A^2/W^2)'};
%     L2CTableVarTypes = repmat({'string'},1,length(L2CTableHeaders));
%     L2CTableSize = [1,length(L2CTableHeaders)];
%     L2CTableData = repmat("",L2CTableSize);
%     L2CTableData(1,1) = app.alphaEditField.Value;
%     L2CTableData(1,2) = app.rpEditField.Value;
%     L2CTableData(1,3) = app.rbardEditField.Value;
%     L2CTableData(1,4) = app.EEditField.Value;
%     L2CTableData(1,5) = app.vulEditField.Value;
%     L2CTableData(1,6) = app.vllEditField.Value;
%     L2CTableData(1,7) = app.iulEditField.Value;
%     L2CTableData(1,8) = app.a_rhoEditField.Value;
%     L2CTableData(1,9) = app.a_lambdaEditField.Value;
%     L2CTableData(1,10) = app.a_muEditField.Value;
%     L2CTableData(1,11) = app.a_gammaEditField.Value;
%     L2CTableData(1,12) = app.a_nuEditField.Value;
%     L2CTableData(1,13) = app.a_zetaEditField.Value;
%     L2CTableData(1,14) = app.c_rhoEditField.Value;
%     L2CTableData(1,15) = app.c_lambdaEditField.Value;
%     L2CTableData(1,16) = app.c_muEditField.Value;
%     L2CTableData(1,17) = app.c_gammaEditField.Value;
%     L2CTableData(1,18) = app.c_nuEditField.Value;
%     L2CTableData(1,19) = app.c_zetaEditField.Value;
%     L2CTable = cell2table(cellstr(L2CTableData));
%     L2CTable.Properties.VariableNames = L2CTableHeaders;
%     writetable(L2CTable,[pwd '\Config\MATDSS_L2C.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
% end
% app.L2CTab.Tag = 'Loaded';
% 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               V & I Tab Setup                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSSApp_Details(app, 'Loading V & I tables');
if size(AllNodesNames,1) == 1 % Error in OpenDSS Link
    app.VTrackingListBox.Items = AllNodesNames;
    app.ITrackingListBox.Items = MyLinesAllPhasesNames;
else
    app.VTrackingListBox.Items = AllNodesNames(3+1:end);
    app.ITrackingListBox.Items = MyLinesAllPhasesNames;

    if exist([pwd '\Config\MATDSS_VI.xlsx'], 'file')
        % Read MATDSS_VI.xlsx file if available to
        % select the saved buses and lines.
        MATDSSVISheetnames =  cellstr(sheetnames('MATDSS_VI.xlsx')); % List of all sheets in the Excel file
        MATDSSVIFlag = MATDSS_StrComp(MATDSSVISheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file

        if MATDSSVIFlag > 0
            VIReadOptions = detectImportOptions('MATDSS_VI.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
            VIReadOptions = setvartype(VIReadOptions,VIReadOptions.VariableNames(:),'string');
            VITable = readtable([pwd '\Config\MATDSS_VI.xlsx'],VIReadOptions);
            VITableData = VITable.Variables;

            Mv = VITableData(:,1);
            Mv(ismissing(Mv)) = [];
            v_size = length(Mv);
            A = app.VTrackingListBox.Items';
            B = convertStringsToChars(VITableData(1:v_size,1));
            [C, iA, iB] = intersect(A,B, 'stable');
            
            k = iA;
            if k >= 1
                app.VTrackingListBox.Value = app.VTrackingListBox.Items(k);
            end

            Mi = VITableData(:,3);
            Mi(ismissing(Mi)) = [];
            i_size = length(Mi);

            A = app.ITrackingListBox.Items';
            B = convertStringsToChars(VITableData(1:i_size,3));
            [C, iA, iB] = intersect(A,B, 'stable');
            k = iA;
            if k >= 1
                app.ITrackingListBox.Value = app.ITrackingListBox.Items(k);
            end
        end
    end
end
app.VITab.Tag = 'Loaded';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Control Areas Tab Setup             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MATDSSApp_Details(app, 'Loading controllers tables');
MATDSSControlAreasExcelSheetMissingFlag = true;
if exist([pwd '\Config\MATDSS_ControlAreas.xlsx'], 'file')
    MATDSSControlAreasSheetnames =  cellstr(sheetnames('MATDSS_ControlAreas.xlsx')); % List of all sheets in the Excel file
    MATDSSControlAreasFlag = MATDSS_StrComp(MATDSSControlAreasSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file
    if MATDSSControlAreasFlag > 0
        MATDSSControlAreasExcelSheetMissingFlag = false;
        ControlAreasReadOptions = detectImportOptions('MATDSS_ControlAreas.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
        ControlAreasReadOptions = setvartype(ControlAreasReadOptions,ControlAreasReadOptions.VariableNames(:),'string');
        clear ControlAreasTable
        ControlAreasTable = readtable('MATDSS_ControlAreas.xlsx',ControlAreasReadOptions);
        ControlAreasTableContent = ControlAreasTable.Variables;
        ControlAreasTableDefaultValues = ["";"";"1";"F";""];
        for i = 1:size(ControlAreasTable.Variables,2)
            ControlAreasTableContent(find(ismissing(ControlAreasTableContent(:,i))),i) = ControlAreasTableDefaultValues(i);
        end

        ControlAreasTable.Variables = ControlAreasTableContent;
        app.ControlAreasTable.Data = ControlAreasTable;
        app.ControlAreasTable.ColumnName = ControlAreasTable.Properties.VariableNames;
        app.ControlAreasTable.ColumnWidth = 'auto';
    end
end

if MATDSSControlAreasExcelSheetMissingFlag
    AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames;

    % Control Areas Settings Section
    ControlAreasTableHeaders = {'Bus Name','Bus Number', 'Assigned Control Area (#)', 'Interface Bus (T/F)','Connected to Bus name'};
    ControlAreasTableVarTypes = repmat({'string'},1,length(ControlAreasTableHeaders));
    ControlAreasTableSize = [length(AllBusNames),length(ControlAreasTableHeaders)];
    ControlAreasTableData = repmat("",ControlAreasTableSize);
    ControlAreasTableData(:,1) = AllBusNames;
    ControlAreasTableData(:,2) = num2str([1:length(AllBusNames)]');
    ControlAreasTableData(:,3) = num2str(1);
    ControlAreasTableData(:,4) = 'F';
    ControlAreasTableData(:,5) = '';
    ControlAreasTable = cell2table(cellstr(ControlAreasTableData));
    ControlAreasTable.Properties.VariableNames = ControlAreasTableHeaders;
    writetable(ControlAreasTable,[pwd '\Config\MATDSS_ControlAreas.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
    app.ControlAreasTable.Data = ControlAreasTable;
    app.ControlAreasTable.ColumnName = ControlAreasTable.Properties.VariableNames;
    app.ControlAreasTable.ColumnWidth = 'auto';
    app.ControlAreasTable.Tag = 'Loaded';
end

app.ControlAreasTableData = ControlAreasTable;


%             if MATDSSControlAreaFlag <= 0 % if Sheet is not available
%                 ControlAreasTable = table('Size', [0, length(ControlAreasTableHeaders)], 'VariableTypes', ControlAreasTableVarTypes, 'VariableNames', ControlAreasTableHeaders);
%                 writetable(ControlAreasTable,[pwd '\Config\MATDSS_ControlAreas.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
%             end
%
%             ControlAreasReadOptions = detectImportOptions('MATDSS_ControlAreas.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
%             ControlAreasReadOptions = setvartype(ControlAreasReadOptions,ControlAreasReadOptions.VariableNames(:),'string');
%             clear ControlAreasTable
%             ControlAreasTable = readtable('MATDSS_ControlAreas.xlsx',ControlAreasReadOptions);
%             ControlAreasTableContent = ControlAreasTable.Variables;
%             ControlAreasTableDefaultValues = ["";"1"];
%             for i = 1:size(ControlAreasTable.Variables,2)
%                 ControlAreasTableContent(find(ismissing(ControlAreasTableContent(:,i))),i) = ControlAreasTableDefaultValues(i);
%             end
%
%             ControlAreasTable.Variables = ControlAreasTableContent;
%             app.ControlAreasTable.Data = ControlAreasTable;
%             app.ControlAreasTable.ColumnName = ControlAreasTable.Properties.VariableNames;
%             app.ControlAreasTable.ColumnWidth = 'auto';
%             app.ControlAreasTable.Tag = 'Loaded';

ControllersSettingsTableDefaultValues = ["1";"LLC";"auto";"1e-4";"1e-3";"1e2";"1.05";"0.95";"Specified in DSS File";...
    "1";"1"; "1e3"; "1e3"; "1e3"; "1e3"; "1e12"; "1e12"; "1e8";...
    "1";"1";"1e-3";"1e-3";"1e-3";"1e-3";"1e-12";"1e-12";"1e-8"];

% Controllers Settings Section
MATDSSControllersSettingsExcelSheetMissingFlag = true;
if exist([pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'file')
    MATDSSControllersSettingsSheetnames =  cellstr(sheetnames('MATDSS_ControllersSettings.xlsx')); % List of all sheets in the Excel file
    MATDSSControllersSettingsFlag = MATDSS_StrComp(MATDSSControllersSettingsSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file
    if MATDSSControllersSettingsFlag > 0
        MATDSSControllersSettingsExcelSheetMissingFlag = false;
        ControllersSettingsReadOptions = detectImportOptions('MATDSS_ControllersSettings.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
        ControllersSettingsReadOptions = setvartype(ControllersSettingsReadOptions,ControllersSettingsReadOptions.VariableNames(:),'string');
        clear ControllersSettingsTable
        ControllersSettingsTable = readtable('MATDSS_ControllersSettings.xlsx',ControllersSettingsReadOptions);
        ControllersSettingsTableContent = ControllersSettingsTable.Variables;
        for i = 1:size(ControllersSettingsTable.Variables,2)
            ControllersSettingsTableContent(find(ismissing(ControllersSettingsTableContent(:,i))),i) = ControllersSettingsTableDefaultValues(i);
        end

        ControllersSettingsTable.Variables = ControllersSettingsTableContent;
        app.ControllersSettingsTable.Data = ControllersSettingsTable;
        app.ControllersSettingsTable.ColumnName = ControllersSettingsTable.Properties.VariableNames;
        app.ControllersSettingsTable.ColumnWidth = 'auto';
    end
end

if MATDSSControllersSettingsExcelSheetMissingFlag
    % Control Areas Settings Section
    ControllersSettingsTableHeaders = {'# Control Area','Controller Type','alpha','r_p','rbar_d','E (W)','v_ul (p.u.)','v_ll (p.u.)','i_ul (p.u.)',...
        'a_rho','a_sigma','a_lambda','a_mu','a_eta','a_psi','a_gamma (W^2/V^2)','a_nu (W^2/V^2)','a_zeta (W^2/A^2)',...
        'c_rho','c_sigma','c_lambda','c_mu','c_eta','c_psi','c_gamma (V^2/W^2)','c_nu (V^2/W^2)','c_zeta (A^2/W^2)'};
    ControllersSettingsTableSize = [1,length(ControllersSettingsTableHeaders)];
    ControllersSettingsTableData = repmat("",ControllersSettingsTableSize);
    ControllersSettingsTableData(1,:) = ControllersSettingsTableDefaultValues;
    ControllersSettingsTable = cell2table(cellstr(ControllersSettingsTableData));
    ControllersSettingsTable.Properties.VariableNames = ControllersSettingsTableHeaders;
    writetable(ControllersSettingsTable,[pwd '\Config\MATDSS_ControllersSettings.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
    app.ControllersSettingsTable.Data = ControllersSettingsTable;
    app.ControllersSettingsTable.ColumnName = ControllersSettingsTable.Properties.VariableNames;
    app.ControllersSettingsTable.ColumnWidth = 'auto';
    app.ControllersSettingsTable.Tag = 'Loaded';
end

% VDER Configs
MATDSSVDERsExcelSheetMissingFlag = true;
if exist([pwd '\Config\MATDSS_VDERs.xlsx'], 'file')
    MATDSSVDERsSheetnames =  cellstr(sheetnames('MATDSS_VDERs.xlsx')); % List of all sheets in the Excel file
    MATDSSVDERsFlag = MATDSS_StrComp(MATDSSVDERsSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file
    if MATDSSVDERsFlag > 0
        MATDSSVDERsExcelSheetMissingFlag = false;
        VDERsReadOptions = detectImportOptions('MATDSS_VDERs.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
        VDERsReadOptions = setvartype(VDERsReadOptions,VDERsReadOptions.VariableNames(:),'string');
        clear VDERsTable
        VDERsTable = readtable('MATDSS_VDERs.xlsx',VDERsReadOptions);
        VDERsTableContent = VDERsTable.Variables;
        VDERsTableDefaultValues = ["";"";"1e-6";"Wye";"DSS_load"; '[10,0,0]'; '[2,0,0]'; '1';'1'];
        if size(VDERsTable,1)
            for i = 1:size(VDERsTable.Variables,2)
                VDERsTableContent(find(ismissing(VDERsTableContent(:,i))),i) = VDERsTableDefaultValues(i);
            end
        end
        VDERsTable.Variables = VDERsTableContent;
        app.VDERsTable.Data = VDERsTable;
        app.VDERsTable.ColumnName = VDERsTable.Properties.VariableNames;
        app.VDERsTable.ColumnWidth = 'auto';
    end
end

if MATDSSVDERsExcelSheetMissingFlag
    AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames;

    % Control Areas Settings Section
    VDERsTableHeaders = {'CA#','VDER Name', 'Tau', ['Y/' char(916)], 'DER Mode', 'P(x)', 'Q(x)', 'ax', 'cx'};
    VDERsTableVarTypes = repmat({'string'},1,length(VDERsTableHeaders));
    VDERsTableSize = [0,length(VDERsTableHeaders)];
    VDERsTable = table('Size', [0, length(VDERsTableHeaders)], 'VariableTypes', VDERsTableVarTypes, 'VariableNames', VDERsTableHeaders);
    writetable(VDERsTable,[pwd '\Config\MATDSS_VDERs.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
    VDERsTable.Properties.VariableNames = VDERsTableHeaders;
    app.VDERsTable.Data = VDERsTable;
    app.VDERsTable.ColumnName = VDERsTable.Properties.VariableNames;
    app.VDERsTable.ColumnWidth = 'auto';
    app.VDERsTable.Tag = 'Loaded';
end



app.VDERsTableData = VDERsTable;

MATDSSApp_Status(app,'Ready', ['Configurations for "' app.OpenDSSFilesListBox.Value '" loaded successfully'])


% clearvars -except app


%}