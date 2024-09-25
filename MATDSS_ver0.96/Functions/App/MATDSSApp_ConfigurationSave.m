function MATDSSApp_ConfigurationSave(app, option)
% MATDSSAppConfiguration_Save(app, option)
% This function saves the specified configurations in MATDSS
% Configuration window of the selected circuit to the corresponding Excel
% files based on the selected tab in the configuration window.
%
% Parameters:
%   - app: The application object containing data to be saved.
%   - option: Specifies the tab or section to be saved (e.g., 'DERs', 'V & I', 'Control Areas').
%
% Last Update: MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

switch option
    case 'DERs'
        % Save DERs table data to an Excel file
        % Define file path and save the table data to the specified sheet
        writetable(app.DERsTable.Data, [pwd '\Config\MATDSS_DER.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        % Notify the user that the DERs table has been saved
        MATDSSApp_Details(app, 'DERs Table Saved');
        
    case 'V & I'
        % Headers for the Voltage and Current Table
        VITableHeaders = {'Mv', 'Mv_index', 'Mi', 'Mi_index', 'VNodes', 'ILines'};
        
        % Retrieve selected nodes and lines from the tracking list boxes
        MvList = app.VTrackingListBox.Value;
        MiList = app.ITrackingListBox.Value;
        
        % Initialize MATDSS structure and check if DSS is loaded successfully
        [MATDSS] = MATDSS_Initialize(app, 1);
        
        if MATDSS.Sim.DSSStartOk
            % Update MATDSS table with DER datadddd
            MATDSS.Table.DER = app.DERsTable.Data;
            MATDSS.TableData.DER = app.DERsTable.Data.Variables;
            [MATDSS, DER] = MATDSS_DER(MATDSS);
            
            % Retrieve all node names and line names from DSS circuit
            AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
            MyLines = MATDSS.Sim.DSSCircuit.Lines;
            AllLinesNames = MyLines.AllNames;
            AllLinesAllPhasesNames = {};
            for i = 1:size(AllLinesNames, 1)
                MyLines.Name = AllLinesNames{i};
                nphases = MyLines.Phases;
                for j = 1:nphases
                    AllLinesAllPhasesNames = [AllLinesAllPhasesNames; [AllLinesNames{i} '.' num2str(j)]];
                end
            end
        end

        % Find indices of selected nodes and lines
        k_v = MATDSS_StrComp(AllNodesNames, MvList);
        k_i = MATDSS_StrComp(AllLinesAllPhasesNames, MiList);

        % Create table data with selected and all nodes and lines
        VITableSize = [max(length(AllNodesNames), length(AllLinesAllPhasesNames)), length(VITableHeaders)];
        VITableData = repmat("", VITableSize);
        VITableData(1:length(MvList), 1) = app.VTrackingListBox.Value';
        VITableData(1:length(MvList), 2) = k_v;
        VITableData(1:length(MiList), 3) = app.ITrackingListBox.Value';
        VITableData(1:length(MiList), 4) = k_i;
        VITableData(1:length(AllNodesNames), 5) = AllNodesNames;
        VITableData(1:length(AllLinesAllPhasesNames), 6) = AllLinesAllPhasesNames;
        VITable = cell2table(cellstr(VITableData));
        VITable.Properties.VariableNames = VITableHeaders;
        
        % Save the Voltage and Current Table data to an Excel file
        writetable(VITable, [pwd '\Config\MATDSS_VI.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        % Notify the user that the VI table has been saved
        MATDSSApp_Details(app, 'VI Table Saved');
        
    case 'Control Areas'
        % Save Control Areas, Controllers Settings, and VDERs tables to their respective Excel files
        writetable(app.ControlAreasTable.Data, [pwd '\Config\MATDSS_ControlAreas.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        writetable(app.ControllersSettingsTable.Data, [pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        writetable(app.VDERsTable.Data, [pwd '\Config\MATDSS_VDERs.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        % Notify the user that the control areas, controllers settings, and VDERs tables have been saved
        MATDSSApp_Details(app, 'Controllers, Areas and VDERs Tables saved');
        
    otherwise
        % If the option is not recognized, do nothing
        % Option can be expanded or handled in future versions
end



%% Old function code
%{

function MATDSSAppConfiguration_Save(app, option)
% MATDSSAppConfiguration_Save(app, option)
% This function will save the specified configurations in MATDSS
% Configuration window of the selected circuit in the corresponding Excel
% files.
%
% option currently is set for the "tab" selected to be saved from the
% configuration window.
%
%   Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
%   MATDSS Application
%   Copyright (c) 2024, Ilyas Farhat
%   by Ilyas Farhat
%
%   This file is part of MATDSS Application
%   Contact the developer at ilyas.farhat@outlook.com

switch option
    case 'DERs'
        writetable(app.DERsTable.Data, [pwd '\Config\MATDSS_DER.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        MATDSSApp_Details(app, 'DERs Table Saved');
    case 'V & I'
        VITableHeaders = {'Mv', 'Mv_index', 'Mi', 'Mi_index', 'VNodes', 'ILines'};
        MvList = app.VTrackingListBox.Value;
        MiList = app.ITrackingListBox.Value;
        
        [MATDSS] = MATDSS_Initialize(app, 1);
        
        if MATDSS.Sim.DSSStartOk % if DSS is loaded and linked successfully
            MATDSS.Table.DER = app.DERsTable.Data;
            MATDSS.TableData.DER = app.DERsTable.Data.Variables;
            [MATDSS, DER] = MATDSS_DER(MATDSS);

            AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
            MyLines = MATDSS.Sim.DSSCircuit.Lines;
            AllLinesNames = MyLines.AllNames;
            AllLinesAllPhasesNames = {};
            for i = 1:size(AllLinesNames, 1)
                MyLines.Name = AllLinesNames{i};
                nphases = MyLines.Phases;
                for j = 1:nphases
                    AllLinesAllPhasesNames = [AllLinesAllPhasesNames; [AllLinesNames{i} '.' num2str(j)]];
                end
            end
        end

        k_v = MATDSS_StrComp(AllNodesNames, MvList);
        k_i = MATDSS_StrComp(AllLinesAllPhasesNames, MiList);

        VITableSize = [max(length(AllNodesNames), length(AllLinesAllPhasesNames)), length(VITableHeaders)];
        VITableData = repmat("", VITableSize);
        VITableData(1:length(MvList), 1) = app.VTrackingListBox.Value';
        VITableData(1:length(MvList), 2) = k_v;
        VITableData(1:length(MiList), 3) = app.ITrackingListBox.Value';
        VITableData(1:length(MiList), 4) = k_i;
        VITableData(1:length(AllNodesNames), 5) = AllNodesNames;
        VITableData(1:length(AllLinesAllPhasesNames), 6) = AllLinesAllPhasesNames;
        VITable = cell2table(cellstr(VITableData));
        VITable.Properties.VariableNames = VITableHeaders;
        writetable(VITable, [pwd '\Config\MATDSS_VI.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        MATDSSApp_Details(app, 'VI Table Saved');
    case 'Control Areas'
        writetable(app.ControlAreasTable.Data, [pwd '\Config\MATDSS_ControlAreas.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        writetable(app.ControllersSettingsTable.Data, [pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        writetable(app.VDERsTable.Data, [pwd '\Config\MATDSS_VDERs.xlsx'], 'Sheet', app.MATDSSExcelSheet, 'WriteMode', 'overwritesheet', 'autofitwidth', 1);
        MATDSSApp_Details(app, 'Controllers, Areas and VDERs Tables saved');
    otherwise
        % Do nothing, option is not recognized.
end

%}