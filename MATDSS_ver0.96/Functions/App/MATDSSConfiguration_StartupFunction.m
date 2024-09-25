function MATDSSConfiguration_StartupFunction(app)
% Tag all tabs as unloaded
for i = 1:size(app.TabGroup.Children,1)
    app.TabGroup.Children(i).Tag = 'Not Loaded';
end
% Fixing the window location and size
myscreen = get(0,'ScreenSize'); % Get screen resolution info
mywindow = [1500,800]; % Set app window size here
sp = [(myscreen(3)-mywindow(1))/2, (myscreen(4)-mywindow(2))/2]; % Calculate the correct starting point to get the app centered on the screen
app.MATDSSTuneAssistantUIFigure.Position = [sp,mywindow]; % set the correct info in the position varible of the app!

% Read DER Information from MATDSS_DER.xlsx (if corresponding
% sheet is not available, the program will make an empty one)
app.MATDSSExcelSheet = app.CallingApp.OpenDSSFilesListBox.Value; % Target sheet in the Excel file




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DER Tab Setup                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DERTableHeaders = {'Index', 'Name', 'Bus (# or ''Bus Name'')', 'Tau', 'DER Type', 'Nodes', 'Connection Type', 'Nphase', 'Mode', 'P(x)', 'Q(x)','Pmin','Pmax','Qmin','Qmax', 'ax', 'cx'};
DERTableVarTypes = repmat({'string'},1,length(DERTableHeaders));

if ~exist([pwd '\Config\MATDSS_DER.xlsx'], 'file')
    DERTable = table('Size', [0, length(DERTableHeaders)], 'VariableTypes', DERTableVarTypes, 'VariableNames', DERTableHeaders);
    writetable(DERTable,[pwd '\Config\MATDSS_DER.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
end
MATDSSDERSheetnames =  cellstr(sheetnames('MATDSS_DER.xlsx')); % List of all sheets in the Excel file
MATDSSDERFlag = MATDSS_StrComp(MATDSSDERSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file

if MATDSSDERFlag <= 0 % if Sheet is not availableha
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

[MATDSS, DER] = MATDSS_Initialize(app.CallingApp,1);
if MATDSS.Sim.DSSStartOk % if DSS is loaded and linked successfully
    % Defining a DER device with default values at bus #
    %     DER_Bus = 52; % IEEE123Master
    [MATDSS, DER] = MATDSS_DER(app.CallingApp,MATDSS,DER);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 L2C Tab Setup                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MATDSSL2CExcelSheetMissingFlag = true;
if exist([pwd '\Config\MATDSS_L2C.xlsx'], 'file')
    MATDSSL2CSheetnames =  cellstr(sheetnames('MATDSS_L2C.xlsx')); % List of all sheets in the Excel file
    MATDSSL2CFlag = MATDSS_StrComp(MATDSSL2CSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file
    if MATDSSL2CFlag > 0
        MATDSSL2CExcelSheetMissingFlag = false;
        L2CReadOptions = detectImportOptions('MATDSS_L2C.xlsx','Sheet',app.MATDSSExcelSheet,'VariableNamingRule','preserve'); % set the import settings
        L2CReadOptions = setvartype(L2CReadOptions,L2CReadOptions.VariableNames(:),'string');
        L2CTable = readtable([pwd '\Config\MATDSS_L2C.xlsx'],L2CReadOptions);%'Sheet',app.OpenDSSFilesListBox.Value,'VariableNamingRule','preserve',);
        L2CTableData = L2CTable.Variables;
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
if MATDSSL2CExcelSheetMissingFlag
    L2CTableHeaders = {'alpha','rp','rbard','E (W)','vul (p.u.)','vll (p.u.)','iul (p.u.)','a_rho','a_lambda','a_mu','a_gamma (W^2/V^2)','a_nu (W^2/V^2)','a_zeta (W^2/A^2)','c_rho','c_lambda','c_mu','c_gamma (V^2/W^2)','c_nu (V^2/W^2)', 'c_zeta (A^2/W^2)'};
    L2CTableVarTypes = repmat({'string'},1,length(L2CTableHeaders));
    L2CTableSize = [1,length(L2CTableHeaders)];
    L2CTableData = repmat("",L2CTableSize);
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
    L2CTable = cell2table(cellstr(L2CTableData));
    L2CTable.Properties.VariableNames = L2CTableHeaders;
    writetable(L2CTable,[pwd '\Config\MATDSS_L2C.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
end
app.L2CTab.Tag = 'Loaded';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               V & I Tab Setup                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            k = MATDSS_StrComp(app.VTrackingListBox.Items,convertStringsToChars(VITableData(1:v_size,1)));
            if k >= 1
                app.VTrackingListBox.Value = app.VTrackingListBox.Items(k);
            end

            Mi = VITableData(:,3);
            Mi(ismissing(Mi)) = [];
            i_size = length(Mi);
            k = MATDSS_StrComp(app.ITrackingListBox.Items,convertStringsToChars(VITableData(1:i_size,3)));
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




end