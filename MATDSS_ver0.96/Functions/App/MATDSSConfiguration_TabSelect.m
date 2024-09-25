function MATDSSConfiguration_TabSelect(app)
selectedTab = app.TabGroup.SelectedTab;
switch selectedTab.Title
    case 'DER'

    case 'L2C'
        if app.L2CTab.Tag ~= "Loaded"
            if exist([pwd '\Config\MATDSS_L2C.xlsx'], 'file')
                MATDSSL2CSheetnames =  cellstr(sheetnames('MATDSS_L2C.xlsx')); % List of all sheets in the Excel file
                MATDSSL2CFlag = MATDSS_StrComp(MATDSSL2CSheetnames,app.MATDSSExcelSheet); % Check if teh target sheet is available in the Excel file
                if MATDSSL2CFlag > 0
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
            app.L2CTab.Tag = 'Loaded';
        end
    case 'V & I'
        if app.VITab.Tag ~= "Loaded"
            if isempty(app.VTrackingListBox.Items)
                [MATDSS, DER] = MATDSS_Initialize(app.CallingApp,1);
                if MATDSS.Sim.DSSStartOk % if DSS is loaded and linked successfully
                    % Defining a DER device with default values at bus #
                    %     DER_Bus = 52; % IEEE123Master
                    [MATDSS, DER] = MATDSS_DER(app.CallingApp,MATDSS,DER);
                    AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
                    app.VTrackingListBox.Items = AllNodesNames(3+1:end);

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
                            v_size = 1;
                            for i = 1:size(VITableData,1)
                                if VITableData(i,1) == "" || ismissing(VITableData(i,1))
                                    break;
                                end
                                v_size = i;
                            end
                            k = MATDSS_StrComp(app.VTrackingListBox.Items,convertStringsToChars(VITableData(1:v_size,1)));
                            if k >= 1
                                app.VTrackingListBox.Value = app.VTrackingListBox.Items(k);
                            end



                            i_size = 1;
                            for i = 1:size(VITableData,1)
                                if VITableData(i,3) == ""  || ismissing(VITableData(i,3))
                                    break;
                                end
                                i_size = i;
                            end
                            k = MATDSS_StrComp(app.ITrackingListBox.Items,convertStringsToChars( VITableData(1:i_size,3)));
                            if k >= 1
                                app.ITrackingListBox.Value = app.ITrackingListBox.Items(k);
                            end
                        end
                    end
                else
                    app.VTrackingListBox.Items = 'Error! Could not initiate the link with OpenDSS';
                    app.ITrackingListBox.Items = 'Error! Could not initiate the link with OpenDSS';
                end
            end
            app.VITab.Tag = 'Loaded';
        end
end
end
