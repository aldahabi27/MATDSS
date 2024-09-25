function MATDSSConfiguration_Save(app,option)
% MATDSSConfiguration_Save(app,option)
% This function will save the specified configurations in MATDSS
% Configuration window of the selected circuit in the corresponding Excel
% files.
%
% option currently is set for the "tab" selected to be saved from the
% configuration window.


switch option
    case 'DERs'
        writetable(app.DERsTable.Data,[pwd '\Config\MATDSS_DER.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
%     case 'L2C'
%         L2CTableHeaders = {'alpha','rp','rbard','E (W)','vul (p.u.)','vll (p.u.)','iul (p.u.)','a_rho','a_lambda','a_mu','a_gamma (W^2/V^2)','a_nu (W^2/V^2)','a_zeta (W^2/A^2)','c_rho','c_lambda','c_mu','c_gamma (V^2/W^2)','c_nu (V^2/W^2)', 'c_zeta (A^2/W^2)'};
%         L2CTableSize = [1,length(L2CTableHeaders)];
%         L2CTableData = repmat("",L2CTableSize);
%         L2CTableData(1,1) = app.alphaEditField.Value;
%         L2CTableData(1,2) = app.rpEditField.Value;
%         L2CTableData(1,3) = app.rbardEditField.Value;
%         L2CTableData(1,4) = app.EEditField.Value;
%         L2CTableData(1,5) = app.vulEditField.Value;
%         L2CTableData(1,6) = app.vllEditField.Value;
%         L2CTableData(1,7) = app.iulEditField.Value;
%         L2CTableData(1,8) = app.a_rhoEditField.Value;
%         L2CTableData(1,9) = app.a_lambdaEditField.Value;
%         L2CTableData(1,10) = app.a_muEditField.Value;
%         L2CTableData(1,11) = app.a_gammaEditField.Value;
%         L2CTableData(1,12) = app.a_nuEditField.Value;
%         L2CTableData(1,13) = app.a_zetaEditField.Value;
%         L2CTableData(1,14) = app.c_rhoEditField.Value;
%         L2CTableData(1,15) = app.c_lambdaEditField.Value;
%         L2CTableData(1,16) = app.c_muEditField.Value;
%         L2CTableData(1,17) = app.c_gammaEditField.Value;
%         L2CTableData(1,18) = app.c_nuEditField.Value;
%         L2CTableData(1,19) = app.c_zetaEditField.Value;
% 
%         L2CTable = cell2table(cellstr(L2CTableData));
%         L2CTable.Properties.VariableNames = L2CTableHeaders;
%         writetable(L2CTable,[pwd '\Config\MATDSS_L2C.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
    case 'V & I'
        VITableHeaders = {'Mv','Mv_index','Mi','Mi_index','VNodes','ILines'};
        MvList = app.VTrackingListBox.Value;
        MiList = app.ITrackingListBox.Value;

        [MATDSS, DER] = MATDSS_Initialize(app.CallingApp,1);
        if MATDSS.Sim.DSSStartOk % if DSS is loaded and linked successfully
            % Defining a DER device with default values at bus #
            %     DER_Bus = 52; % IEEE123Master
            [MATDSS, DER] = MATDSS_DER(app.CallingApp,MATDSS,DER);

            AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
            % Get all lines/branches names
            MyLines = MATDSS.Sim.DSSCircuit.Lines;
            AllLinesNames = MyLines.AllNames;
            AllLinesAllPhasesNames = {};
            for i = 1:size(AllLinesNames,1)
                MyLines.Name = AllLinesNames{i};
                nphases = MyLines.Phases;
                for j = 1:nphases
                    AllLinesAllPhasesNames = [AllLinesAllPhasesNames;[AllLinesNames{i} '.' num2str(j)]];
                end
            end
        end

        k_v = MATDSS_StrComp(AllNodesNames,MvList);
        k_i = MATDSS_StrComp(AllLinesAllPhasesNames,MiList);


        VITableSize = [max(length(AllNodesNames),length(AllLinesAllPhasesNames)),length(VITableHeaders)];
        VITableData = repmat("",VITableSize);
        VITableData(1:length(MvList),1) = app.VTrackingListBox.Value';
        VITableData(1:length(MvList),2) = k_v;
        VITableData(1:length(MiList),3) = app.ITrackingListBox.Value';
        VITableData(1:length(MiList),4) = k_i;
        VITableData(1:length(AllNodesNames),5) = AllNodesNames;
        VITableData(1:length(AllLinesAllPhasesNames),6) = AllLinesAllPhasesNames;
        VITable = cell2table(cellstr(VITableData));
        VITable.Properties.VariableNames = VITableHeaders;
        writetable(VITable,[pwd '\Config\MATDSS_VI.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
    case 'Control Areas'
        writetable(app.ControlAreasTable.Data,[pwd '\Config\MATDSS_ControlAreas.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
        writetable(app.ControllersSettingsTable.Data,[pwd '\Config\MATDSS_ControllersSettings.xlsx'],'Sheet',app.MATDSSExcelSheet,'WriteMode', 'overwritesheet','autofitwidth', 1);
    otherwise
        %do nothing, it is not possible!
end