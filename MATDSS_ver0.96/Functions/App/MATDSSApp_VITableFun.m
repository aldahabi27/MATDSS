function VITable = MATDSSApp_VITableFun(app,MATDSS)

VITableHeaders = {'Mv','Mv_index','Mi','Mi_index','VNodes','ILines'};
MvList = app.VTrackingListBox.Value;
MiList = app.ITrackingListBox.Value;

MATDSStemp = MATDSS_Initialize(app,1);
if MATDSStemp.Sim.DSSStartOk % if DSS is loaded and linked successfully
    % Defining a DER device with default values at bus #
    %     DER_Bus = 52; % IEEE123Master
    MATDSStemp.Table = MATDSS.Table;
    MATDSStemp.TableData = MATDSS.TableData;
    [MATDSStemp, DERtemp] = MATDSS_DER(MATDSStemp);

    AllNodesNames = MATDSStemp.Sim.DSSCircuit.AllNodeNames;
    % Get all lines/branches names
    MyLines = MATDSStemp.Sim.DSSCircuit.Lines;
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

end