function MATDSSApp_UpdateCASettings(app)
% This function will update the tables in the UI when the user define or
% change the defined control areas buses

% Get a list of control areas defined by the user
CANumbers = str2double(app.ControlAreasTable.Data(:,3).Variables);
UniqueCANumbers = unique(CANumbers);

% app.NCALabel.Text = num2str(length(UniqueCANumbers));
ControllersSettingsTable = app.ControllersSettingsTable.Data;
VDERsTable = app.VDERsTable.Data;
NumControllersSettings = str2double(app.ControllersSettingsTable.Data(:,1).Variables);
for i = 1:length(UniqueCANumbers)
    CAIndex = find(UniqueCANumbers(i) == NumControllersSettings);
    if isempty(CAIndex)
        CSIndex = 1;
        for j = 1:length(NumControllersSettings)
            if UniqueCANumbers(i) > NumControllersSettings(j)
                CSIndex = j+1;
            end
        end
        NCSEntry = [num2str(UniqueCANumbers(i));"LLC";"auto";"1e-4";"1e-3";"1e2";"1.05";"0.95";"Specified in DSS File";...
            "1";"1"; "1e3"; "1e3"; "1e3"; "1e3"; "1e12"; "1e12"; "1e8";...
            "1";"1";"1e-3";"1e-3";"1e-3";"1e-3";"1e-12";"1e-12";"1e-8"];
        ControllersSettingsTable = [ControllersSettingsTable; cellstr(NCSEntry')];



        % Get buses information
        if UniqueCANumbers(i) > 1 % not CA1, then add VDER configs.
            VDEREntry = [num2str(UniqueCANumbers(i));...
                string(['VDER_CA' num2str(UniqueCANumbers(i))]);... DER Name
                '1e-6';...                                                                              Tau
                "Wye";...                                                                               Type of connection (for now will go with Wye, but we might need to have a smarter way of selecting why or delta)
                "DSS_load";...                                                                          DER Mode
                '[10,0,0]';...                                                                          P(x)
                '[2,0,0]';...                                                                          Q(x)
                '1'; '1'];
            VDERsTable = [VDERsTable; cellstr(VDEREntry')];
        end
    end
end
app.ControllersSettingsTable.Data = ControllersSettingsTable;
app.VDERsTable.Data = VDERsTable;

for j = 1:length(NumControllersSettings)
    CAIndex = find(UniqueCANumbers == NumControllersSettings(j));
    if isempty(CAIndex)
        ControllersSettingsTable = app.ControllersSettingsTable.Data;
        ControllersSettingsTable(j,:) = [];
        app.ControllersSettingsTable.Data = ControllersSettingsTable;

        VDERsTable = app.VDERsTable.Data;
        if j ~= 1
            VDERsTable(j-1,:) = [];
            app.VDERsTable.Data = VDERsTable;
        end
        break;
    end
end
end