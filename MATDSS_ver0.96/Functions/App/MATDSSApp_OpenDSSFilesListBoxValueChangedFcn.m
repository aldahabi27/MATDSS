function MATDSSApp_OpenDSSFilesListBoxValueChangedFcn(app)
% MATDSSApp_OpenDSSFilesListBoxValueChangedFcn(app)
% This function will read the OpenDSS File selected and check the
% corresponding configuration files/spreadsheets.
% The function will notify the user if there are any major/fatal issues
% with the configurations.
%
%   Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
%   MATDSS Application
%   Copyright (c) 2024, Ilyas Farhat
%   by Ilyas Farhat
%
%   This file is part of MATDSS Application
%   Contact the developer at ilyas.farhat@outlook.com

if isempty(app.OpenDSSFilesListBox.Value)
    app.OpenDSSFilesListBox.Value = app.OpenDSSFilesListBox.Items(1);
end

if ~strcmpi(app.OpenDSSFilesListBox.Value, "(no OpenDSS file is located)")
    MATDSSApp_Status(app, ['Selecting ' app.OpenDSSFilesListBox.Value]);
    app.SimulateButton.Enable = true;
    DSSFileName = [app.OpenDSSFilesListBox.Tag app.OpenDSSFilesListBox.Value];

    % Open the new OpenDSS File in OpenDSS File tab
    MATDSSOpenDSSFileText = fileread(DSSFileName);
    app.MATDSSOpenDSSText.Value = MATDSSOpenDSSFileText;
    app.MATDSSExcelSheet = app.OpenDSSFilesListBox.Value; % Target sheet in the Excel file
    pause(0.1);

    % Save the OpenDSS File selected in MyRun variable for later use/reference
    app.MyRun.MATDSSOpenDSSFile.DSSFileName = DSSFileName;
    app.MyRun.MATDSSOpenDSSFile.Text = MATDSSOpenDSSFileText;

    MATDSSApp_Status(app, 'Ready', ['"' app.OpenDSSFilesListBox.Value '" Selected successfully!']);
end
