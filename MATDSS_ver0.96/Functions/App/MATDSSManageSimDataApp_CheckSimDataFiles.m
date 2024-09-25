function MATDSSManageSimDataApp_CheckSimDataFiles(app)
% MATDSSManageSimDataApp_CheckSimDataFiles(app)
% This function reads all available Simulation Data files in the folder
% ".\Config\SimSetupMats" and displays them in the MATDSSManageSimDataApp window.
% The files are loaded, and their data is added to the app's SimDataListBox.
%
% Parameters:
%   - app: The app object used to store application data.
%
% Last Update for this function was on MATDSS App Ver 0.96 (20 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Define the folder where the simulation data files are stored
SimDataFolder = [pwd '\Config\SimSetupMats'];

% Use the OS 'dir' command to list all .mat files in the SimDataFolder and subdirectories
[~, SimDataFiles] = dos(['dir /s /b ' '"' fullfile(SimDataFolder, '*.mat') '"']);

% Convert the list of files into a cell array
SimDataFiles = textscan(SimDataFiles, '%s', 'delimiter', '\n');
SimDataFiles = SimDataFiles{:};

% Initialize arrays for storing loaded simulation data and file names
SimDataFilesData = [];
FilesNameList = {};

% Loop through each file, load its contents, and extract the SimSetupData structure
for i = 1:length(SimDataFiles)
    % Load the .mat file
    temp = load(SimDataFiles{i});

    % Concatenate the SimSetupData from each file into a single array
    SimDataFilesData = [SimDataFilesData; temp.SimSetupData];

    % Create a name for each file combining its name and the current time field
    FilesNameList = [FilesNameList; {[temp.SimSetupData.Name '_' temp.SimSetupData.currenttime]}];
end

% Update the app's SimDataListBox with the names of the simulation data files
app.SimDataListBox.Items = FilesNameList;

% Store the loaded data and full file paths in the app properties
app.SimDataFilesData = SimDataFilesData;
app.FilesNameList = FilesNameList;
app.SimDataFilesFullAddress = SimDataFiles;

% If the SimDataListBox is empty, but there are items, set the default selected item
if isempty(app.SimDataListBox.Value) && ~isempty(app.SimDataListBox.Items)
    app.SimDataListBox.Value = app.SimDataListBox.Items(1);

    % Call the function to handle any updates required when the selection changes
    MATDSSManageSimDataApp_SimDataListBoxValueChanged(app);
end

end


%% Old function code

%{

function MATDSSManageSimDataApp_CheckSimDataFiles(app)
% This function reads all avaiable Simulation Data files in the folder
% \Config\SimSetupMats and display them in MATDSSManageSimDataApp window



% check all files in the folder
SimDataFolder = [pwd '\Config\SimSetupMats'];
[~,SimDataFiles] = dos(['dir /s /b ' '"' fullfile(SimDataFolder,'*.mat') '"' ]);  % OS dir command
SimDataFiles=textscan(SimDataFiles,'%s','delimiter','\n'); SimDataFiles=SimDataFiles{:};
SimDataFilesData = [];
FilesNameList = {};
for i = 1:length(SimDataFiles)
    temp = load(SimDataFiles{i});
    SimDataFilesData = [SimDataFilesData;temp.SimSetupData];
    FilesNameList = [FilesNameList; {[temp.SimSetupData.Name '_' temp.SimSetupData.currenttime]}];
end



app.SimDataListBox.Items = FilesNameList;
app.SimDataFilesData = SimDataFilesData;
app.FilesNameList = FilesNameList;
app.SimDataFilesFullAddress = SimDataFiles;

if isempty(app.SimDataListBox.Value) && ~isempty(app.SimDataListBox.Items)
    app.SimDataListBox.Value = app.SimDataListBox.Items(1);
    MATDSSManageSimDataApp_SimDataListBoxValueChanged(app);
end


end

%}