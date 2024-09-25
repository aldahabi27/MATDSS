function MATDSSApp_Export(app)
% MATDSSApp_Export(app)
% This function exports figures and saves them to a folder inside
% ".\Exports\". The function will export each plot in the following formats:
%  - fig
%  - png
%  - pdf (not working yet)
%  - eps (not working yet)
%
% Additionally, the function will export the MATDSS run file ".mat".
%
%   Last Update for this function was on MATDSS App Ver 0.96 (18 Sept. 2024)
%
%   MATDSS Application
%   Copyright (c) 2024, Ilyas Farhat
%   by Ilyas Farhat
%
%   This file is part of MATDSS Application
%   Contact the developer at ilyas.farhat@outlook.com

% Extract necessary data from the MATDSS application instance
MATDSS = app.MyRun.MATDSS;   % The main MATDSS run object
DER = app.MyRun.DER;         % DER (Distributed Energy Resources) data
ContE = MATDSS.Cont.CA(1).E; % Control Area data

% Select the first tab in both TabGroup and DetailsDERSTabGroup to reset the UI
app.TabGroup.SelectedTab = app.TabGroup.Children(1); % Reset the main tab
app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(1); % Reset DER details tab

% Update the status label with a message indicating the export process has started
MATDSSApp_Status(app, 'Exporting Results and Simulation Data', 'Exporting in Progress, Please wait!');

try
    % Preparing Export Folder
    % Get the current directory path of the script and prepare to create the "Exports" folder
    OutputFileDir = mfilename('fullpath');
    bi = strfind(OutputFileDir,'\'); % Locate folder separators
    OutputFileDir = [OutputFileDir(1:bi(end-2)) 'Exports']; % Define the "Exports" folder path
    mkdir(OutputFileDir); % Create the "Exports" directory if it doesn't exist

    % Prompt the user to input a custom name if the custom checkbox is checked
    if app.CustomCheckBox.Value
        UserCustomName = convertStringsToChars(string(inputdlg("What would you like to add to the name?")));
    else
        UserCustomName = []; % Set to empty if no custom name is provided
    end

    % Determine the output folder name based on user input or default configuration
    if isempty(UserCustomName)
        OutputFolder = app.MATDSSExcelSheet(1:end-4); % Default name based on Excel sheet
    else
        OutputFolder = [app.MATDSSExcelSheet(1:end-4) '_' UserCustomName]; % Add custom name if provided
    end

    % Append current date and time to the folder name for uniqueness
    currenttime = datetime;
    currenttime.Format = 'MMddyyHHmm';
    currenttime = char(currenttime); % Convert datetime to char
    OutputFileDir = [OutputFileDir '\' OutputFolder '_' currenttime]; % Final folder path
    mkdir(OutputFileDir); % Create the final output folder
    mkdir([OutputFileDir '\pdf']); % Create subdirectory for PDFs
    mkdir([OutputFileDir '\eps']); % Create subdirectory for EPS files
    mkdir([OutputFileDir '\fig']); % Create subdirectory for FIG files



    % Exporting mat file
    MATDSSApp_Details(app, ['Exporting "' app.MATDSSExcelSheet(1:end-4) '_' currenttime '.dat" file']);
    % Update the details text area, informing the user that the .dat file is being exported.
    % This will show the name of the .dat file based on the current date and time.

    myData = {}; % Initialize an empty cell array to store app data for export
    myDatafieldnames = string(fieldnames(app)); % Get all the field names (properties) of the app object

    % Loop through each field in the app object
    for i = 1:length(myDatafieldnames)
        ClassName = class(app.(myDatafieldnames(i))); % Get the class type of the current field
        idot = strfind(ClassName, '.'); % Check if there is a '.' in the class name (for nested or custom classes)

        % Determine the class type name (after the last dot, if applicable)
        if ~isempty(idot)
            Type = ClassName(idot(end)+1:end);
        else
            Type = ClassName;
        end

        % Switch-case to handle different class types and store relevant data in 'myData'
        switch lower(Type)
            case 'figure' % Handle figure objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder as figures are not exported

            case 'table' % Handle table objects
                % If the table has a 'Data' field, export the table's data
                if find(ismember(lower(fieldnames(app.(myDatafieldnames(i)))), 'data'))
                    myData = [myData; {myDatafieldnames(i)}, {{app.(myDatafieldnames(i)).Data}}];
                else % If no 'Data' field exists, store the whole table object
                    myData = [myData; {myDatafieldnames(i)}, {{app.(myDatafieldnames(i))}}];
                end

            case 'textarea' % Handle textarea objects
                % Store the value of the textarea field as a cell string
                myData = [myData; {myDatafieldnames(i)}, {cellstr(app.(myDatafieldnames(i)).Value)}];

            case 'checkboxtree' % Handle checkboxtree objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder, no specific export for checkboxtree

            case 'button' % Handle button objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder, buttons do not hold exportable data

            case 'statebutton' % Handle statebutton objects
                % Store the button's state (value) in the exported data
                myData = [myData; {myDatafieldnames(i)}, {app.(myDatafieldnames(i)).Value}];

            case 'editfield' % Handle editfield objects
                % Store the value entered in the edit field
                myData = [myData; {myDatafieldnames(i)}, {app.(myDatafieldnames(i)).Value}];

            case 'struct' % Handle struct objects
                % Store the entire structure in the exported data
                myData = [myData; {myDatafieldnames(i)}, {app.(myDatafieldnames(i))}];

            case 'menu' % Handle menu objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder for menu objects

            case 'gridlayout' % Handle gridlayout objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder for gridlayout objects

            case 'tabgroup' % Handle tabgroup objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder for tabgroup objects

            case 'tab' % Handle tab objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder for tab objects

            case 'panel' % Handle panel objects
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Placeholder for panel objects

            otherwise % Handle any other unlisted object types
                myData = [myData; {myDatafieldnames(i)}, {""}]; % Default placeholder for unhandled object types
        end

        % Future additional processing or adjustments to export logic can go here

    end


    % Write the collected app data to a .dat file
    writetable(cell2table(myData), [OutputFileDir '\' app.MATDSSExcelSheet(1:end-4) '_' currenttime '.dat']);
    % Convert the 'myData' cell array to a table and write it to a .dat file in the specified directory.
    % The filename is based on the Excel sheet name and the current time.

    % Update details text area about exporting screenshots
    MATDSSApp_Details(app, ['Exporting screenshots for configurations']);

    % Loop through each tab in the CircuitConfigurationsTabGroup to capture and export screenshots
    for j = 1:size(app.CircuitConfigurationsTabGroup.Children,1)
        app.TabGroup.SelectedTab = app.TabGroup.Children(2); % Switch to the second tab in the main TabGroup
        app.CircuitConfigurationsTabGroup.SelectedTab = app.CircuitConfigurationsTabGroup.Children(j); % Switch to the specific tab in the configurations TabGroup

        pause(0.5); % Short pause to allow the UI to update

        screens = captureScreen; % Capture the current screen (returns a cell array of screenshots)

        % Save each captured screen to an image file
        for i = 1:length(screens)
            imwrite(screens{i}, [OutputFileDir '\' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'SC' num2str(i) ...
                '_' app.CircuitConfigurationsTabGroup.Children(j).Title  '.png']);
            % Save the screenshot as a .png file with a unique name based on the tab title, time, and index.
        end
    end

    % Reset the selected tabs back to their initial state
    app.TabGroup.SelectedTab = app.TabGroup.Children(1); % Switch back to the first tab in the main TabGroup
    app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(1); % Switch back to the first tab in the DERSTabGroup


    % Exporting plot data to a .mat file
    MATDSSApp_Details(app, ['Exporting "' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'PlotsData.mat" file']);

    % Initialize an empty cell array for the plot data
    MyPlotsData = {};

    % Export active power (P) measurement and simulation data
    MyPlotsData = [MyPlotsData;...
        {[char(916) 'P0,meas']}, MATDSS.Meas.at - MATDSS.Time.Sim.ST, (sum(MATDSS.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,sim']}, MATDSS.Sim.Meas.at - MATDSS.Time.Sim.ST, (sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,set-meas']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,set-sim']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,set-ll']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3;...
        {[char(916) 'P0,set-ul']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.P0Set + MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3;...
        {[char(916) 'P0,set-Ell']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3;...
        {[char(916) 'P0,set-Eul']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3];

    % Adding DER power curves
    for i = 1:size(DER, 2)
        MyPlotsData = [MyPlotsData;...
            {['P_' DER(i).DSSName]}, MATDSS.Sim.Meas.at - MATDSS.Time.Sim.ST, DER(i).W(1:end-1)./1e3];
        % Export each DER's power (W) as a curve, converted from W to kW.
    end

    % Check if reactive power (Q) data exists in ControlSignals and export if available
    if isfield(MATDSS.ControlSignals, 'Q0Set')
        MyPlotsData = [MyPlotsData;...
            {[char(916) 'Q0,meas']}, MATDSS.Meas.at - MATDSS.Time.Sim.ST, (sum(MATDSS.Meas.Q0) - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,sim']}, MATDSS.Sim.Meas.at - MATDSS.Time.Sim.ST, (sum(MATDSS.Sim.Meas.Q0) - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,set-meas']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,set-sim']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,set-ll']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1) - 10.*ContE)./1e3;...
            {[char(916) 'Q0,set-ul']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.Q0Set + MATDSS.ControlSignals.Q0Set(1) - 10.*ContE)./1e3;...
            {[char(916) 'Q0,set-Ell']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1) - ContE)./1e3;...
            {[char(916) 'Q0,set-Eul']}, MATDSS.Time.Sim.TimeSpan - MATDSS.Time.Sim.ST, (MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1) + ContE)./1e3];
        % Export reactive power (Q) data if it exists, similar to the active power (P) data.
    end


    % Adding DER Q curves
    for i = 1:size(DER, 2)
        % Store reactive power (Q) data for each DER in MyPlotsData
        MyPlotsData = [MyPlotsData; {['Q_' DER(i).DSSName]}, ...
            MATDSS.Sim.Meas.at - MATDSS.Time.Sim.ST, ...
            DER(i).Var(1:end-1)./1e3];
    end

    % Loop through Control Areas (CAs) and store P and Q data
    for i = 1:size(MATDSS.Cont.CA, 1)
        CA = MATDSS.Cont.CA(i);

        % Store real power (P) data for each Control Area in MyPlotsData
        MyPlotsData = [MyPlotsData; ...
            {[char(916) 'P0_CA' num2str(CA.Area)]}, ...
            MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
            (sum(CA.P0) - sum(CA.P0(:, 1)))./1e3; ...
            {[char(916) 'P0,set_CA' num2str(CA.Area)]}, ...
            MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
            (CA.P0Set - CA.P0Set(1))./1e3; ...
            {[char(916) 'P0,set-Ell_CA' num2str(CA.Area)]}, ...
            MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
            (CA.P0Set - CA.P0Set(1) - CA.E)./1e3; ...
            {[char(916) 'P0,set-Eul_CA' num2str(CA.Area)]}, ...
            MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
            (CA.P0Set - CA.P0Set(1) + CA.E)./1e3];

        % Store reactive power (Q) data for each Control Area
        MyPlotsData = [MyPlotsData; ...
            {[char(916) 'Q0_CA' num2str(CA.Area)]}, ...
            MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
            (sum(CA.Q0) - sum(CA.Q0(:, 1)))./1e3];

        % If reactive power control signals exist, store them in MyPlotsData
        if isfield(MATDSS.ControlSignals, 'Q0Set')
            MyPlotsData = [MyPlotsData; ...
                {[char(916) 'Q0,set_CA' num2str(CA.Area)]}, ...
                MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
                (CA.Q0Set - CA.Q0Set(1))./1e3; ...
                {[char(916) 'Q0,set-Ell_CA' num2str(CA.Area)]}, ...
                MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
                (CA.Q0Set - CA.Q0Set(1) - CA.E)./1e3; ...
                {[char(916) 'Q0,set-Eul_CA' num2str(CA.Area)]}, ...
                MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
                (CA.Q0Set - CA.Q0Set(1) + CA.E)./1e3];
        end
    end

    % Extract DER and node/branch names from MATDSS structures
    AllDERNames = {DER.DSSName};
    k_v = MATDSS.Meas.k_v;
    AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
    AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

    k_i = MATDSS.Meas.k_i;
    AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
    AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);

    % Store voltage magnitude profile (measured and simulated) for each node
    for i = 1:size(AllMeasNodeNames, 1)
        MyPlotsData = [MyPlotsData; ...
            {['V_meas_' convertStringsToChars(AllMeasNodeNames{i})]}, ...
            MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
            MATDSS.Meas.VMagProfilePu(i, :)];
    end

    for i = 1:size(AllNodeNames, 1)
        MyPlotsData = [MyPlotsData; ...
            {['V_sim_' convertStringsToChars(AllNodeNames{i})]}, ...
            MATDSS.Sim.Meas.at - MATDSS.Time.Sim.ST, ...
            MATDSS.Sim.Meas.VMagProfilePu(i, :)];
    end

    % Store current profile (measured and simulated) for each branch
    for i = 1:size(AllMeasBranchNames, 1)
        MyPlotsData = [MyPlotsData; ...
            {['I_meas_' convertStringsToChars(AllMeasBranchNames{i})]}, ...
            MATDSS.Meas.at - MATDSS.Time.Sim.ST, ...
            MATDSS.Meas.IProfile(i, :)];
    end

    for i = 1:size(AllBranchNames, 1)
        MyPlotsData = [MyPlotsData; ...
            {['I_sim_' convertStringsToChars(AllBranchNames{i})]}, ...
            MATDSS.Sim.Meas.at - MATDSS.Time.Sim.ST, ...
            MATDSS.Sim.Meas.IProfile(i, :)];
    end

    % Backup configuration data
    ConfigDataBackup = MyPlotsData;
    MyPlotsData = [];
    save([OutputFileDir '\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'ConfigDataBackup.dat'], "ConfigDataBackup")

    % Prepare to save plot data
    MyPlotsData = app.MyRun;

    % Clear unneeded DSS simulation objects from MyPlotsData
    MyPlotsData.MATDSS.Sim.DSSObj = [];
    MyPlotsData.MATDSS.Sim.DSSText = [];
    MyPlotsData.MATDSS.Sim.DSSCircuit = [];
    MyPlotsData.MATDSS.Sim.DSSSolution = [];

    % Pause for 1 second to ensure data is saved properly
    pause(1);

    % Save plot data in a .mat file
    save([OutputFileDir '\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'PlotsData.mat'], '-struct', 'MyPlotsData', '-v7.3');


    % Copying configuration files for reference
    MATDSSApp_Details(app, "Copying configuration files (excel files)");

    % Create a directory for configuration files backup
    mkdir([OutputFileDir '\Config']);

    % Check if the MATDSS_ControlAreas.xlsx file exists and copy it to the backup directory
    if exist([pwd '\Config\MATDSS_ControlAreas.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_ControlAreas.xlsx'], ...
            [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_ControlAreas_bkp.xlsx']);
    end

    % Check if the MATDSS_ControllersSettings.xlsx file exists and copy it to the backup directory
    if exist([pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_ControllersSettings.xlsx'], ...
            [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_ControllersSettings_bkp.xlsx']);
    end

    % Check if the MATDSS_VI.xlsx file exists and copy it to the backup directory
    if exist([pwd '\Config\MATDSS_VI.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_VI.xlsx'], ...
            [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_VI_bkp.xlsx']);
    end

    % Check if the MATDSS_DER.xlsx file exists and copy it to the backup directory
    if exist([pwd '\Config\MATDSS_DER.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_DER.xlsx'], ...
            [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_DER_bkp.xlsx']);
    end

    % Check if the MATDSS_VDERs.xlsx file exists and copy it to the backup directory
    if exist([pwd '\Config\MATDSS_VDERs.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_VDERs.xlsx'], ...
            [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_VDERs_bkp.xlsx']);
    end

    % Export figures from each tab in the PlotTabGroup
    for i = 1:size(app.PlotTabGroup.Children, 1)
        % Create a new figure for exporting with specified properties
        ExpFig = figure('Name', 'ExpFig', 'Visible', 'off');
        set(ExpFig, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);
        set(ExpFig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 12, 7], 'Visible', 'off');

        % Switch case to handle different plot panels
        switch i
            case 1
                % Copy the last child of the P0Panel to the new figure, or create new axes if empty
                if ~isempty(app.P0Panel.Children)
                    ExpAxes = copyobj(app.P0Panel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.P0Tiles.Children(2), ExpFig);
            case 2
                % Copy the last child of the CA_P0Panel to the new figure, or create new axes if empty
                if ~isempty(app.CA_P0Panel.Children)
                    ExpAxes = copyobj(app.CA_P0Panel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.CA_P0Tiles.Children(2), ExpFig);
            case 3
                % Copy the last child of the P_DERPanel to the new figure, or create new axes if empty
                if ~isempty(app.P_DERPanel.Children)
                    ExpAxes = copyobj(app.P_DERPanel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
            case 4
                % Copy the last child of the Q0Panel to the new figure, or create new axes if empty
                if ~isempty(app.Q0Panel.Children)
                    ExpAxes = copyobj(app.Q0Panel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.Q0Tiles.Children(2), ExpFig);
            case 5
                % Copy the last child of the CA_Q0Panel to the new figure, or create new axes if empty
                if ~isempty(app.CA_Q0Panel.Children)
                    ExpAxes = copyobj(app.CA_Q0Panel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.CA_Q0Tiles.Children(2), ExpFig);
            case 6
                % Copy the last child of the Q_DERPanel to the new figure, or create new axes if empty
                if ~isempty(app.Q_DERPanel.Children)
                    ExpAxes = copyobj(app.Q_DERPanel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.Q_DERTiles.Children(2), ExpFig);
            case 7
                % Copy the last child of the VoltagePanel to the new figure, or create new axes if empty
                if ~isempty(app.VoltagePanel.Children)
                    ExpAxes = copyobj(app.VoltagePanel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.VoltageTiles.Children(2), ExpFig);
            case 8
                % Copy the last child of the CurrentPanel to the new figure, or create new axes if empty
                if ~isempty(app.CurrentPanel.Children)
                    ExpAxes = copyobj(app.CurrentPanel.Children(end), ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.CurrentTiles.Children(2), ExpFig);
            otherwise
                break; % Exit the loop if the case is not handled
        end



        % Exporting plots and saving them in different formats
        MATDSSApp_Details(app, ['Exporting plot ' num2str(i) ' of ' num2str(size(app.PlotTabGroup.Children,1)-1)]);

        % Set the title of the figure to include the Excel sheet name and current time
        title([app.MATDSSExcelSheet(1:end-4) '_' currenttime], 'Interpreter', 'none');

        % Make the legend visible
        legend('Visible', 'on');

        % Print the figure to a PNG file with a resolution of 300 DPI
        print(ExpFig, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png']);

        % Ensure the figure is visible before saving
        set(ExpFig, 'CreateFcn', 'set(gcf,''Visible'',''on'')');

        % Save the figure as a .fig file
        saveas(ExpFig, [OutputFileDir '\fig\' app.PlotTabGroup.Children(i).Title '.fig']);

        % Save the figure as an .eps file
        saveas(ExpFig, [OutputFileDir '\eps\' app.PlotTabGroup.Children(i).Title '.eps']);

        % Set the figure orientation to landscape
        orient(ExpFig, 'landscape');

        % Save the figure as a PDF file
        saveas(ExpFig, [OutputFileDir '\pdf\' app.PlotTabGroup.Children(i).Title '.pdf']);

        % Close the figure to free up memory
        close(ExpFig);

    end

    % Display a message indicating that the export is completed successfully
    MATDSSApp_Details(app, 'Export is completed successfully.');

    % Update the status label to indicate the export is complete
    MATDSSApp_Status(app, 'Export Complete!');

catch MATDSS_Error
    % Initialize variables for capturing error details
    ErrorString = [];
    ErrorStack = MATDSS_Error.stack;

    % Loop through the error stack to compile error details
    for i = 1:length(ErrorStack)
        ErrorFiles = ErrorStack(i).file;
        ErrorLines = ErrorStack(i).line;

        % Concatenate error file and line number information
        ErrorString = [ErrorString; strcat(ErrorFiles, " -> Line ", num2str(ErrorLines))];
    end

    % Format the error message to include identifier, message, and stack trace
    ErrorMsg = sprintf('%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', ...
        MATDSS_Error.identifier, MATDSS_Error.message, ErrorString);

    % Display the error message and additional information
    MATDSSApp_Details(app, {ErrorMsg; ' '; '******************************************************'; ' '});

    % Update the status label to indicate an error occurred
    app.StatusLabel.Text = 'Error occurred!';
    app.StatusLabel.FontColor = app.MyRun.Plot.Colors(12, :);

    % Display a message indicating that an error occurred in the export function
    MATDSSApp_Details(app, {'Error in MATDSSApp_Export function occurred. Check the details above for more information'; 'The program will stop the export.'});
end


%% Old function Code

%{

function MATDSSApp_Export(app)
% MATDSSApp_Export(app)
% This function exports figures and saves them to a folder inside
% ".\Exports\". The function will export each plot in the following formats:
%  - fig
%  - png
%  - pdf (not working yet)
%  - eps (not working yet)
%
% Additionally, the function will export the MATDSS run file ".mat".
%
%   Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
%   MATDSS Application
%   Copyright (c) 2024, Ilyas Farhat
%   by Ilyas Farhat
%
%   This file is part of MATDSS Application
%   Contact the developer at ilyas.farhat@outlook.com

MATDSS = app.MyRun.MATDSS;
DER = app.MyRun.DER;
ContE = MATDSS.Cont.CA(1).E;


app.TabGroup.SelectedTab = app.TabGroup.Children(1);
app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(1);

MATDSSApp_Status(app,'Exporting Results and Simulation Data', 'Exporting in Progress, Please wait!');
try
    % Preparing Export Folder
    OutputFileDir = mfilename('fullpath'); bi = strfind(OutputFileDir,'\');
    OutputFileDir = [OutputFileDir(1:bi(end-2)) 'Exports'];
    mkdir(OutputFileDir);
    if app.CustomCheckBox.Value
        UserCustomName = convertStringsToChars(string(inputdlg("What would you like to add to the name?")));
    else
        UserCustomName = [];
    end

    % if isempty(UserCustomName)
    %     return;
    % end
    if isempty(UserCustomName)
        OutputFolder = app.MATDSSExcelSheet(1:end-4);
    else
        OutputFolder = [app.MATDSSExcelSheet(1:end-4) '_' UserCustomName];
    end
    currenttime = datetime;
    currenttime.Format = 'MMddyyHHmm';
    currenttime = char(currenttime);
    OutputFileDir = [OutputFileDir '\' OutputFolder '_' currenttime];
    mkdir(OutputFileDir);
    mkdir([OutputFileDir '\pdf']);
    mkdir([OutputFileDir '\eps']);
    mkdir([OutputFileDir '\fig']);


 


    % Exporting mat file
    MATDSSApp_Details(app,['Exporting "' app.MATDSSExcelSheet(1:end-4) '_' currenttime '.dat" file']);
    myData = {};
    myDatafieldnames = string(fieldnames(app));

    for i = 1:length(myDatafieldnames)
        ClassName = class(app.(myDatafieldnames(i)));
        idot = strfind(ClassName,'.');
        if ~isempty(idot)
            Type = ClassName(idot(end)+1:end);
        else
            Type = ClassName;
        end
        switch lower(Type)
            case 'figure'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            case 'table'
                if find(ismember(lower(fieldnames(app.(myDatafieldnames(i)))),'data'))
                    myData = [myData;...
                        {myDatafieldnames(i)}, {{app.(myDatafieldnames(i)).Data}}];
                else
                    myData = [myData;...
                        {myDatafieldnames(i)}, {{app.(myDatafieldnames(i))}}];
                end
            case 'textarea'
                myData = [myData;...
                    {myDatafieldnames(i)}, {cellstr(app.(myDatafieldnames(i)).Value)}];
            case 'checkboxtree'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            case 'button'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            case 'statebutton'
                myData = [myData;...
                    {myDatafieldnames(i)}, {app.(myDatafieldnames(i)).Value}];
            case 'editfield'
                myData = [myData;...
                    {myDatafieldnames(i)}, {app.(myDatafieldnames(i)).Value}];
            case 'struct'
                myData = [myData;...
                    {myDatafieldnames(i)}, {app.(myDatafieldnames(i))}];
            case 'menu'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            case 'gridlayout'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            case 'tabgroup'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            case 'tab'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            case 'panel'
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
            otherwise
                myData = [myData;...
                    {myDatafieldnames(i)},{""}];% {app.(myDatafieldnames(i)).Value}];
        end


%{
               
        
%}
    end


    writetable(cell2table(myData), [OutputFileDir '\' app.MATDSSExcelSheet(1:end-4) '_' currenttime '.dat']);
    MATDSSApp_Details(app,['Exporting screenshots for configurations']);
    for j = 1:size(app.CircuitConfigurationsTabGroup.Children,1)
        app.TabGroup.SelectedTab = app.TabGroup.Children(2);
        app.CircuitConfigurationsTabGroup.SelectedTab = app.CircuitConfigurationsTabGroup.Children(j);
        pause(0.5);
        screens = captureScreen;
        for i = 1:length(screens)
            imwrite(screens{i},[OutputFileDir '\' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'SC' num2str(i) ...
                '_' app.CircuitConfigurationsTabGroup.Children(j).Title  '.png'])
        end
    end
    app.TabGroup.SelectedTab = app.TabGroup.Children(1);
    app.DetailsDERSTabGroup.SelectedTab = app.DetailsDERSTabGroup.Children(1);
    
    



MATDSSApp_Details(app,['Exporting "' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'PlotsData.mat" file']);
    % Exporting Plots Data to excel
    MyPlotsData = {};
    MyPlotsData = [MyPlotsData;...
        {[char(916) 'P0,meas']}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,sim']}, MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,set-meas']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,set-sim']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3;...
        {[char(916) 'P0,set-ll']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3;...
        {[char(916) 'P0,set-ul']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set + MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3;...
        {[char(916) 'P0,set-Ell']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3;...
        {[char(916) 'P0,set-Eul']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3];

    % Adding DER P curves
    for i = 1:size(DER,2)
        MyPlotsData = [MyPlotsData;...
        {['P_' DER(i).DSSName]}, MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(i).W(1:end-1)./1e3;];
    end
    
    if isfield(MATDSS.ControlSignals,'Q0Set')
        MyPlotsData = [MyPlotsData;...
            {[char(916) 'Q0,meas']}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.Q0) - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,sim']}, MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.Q0) - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,set-meas']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,set-sim']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1))./1e3;...
            {[char(916) 'Q0,set-ll']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1) - 10.*ContE)./1e3;...
            {[char(916) 'Q0,set-ul']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set + MATDSS.ControlSignals.Q0Set(1) - 10.*ContE)./1e3;...
            {[char(916) 'Q0,set-Ell']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1) - ContE)./1e3;...
            {[char(916) 'Q0,set-Eul']}, MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1) + ContE)./1e3];
    end
    % Adding DER Q curves
    for i = 1:size(DER,2)
        MyPlotsData = [MyPlotsData;...
        {['Q_' DER(i).DSSName]}, MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(i).Var(1:end-1)./1e3;];
    end

    for i = 1:size(MATDSS.Cont.CA,1)
        CA = MATDSS.Cont.CA(i);
        MyPlotsData = [MyPlotsData;...
            {[char(916) 'P0_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(CA.P0) - sum(CA.P0(:,1)))./1e3;...
            {[char(916) 'P0,set_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1))./1e3;...
            {[char(916) 'P0,set-Ell_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1) - CA.E)./1e3;...
            {[char(916) 'P0,set-Eul_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1) + CA.E)./1e3;];

        MyPlotsData = [MyPlotsData;...
            {[char(916) 'Q0_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(CA.Q0) - sum(CA.Q0(:,1)))./1e3;];
        if isfield(MATDSS.ControlSignals,'Q0Set')
            MyPlotsData = [MyPlotsData;...
                {[char(916) 'Q0,set_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1))./1e3;...
                {[char(916) 'Q0,set-Ell_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1) - CA.E)./1e3;...
                {[char(916) 'Q0,set-Eul_CA' num2str(CA.Area)]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1) + CA.E)./1e3;];
        end
    end

    AllDERNames = {DER.DSSName};
    k_v = MATDSS.Meas.k_v;
    AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
    AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

    k_i = MATDSS.Meas.k_i;
    AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
    AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);
    for i = 1:size(AllMeasNodeNames,1)
        MyPlotsData = [MyPlotsData;...
            {['V_meas_' convertStringsToChars(AllMeasNodeNames{i})]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.VMagProfilePu(i,:);];
    end

    for i = 1:size(AllNodeNames,1)
        MyPlotsData = [MyPlotsData;...
            {['V_sim_' convertStringsToChars(AllNodeNames{i})]}, MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.VMagProfilePu(i,:);];
    end

    for i = 1:size(AllMeasBranchNames,1)
        MyPlotsData = [MyPlotsData;...
            {['I_meas_' convertStringsToChars(AllMeasBranchNames{i})]}, MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.IProfile(i,:);];
    end

    for i = 1:size(AllBranchNames,1)
        MyPlotsData = [MyPlotsData;...
            {['I_sim_' convertStringsToChars(AllBranchNames{i})]}, MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.IProfile(i,:);];
    end
    
    ConfigDataBackup = MyPlotsData;
    MyPlotsData = [];
    save([OutputFileDir '\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'ConfigDataBackup.dat'],"ConfigDataBackup")
    MyPlotsData = app.MyRun;
    MyPlotsData.MATDSS.Sim.DSSObj = [];
    MyPlotsData.MATDSS.Sim.DSSText = [];
    MyPlotsData.MATDSS.Sim.DSSCircuit = [];
    MyPlotsData.MATDSS.Sim.DSSSolution = [];
    
    pause(1);
    save([OutputFileDir '\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'PlotsData.mat'],'-struct', 'MyPlotsData','-v7.3')
    
    

    % Copying config. files for refs.
    MATDSSApp_Details(app,"Copying configuration files (excel files)")
    mkdir([OutputFileDir '\Config']);

    if exist([pwd '\Config\MATDSS_ControlAreas.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_ControlAreas.xlsx'], [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_ControlAreas_bkp.xlsx'])
    end

    if exist([pwd '\Config\MATDSS_ControllersSettings.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_ControllersSettings.xlsx'], [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_ControllersSettings_bkp.xlsx'])
    end

    if exist([pwd '\Config\MATDSS_VI.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_VI.xlsx'], [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_VI_bkp.xlsx'])
    end

    if exist([pwd '\Config\MATDSS_DER.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_DER.xlsx'], [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_DER_bkp.xlsx'])
    end

    if exist([pwd '\Config\MATDSS_VDERs.xlsx'], 'file')
        copyfile([pwd '\Config\MATDSS_VDERs.xlsx'], [OutputFileDir '\Config\' UserCustomName '_' app.MATDSSExcelSheet(1:end-4) '_' currenttime 'MATDSS_VDERs_bkp.xlsx'])
    end

    for i = 1:size(app.PlotTabGroup.Children,1)
        ExpFig = figure('Name', 'ExpFig','Visible', 'off');
        set(ExpFig, 'Units', 'Normalized', 'Position', [0,0,1,1]);
        set(ExpFig, 'PaperUnits','inches', 'PaperPosition', [0,0,12,7], 'Visible', 'off');
        switch i
            case 1
                if ~isempty(app.P0Panel.Children)
                    ExpAxes = copyobj(app.P0Panel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.P0Tiles.Children(2),ExpFig);
            case 2
                if ~isempty(app.CA_P0Panel.Children)
                    ExpAxes = copyobj(app.CA_P0Panel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.CA_P0Tiles.Children(2),ExpFig);
            case 3
                if ~isempty(app.P_DERPanel.Children)
                    ExpAxes = copyobj(app.P_DERPanel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
            case 4
                 if ~isempty(app.Q0Panel.Children)
                    ExpAxes = copyobj(app.Q0Panel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.Q0Tiles.Children(2),ExpFig);
            case 5
                 if ~isempty(app.CA_Q0Panel.Children)
                    ExpAxes = copyobj(app.CA_Q0Panel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.CA_Q0Tiles.Children(2),ExpFig);
            case 6
                 if ~isempty(app.Q_DERPanel.Children)
                     ExpAxes = copyobj(app.Q_DERPanel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.Q_DERTiles.Children(2),ExpFig);
            case 7
                 if ~isempty(app.VoltagePanel.Children)
                    ExpAxes = copyobj(app.VoltagePanel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.VoltageTiles.Children(2),ExpFig);
            case 8
                 if ~isempty(app.CurrentPanel.Children)
                    ExpAxes = copyobj(app.CurrentPanel.Children(end),ExpFig);
                else
                    ExpAxes = axes(ExpFig);
                end
                % ExpAxes = copyobj(app.MyRun.Plot.PlotHandle.CurrentTiles.Children(2),ExpFig);
            otherwise
                break;
        end


        MATDSSApp_Details(app,['Exporting plot ' num2str(i) ' of ' num2str(size(app.PlotTabGroup.Children,1)-1)]);
        title([app.MATDSSExcelSheet(1:end-4) '_' currenttime], 'Interpreter','none');
        legend('Visible','on')
        % set(findall(ExpFig, '-property', 'FontSize'), 'FontSize', 16);
        print(ExpFig, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
        set(ExpFig, 'CreateFcn', 'set(gcf,''Visible'',''on'')'); 
        saveas(ExpFig, [OutputFileDir '\fig\' app.PlotTabGroup.Children(i).Title '.fig'])
        saveas(ExpFig, [OutputFileDir '\eps\' app.PlotTabGroup.Children(i).Title '.eps'])

        orient(ExpFig, 'landscape');
        saveas(ExpFig, [OutputFileDir '\pdf\' app.PlotTabGroup.Children(i).Title '.pdf'])

        close(ExpFig)


    end

    MATDSSApp_Details(app,'Export is completed successfully.');
    MATDSSApp_Status(app,'Export Complete!');

catch MATDSS_Error
    ErrorString = [];
    ErrorStack = MATDSS_Error.stack;
    for i = 1:length(ErrorStack)
        ErrorFiles = ErrorStack(i).file;
        ErrorLines = ErrorStack(i).line;

        ErrorString = [ErrorString;strcat(ErrorFiles," -> Line ", num2str(ErrorLines))];
    end
    ErrorMsg = sprintf('%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', MATDSS_Error.identifier,MATDSS_Error.message,ErrorString);
    %                 app.DetailsTextArea.Value = [app.DetailsTextArea.Value; ErrorMsg;' ';'******************************************************';' '];
    MATDSSApp_Details(app,{ErrorMsg;' ';'******************************************************';' '});
    app.StatusLabel.Text = 'Error occured!';
    app.StatusLabel.FontColor = app.MyRun.Plot.Colors(12,:);
    MATDSSApp_Details(app,{'Error in MATDSSApp_Export function occured. Check the details above for more information';'The program will stop the export.'});
end

%}