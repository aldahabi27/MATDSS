function MATDSSManageSimDataApp_SimDataListBoxValueChanged(app)
% MATDSSManageSimDataApp_SimDataListBoxValueChanged(app)
% This function displays the content of the selected simulation configuration's
% SimDataFile in the text area of the MATDSS application.
%
% Parameters:
%   - app: The MATDSS application instance.
%
% Last Update for this function was on MATDSS App Ver 0.96 (20 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Find the index of the selected simulation data file
iSimData = find(strcmpi([app.SimDataListBox.Items], app.SimDataListBox.Value));

% Display summary of the configuration in the text area
app.SimDataTextArea.Value = ['OpenDSS File = ' app.SimDataFilesData(iSimData).Name];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['Date & Time = ' app.SimDataFilesData(iSimData).currenttime]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nDERs = ' num2str(length(app.SimDataFilesData(iSimData).DER))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nCAs = ' num2str(length(app.SimDataFilesData(iSimData).MATDSS.Cont.CA))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['Global Gain = ' num2str(app.SimDataFilesData(iSimData).MATDSS.Cont.Gain)]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nV_buses = ' num2str(length(app.SimDataFilesData(iSimData).Table.VTrackingListBoxValues))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nI_Lines = ' num2str(length(app.SimDataFilesData(iSimData).Table.ITrackingListBoxValues))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line

% Display DERs Table
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['DERs Table']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.DER)];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line

% Display Controllers Settings
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['ControllersSettings']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.ControllersSettings)];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line

% Display VDERs
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['VDERs']];
if isempty(app.SimDataFilesData(iSimData).Table.VDERSettings)
    app.SimDataTextArea.Value = [app.SimDataTextArea.Value; 'No VDERs Defined'];
else
    app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.VDERSettings)];
end
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line

% Display Control Areas
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['ControlAreas']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.CASettings)];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']]; % Empty line

% Display Controlled Buses and Lines
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['Controlled Buses and Lines']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(createTable(app.SimDataFilesData(iSimData).Table.VTrackingListBoxValues', app.SimDataFilesData(iSimData).Table.ITrackingListBoxValues'))];

end

% Helper function to format tables for display
function formattedText = formatTable(tableVariable)
% Get the number of columns in the table
numCols = size(tableVariable, 2);

% Initialize an array to store the maximum width of each column
colWidths = zeros(1, numCols);

% Iterate over each column to find the maximum width of cell values
for col = 1:numCols
    colWidths(col) = max(cellfun(@(x) numel(char(x)), tableVariable{:, col}));
end

% Initialize an array to store the maximum width of each header
headerWidths = cellfun(@(x) numel(char(x)), tableVariable.Properties.VariableNames);

% Calculate the maximum width for each column
colWidths = max(colWidths, headerWidths);

% Initialize an empty string to store the formatted text
formattedText = '';

% Format and append the table headers
for col = 1:numCols
    header = char(tableVariable.Properties.VariableNames{col});
    numSpacesBefore = floor((colWidths(col) - numel(header)) / 2);
    numSpacesAfter = colWidths(col) - numel(header) - numSpacesBefore;
    formattedText = [formattedText, repmat(' ', 1, numSpacesBefore), header, repmat(' ', 1, numSpacesAfter), '   '];
end
formattedText = [formattedText, newline];

% Add separator line
formattedText = [formattedText, repmat('-', 1, sum(colWidths) + 3 * numCols), newline];

% Iterate over each row in the table
for row = 1:size(tableVariable, 1)
    % Iterate over each column in the table
    for col = 1:numCols
        % Get the cell value as a string
        cellValue = char(tableVariable{row, col});

        % Calculate the number of spaces needed before and after the cell value
        numSpacesBefore = floor((colWidths(col) - numel(cellValue)) / 2);
        numSpacesAfter = colWidths(col) - numel(cellValue) - numSpacesBefore;

        % Append the required spaces, the cell value, and 3 additional spaces to the formatted text
        formattedText = [formattedText, repmat(' ', 1, numSpacesBefore), cellValue, repmat(' ', 1, numSpacesAfter), '   '];
    end

    % Add a new line character after each row
    formattedText = [formattedText, newline];
end
end

% Helper function to create a table from V_buses and I_Lines
function tableResult = createTable(V_buses, I_Lines)
numRows = max(length(V_buses), length(I_Lines));

% Pad the shorter cell array with empty strings
if length(V_buses) < numRows
    V_buses = [V_buses(:); repmat({''}, numRows - length(V_buses), 1)];
end

if length(I_Lines) < numRows
    I_Lines = [I_Lines(:); repmat({''}, numRows - length(I_Lines), 1)];
end

% Create the table
tableResult = table(V_buses, I_Lines, 'VariableNames', {'V_buses', 'I_Lines'});
end


%% Old Function code

%{

function MATDSSManageSimDataApp_SimDataListBoxValueChanged(app)
% This function will display the content of the SimDataFile of the selected
% simulation configuration, in the text area.

% Which data file is selected?
iSimData = find(strcmpi([app.SimDataListBox.Items],app.SimDataListBox.Value));


% Display summary of the configruation
app.SimDataTextArea.Value = ['OpenDSS File = ' app.SimDataFilesData(iSimData).Name];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['Date & Time = ' app.SimDataFilesData(iSimData).currenttime]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nDERs = ' num2str(length(app.SimDataFilesData(iSimData).DER))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nCAs = ' num2str(length(app.SimDataFilesData(iSimData).MATDSS.Cont.CA))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['Global Gain = ' num2str(app.SimDataFilesData(iSimData).MATDSS.Cont.Gain)]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nV_buses = ' num2str(length(app.SimDataFilesData(iSimData).Table.VTrackingListBoxValues))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['nI_Lines = ' num2str(length(app.SimDataFilesData(iSimData).Table.ITrackingListBoxValues))]];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['DERs Table']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.DER)];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['ControllersSettings']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.ControllersSettings)];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['VDERs']];
if isempty(app.SimDataFilesData(iSimData).Table.VDERSettings)
    app.SimDataTextArea.Value = [app.SimDataTextArea.Value; 'No VDERs Defined'];
else
    app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.VDERSettings)];    
end
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['ControlAreas']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(app.SimDataFilesData(iSimData).Table.CASettings)];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; [' ']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; ['Controlled Buses and Lines']];
app.SimDataTextArea.Value = [app.SimDataTextArea.Value; formatTable(createTable(app.SimDataFilesData(iSimData).Table.VTrackingListBoxValues', app.SimDataFilesData(iSimData).Table.ITrackingListBoxValues'))];

end


% CHAT GPT GENERATED FUNCTION
function formattedText = formatTable(tableVariable)
    % Get the number of columns in the table
    numCols = size(tableVariable, 2);
    
    % Initialize an array to store the maximum width of each column
    colWidths = zeros(1, numCols);
    
    % Iterate over each column to find the maximum width of cell values
    for col = 1:numCols
        colWidths(col) = max(cellfun(@(x) numel(char(x)), tableVariable{:, col}));
    end
    
    % Initialize an array to store the maximum width of each header
    headerWidths = cellfun(@(x) numel(char(x)), tableVariable.Properties.VariableNames);
    
    % Calculate the maximum width for each column
    colWidths = max(colWidths, headerWidths);
    
    % Initialize an empty string to store the formatted text
    formattedText = '';
    
    % Format and append the table headers
    for col = 1:numCols
        header = char(tableVariable.Properties.VariableNames{col});
        numSpacesBefore = floor((colWidths(col) - numel(header)) / 2);
        numSpacesAfter = colWidths(col) - numel(header) - numSpacesBefore;
        formattedText = [formattedText, repmat(' ', 1, numSpacesBefore), header, repmat(' ', 1, numSpacesAfter), '   '];
    end
    formattedText = [formattedText, newline];
    
    % Add separator line
    formattedText = [formattedText, repmat('-', 1, sum(colWidths) + 3 * numCols), newline];
    
    % Iterate over each row in the table
    for row = 1:size(tableVariable, 1)
        % Iterate over each column in the table
        for col = 1:numCols
            % Get the cell value as a string
            cellValue = char(tableVariable{row, col});
            
            % Calculate the number of spaces needed before and after the cell value
            numSpacesBefore = floor((colWidths(col) - numel(cellValue)) / 2);
            numSpacesAfter = colWidths(col) - numel(cellValue) - numSpacesBefore;
            
            % Append the required spaces, the cell value, and 3 additional spaces to the formatted text
            formattedText = [formattedText, repmat(' ', 1, numSpacesBefore), cellValue, repmat(' ', 1, numSpacesAfter), '   '];
        end
        
        % Add a new line character after each row
        formattedText = [formattedText, newline];
    end
end

function tableResult = createTable(V_buses, I_Lines)
    numRows = max(length(V_buses), length(I_Lines));
    
    % Pad the shorter cell array with empty strings
    if length(V_buses) < numRows
        V_buses = [V_buses(:); repmat({''}, numRows - length(V_buses), 1)];
    end
    
    if length(I_Lines) < numRows
        I_Lines = [I_Lines(:); repmat({''}, numRows - length(I_Lines), 1)];
    end
    
    % Create the table
    tableResult = table(V_buses, I_Lines, 'VariableNames', {'V_buses', 'I_Lines'});
end




%}