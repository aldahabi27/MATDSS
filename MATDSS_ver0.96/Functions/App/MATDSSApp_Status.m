function MATDSSApp_Status(app, StatusMsg, DetailsMsg, Per)
% MATDSSApp_Status(app, StatusMsg, DetailsMsg, Per)
% This function handles all status messages in the MATDSS Application.
% It allows updating the status information, optionally including details
% and percentage progress. The behavior varies based on the number and type 
% of input arguments provided.
%
% Parameters:
%   - app: The application object containing the status label to be updated.
%   - StatusMsg: The status message to display. Can be a string or numeric.
%   - DetailsMsg: Optional additional details message to display.
%   - Per: Optional percentage value for progress updates.
%
% MsgFlag values determine the type of update:
%   - 1: Update status text only.
%   - 2: Update percentage only.
%   - 3: Update status text and percentage.
%   - 4: Update status text and details message without percentage.
%   - 5: Update status text, percentage, and details message.
%
% Last Update: MATDSS App Ver 0.96 (12 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

MsgFlag = 5; % Default message flag

% Check the number of input arguments and set MsgFlag and Per accordingly
if nargin == 2
    if isnumeric(StatusMsg)
        Per = StatusMsg * 100; % Convert to percentage
        MsgFlag = 2; % Update only percentage
        DetailsMsg = ''; % No details message
    else
        MsgFlag = 1; % Update status text only
        DetailsMsg = ''; % No details message
        Per = -1; % No percentage
    end
end

if nargin == 3
    if ischar(DetailsMsg) || isstring(DetailsMsg)
        MsgFlag = 4; % Update status and details without percentage
        Per = -1; % No percentage
    elseif isnumeric(DetailsMsg)
        Per = DetailsMsg * 100; % Convert to percentage
        MsgFlag = 3; % Update status with percentage
    end
end

% Convert char inputs to string for consistency
if ischar(StatusMsg)
    StatusMsg = convertCharsToStrings(StatusMsg);
end

if ischar(DetailsMsg)
    DetailsMsg = convertCharsToStrings(DetailsMsg);
end

% Update the status label based on MsgFlag
switch MsgFlag
    case 1
        % Update only the status text
        % app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Text = StatusMsg;
        app.StatusLabel.Tag = StatusMsg;
    case 2
        % Update only the percentage
        % app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        Addpercentage(app, Per);
    case 3
        % Update status text and percentage
        % app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Tag = StatusMsg;
        Addpercentage(app, Per);
    case 4
        % Update status text and details message without percentage
        % app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Text = StatusMsg;
        app.StatusLabel.Tag = StatusMsg;
        MATDSSApp_Details(app, DetailsMsg);
    case 5
        % Update status text, percentage, and details message
        % app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Tag = StatusMsg;
        Addpercentage(app, Per);
        MATDSSApp_Details(app, DetailsMsg);
    otherwise
        % If MsgFlag is not recognized, do nothing
        % Option can be expanded or handled in future versions
end

% Brief pause to ensure UI updates are reflected
pause(0.001);

end

function Addpercentage(app, Per)
    % Addpercentage(app, Per)
    % Helper function to append percentage to the status label text.
    %
    % Parameters:
    %   - app: The application object containing the status label.
    %   - Per: The percentage value to be displayed.
    
    app.StatusLabel.Text = strcat(app.StatusLabel.Tag, " (" , num2str(round(Per)), "%)");
end



%% Old function Code
% 
%{

function MATDSSApp_Status(app, StatusMsg, DetailsMsg, Per)
% MATDSSApp_Status(app, StatusMsg, MsgFlag, DetailsMsg)
% This function will handle all status messages in MATDSS Application. You
% can updte status information, and add details message too if needed by
% calling this function.

MsgFlag = 5;

if nargin == 2
    if isnumeric(StatusMsg)
        Per = StatusMsg*100;
        MsgFlag = 2; %Case 2 --> update only the percentage number
        DetailsMsg = '';
    else
        MsgFlag = 1; %Case 1 --> change the whole status text
        DetailsMsg = '';
        Per = -1;
    end
end

if nargin == 3
    if ischar(DetailsMsg) || isstring(DetailsMsg)
        MsgFlag = 4; % Case 4 --> update status and details without percentage
        Per = -1;
    elseif isnumeric(DetailsMsg) % Case 3 --> update status with percentage
        Per = DetailsMsg*100;
        MsgFlag = 3;
    end
end


if ischar(StatusMsg)
    StatusMsg = convertCharsToStrings(StatusMsg);
end

if ischar(DetailsMsg)
    DetailsMsg = convertCharsToStrings(DetailsMsg);
end

switch MsgFlag
    case 1
        app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Text = StatusMsg;
        app.StatusLabel.Tag = StatusMsg;
    case 2
        app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        Addpercentage(app,Per);
    case 3
        app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Tag = StatusMsg;
        Addpercentage(app,Per);
    case 4
        app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Text = StatusMsg;
        app.StatusLabel.Tag = StatusMsg;
        MATDSSApp_Details(app, DetailsMsg);
    case 5
        app.StatusLabel.FontColor = app.MyRun.Plot.Colors(10,:);
        app.StatusLabel.Tag = StatusMsg;
        Addpercentage(app,Per);
        MATDSSApp_Details(app, DetailsMsg);
    otherwise
        % Do nothing :D
end

pause(0.001);
end


function Addpercentage(app,Per)
    app.StatusLabel.Text = strcat(app.StatusLabel.Tag, " (" ,num2str(round(Per)),"%)");
end

%}