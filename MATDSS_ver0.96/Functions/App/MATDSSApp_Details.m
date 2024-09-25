function MATDSSApp_Details(app, msg, ScrollFlag)
% MATDSSApp_Details(app, msg, ScrollFlag)
% This function appends a message to the end of the Details textarea in the 
% MATDSS Application's main window. The message can be multiline text. The 
% function will scroll the textarea to the bottom unless ScrollFlag is set to 
% false.
%
% Parameters:
%   - app: The application object containing the Details textarea.
%   - msg: The message text to be appended. Can be multiline text.
%   - ScrollFlag: Optional flag indicating whether to scroll to the bottom.
%                 Default is true (scroll to bottom).
%
% Last Update: MATDSS App Ver 0.96 (12 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Set default value for ScrollFlag if not provided
if nargin < 3
    ScrollFlag = true;
end

% Check if the application version has changed
if ~strcmp(app.AppVer, app.MATDSSApplicationUIFigure.Tag)
    % If the version has changed, reset the Details textarea with the new message
    app.DetailsTextArea.Value = msg;
else
    % If the version has not changed, append the new message to the existing content
    app.DetailsTextArea.Value = [app.DetailsTextArea.Value; msg];
end

% Scroll to the bottom of the textarea if ScrollFlag is true
if ScrollFlag
    scroll(app.DetailsTextArea, 'bottom');
end

end



%% Old Function Code
%{
function MATDSSApp_Details(app,msg,ScrollFlag)
% MATDSSApp_Details(app,msg)
% This function will append msg at the end of Details textarea in MATDSS
% App main window. The msg can be multiline text. The function will scroll
% the textarea to bottom, unless ScrollFlag is given false.

if nargin < 3
    ScrollFlag = true;
end


if ~strcmp(app.AppVer, app.MATDSSApplicationUIFigure.Tag)
    app.DetailsTextArea.Value = msg;
else
    app.DetailsTextArea.Value = [app.DetailsTextArea.Value; msg];
end
scroll(app.DetailsTextArea,'bottom');

%}