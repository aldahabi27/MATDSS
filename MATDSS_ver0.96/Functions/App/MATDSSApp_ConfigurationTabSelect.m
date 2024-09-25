function MATDSSApp_ConfigurationTabSelect(app)
% MATDSSAppConfiguration_TabSelect(app)
% This function is triggered when a tab is selected in the Circuit 
% Configurations TabGroup. It enables or disables buttons based on the 
% selected tab and performs actions specific to each tab.
%
% Parameters:
%   - app: The application object containing the TabGroup and buttons.
%
% Last Update: MATDSS App Ver 0.96 (12 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Retrieve the currently selected tab from the TabGroup
selectedTab = app.CircuitConfigurationsTabGroup.SelectedTab;

% Determine which tab is selected and enable/disable buttons accordingly
switch selectedTab.Title
    case 'DERs'
        % If the 'DERs' tab is selected
        app.NewDERButton.Enable = true; % Enable the 'New DER' button
        
        % Enable the 'Delete DER' button if there is a selection in the DERs table
        if ~isempty(app.DERsTable.Selection)
            app.DeleteDERButton.Enable = true;
        else
            app.DeleteDERButton.Enable = false; % Disable 'Delete DER' button if no selection
        end
    case 'V & I'
        % If the 'V & I' tab is selected
        app.NewDERButton.Enable = false; % Disable the 'New DER' button
        app.DeleteDERButton.Enable = false; % Disable the 'Delete DER' button
end

end


%% Old Function Code

%{
function MATDSSAppConfiguration_TabSelect(app)
% MATDSSAppConfiguration_TabSelect(app)
% This function is triggered when a tab is selected in the Circuit
% Configurations TabGroup. It enables or disables buttons based on the
% selected tab and performs actions specific to each tab.

%   Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
%   MATDSS Application
%   Copyright (c) 2024, Ilyas Farhat
%   by Ilyas Farhat
%
%   This file is part of MATDSS Application
%   Contact the developer at ilyas.farhat@outlook.com

selectedTab = app.CircuitConfigurationsTabGroup.SelectedTab;
switch selectedTab.Title
    case 'DERs'
        app.NewDERButton.Enable = true;
        if ~isempty(app.DERsTable.Selection)
            app.DeleteDERButton.Enable = true;
        end
    case 'V & I'
        app.NewDERButton.Enable = false;
        app.DeleteDERButton.Enable = false;
end
end

%}