function MATDSSApp_ProgressBar(Percent, ProgBar, BarText, ProgressBarStatus)
% MATDSSApp_ProgressBar(Percent,ProgBar,BarText) function is developed to
% control the progress bars in MATDSS Application.
%
% This function will update the progress bar indicator to the percentage
% provided in "Percent" (normalized to 1).
%
% ProgBar is the handle to the progress bar element (button in mlapp of
% MATLAB).
%
% BarText is an extra text that can be set (label of the button in mlapp)
% to show additional information. (Default is empty text)
%
% MATDSSApp_ProgressBar function serves as the progress bar controller 
% for the MATDSS application. It updates the visual representation 
% of the progress during simulations and tasks. The function handles 
% input parameters for progress percentage, a handle to the progress 
% bar UI element, an optional text label for additional information, 
% and an optional status display for current progress.
%
% Parameters:
%   - Percent: A normalized value (between 0 and 1) representing the 
%              progress completion percentage.
%   - ProgBar: The handle to the progress bar element in the MATLAB 
%              application.
%   - BarText: Optional string to display additional information on 
%              the progress bar.
%   - ProgressBarStatus: Optional handle for a status text area that 
%                        shows the progress percentage in text form.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

if nargin < 3
    BarText = ''; % Set BarText to empty if not provided
    ProgressBarStatus = []; % Initialize ProgressBarStatus to empty
end

if Percent <= 0
    % If Percent is less than or equal to 0, reset the progress bar color 
    % to its background color
    ProgBar.Icon = permute(repmat(ProgBar.BackgroundColor, [size(ProgBar.Icon, 1), 1, size(ProgBar.Icon, 2)]), [1, 3, 2]);
else
    % Update the progress bar when Percent is greater than 0
    ProgBar.Text = BarText; % Set the text label of the progress bar
    % Calculate the current progress based on the size of the icon
    CurrentProgress = min(round((size(ProgBar.Icon, 2) - 2) * (Percent)), size(ProgBar.Icon, 2) - 2);
    ProgBarColor = ProgBar.Icon; % Get the current icon color
    
    % Update the color of the progress section of the bar to royal blue
    ProgBarColor(2:end-1, 2:CurrentProgress + 1, 1) = 0.8500; % Red channel
    ProgBarColor(2:end-1, 2:CurrentProgress + 1, 2) = 0.3250; % Green channel
    ProgBarColor(2:end-1, 2:CurrentProgress + 1, 3) = 0.0980; % Blue channel
    ProgBar.Icon = ProgBarColor; % Update the progress bar icon with new color

    if ~isempty(ProgressBarStatus)
        % If ProgressBarStatus is provided, update its text to show 
        % the current percentage as a string
        ProgressBarStatus.Text = [num2str(round(Percent * 100, 0)) '%'];
    end
end

% clearvars % Clear all variables from the workspace

end
