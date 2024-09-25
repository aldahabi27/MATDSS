function MATDSSApp_P0SetFunButtonGroup(app)
    % MATDSSApp_P0SetFunButtonGroup(app)
    % This function modifies the GUI for the P0Set function according to the selected option.
    % To change the default values, modify them within this function code.
    %
    % Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
    %
    % MATDSS Application
    % Copyright (c) 2023, Ilyas Farhat
    % by Ilyas Farhat
    % This file is part of the MATDSS Application.
    % Contact the developer at ilyas.farhat@outlook.com

    % Read the selected option
    selectedButton = app.P0SetFunctionButtonGroup.SelectedObject;

    % Determine the selected P0Set function
    P0SetSelected = MATDSS_StrComp({'constant', 'unit step', 'ramp'}, lower(selectedButton.Text));
    
    % Update GUI elements based on the selected function
    switch P0SetSelected
        case 1 % Constant function (no change in P0Setpoint)
            app.ChangekWEditField.Enable = 'off';
            app.ChangekWEditFieldLabel.Enable = 'off';
            app.ChangekWEditFieldLabel.Text = 'Change (kW)';
            app.ChangekWEditFieldLabel.Interpreter = 'none';
            app.ChangekWEditField.Value = '0';
            
        case 2 % Unit Step function
            app.ChangekWEditField.Enable = 'on';
            app.ChangekWEditFieldLabel.Enable = 'on';
            app.ChangekWEditFieldLabel.Text = '$\Delta \mathrm{P}_{0} (kW)$';
            app.ChangekWEditFieldLabel.Interpreter = 'latex';
            app.ChangekWEditField.Value = '-200';
            
        case 3 % Ramp function
            app.ChangekWEditField.Enable = 'on';
            app.ChangekWEditFieldLabel.Enable = 'on';
            app.ChangekWEditFieldLabel.Text = '$\Delta \mathrm{P}_{0,\mathrm{step}} (kW)$';
            app.ChangekWEditFieldLabel.Interpreter = 'latex';
            app.ChangekWEditField.Value = '-20';
            
        otherwise % Error case (should not happen)
            app.ChangekWEditField.Enable = 'off';
            app.ChangekWEditFieldLabel.Enable = 'off';
            app.ChangekWEditFieldLabel.Text = 'Change (kW)';
            app.ChangekWEditFieldLabel.Interpreter = 'none';
            app.ChangekWEditField.Value = 'Error';
    end

end
