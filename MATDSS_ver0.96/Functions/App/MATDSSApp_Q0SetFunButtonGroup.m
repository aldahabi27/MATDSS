function MATDSSApp_Q0SetFunButtonGroup(app)
    % MATDSSApp_Q0SetFunButtonGroup(app)
    % This function modifies the GUI for the Q0Set function according to the selected option.
    % To change the default values, modify them within this function code.
    %
    % Last Update for this function was on MATDSS App Ver 0.96 (30 Aug. 2024)
    %
    % MATDSS Application
    % Copyright (c) 2024, Ilyas Farhat
    % by Ilyas Farhat
    % This file is part of the MATDSS Application.
    % Contact the developer at ilyas.farhat@outlook.com

    % Read the selected option
    selectedButton = app.Q0SetFunctionButtonGroup.SelectedObject;

    % Determine the selected Q0Set function
    Q0SetSelected = MATDSS_StrComp({'constant', 'unit step', 'ramp'}, lower(selectedButton.Text));
    
    % Update GUI elements based on the selected function
    switch Q0SetSelected
        case 1 % Constant function (no change in Q0Setpoint)
            app.ChangekVareditField.Enable = 'off';
            app.ChangekVareditFieldLabel.Enable = 'off';
            app.ChangekVareditFieldLabel.Text = 'Change (kVar)';
            app.ChangekVareditFieldLabel.Interpreter = 'none';
            app.ChangekVareditField.Value = '0';
            
        case 2 % Unit Step function
            app.ChangekVareditField.Enable = 'on';
            app.ChangekVareditFieldLabel.Enable = 'on';
            app.ChangekVareditFieldLabel.Text = '$\Delta \mathrm{Q}_{0} (kVar)$';
            app.ChangekVareditFieldLabel.Interpreter = 'latex';
            app.ChangekVareditField.Value = '0';
            
        case 3 % Ramp function
            app.ChangekVareditField.Enable = 'on';
            app.ChangekVareditFieldLabel.Enable = 'on';
            app.ChangekVareditFieldLabel.Text = '$\Delta \mathrm{Q}_{0,\mathrm{step}} (kVar)$';
            app.ChangekVareditFieldLabel.Interpreter = 'latex';
            app.ChangekVareditField.Value = '0';
            
        otherwise % Error case (should not happen)
            app.ChangekVareditField.Enable = 'off';
            app.ChangekVareditFieldLabel.Enable = 'off';
            app.ChangekVareditFieldLabel.Text = 'Change (kVar)';
            app.ChangekVareditFieldLabel.Interpreter = 'none';
            app.ChangekVareditField.Value = 'Error';
    end

    % Clean up unnecessary variables (optional)
    % clearvars -except app
end
