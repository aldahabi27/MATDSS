% function MATDSSApp_GUIChanges(app,Action)
% % This function will enable/disable (and possibly do more stuff later)
% % related to managing GUI
% if false
%     if MATDSS_StrComp(fieldnames(app),'MATDSSApplicationUIFigure') > 0
%         switch Action
%             case 'Disable_GUI'
%                 % app.DefaultButton.Enable = 'off';
%                 % app.BrowseButton.Enable = 'off';
%                 % app.OpenDSSFilesListBox.Enable = 'off';
%                 % app.SimulateButton.Enable = 'off';
%                 % app.P0SetFunctionButtonGroup.Enable = 'off';
%                 % app.TimeSettingsPanel.Enable = 'off';
%                 % app.SaveButton.Enable = 'off';
%                 % app.ConfigurationButton.Enable = 'off';
%             otherwise
%                 % app.DefaultButton.Enable = 'on';
%                 app.BrowseButton.Enable = 'on';
%                 app.OpenDSSFilesListBox.Enable = 'on';
%                 app.SimulateButton.Enable = 'on';
%                 app.P0SetFunctionButtonGroup.Enable = 'on';
%                 app.TimeSettingsPanel.Enable = 'on';
%                 app.SaveButton.Enable = 'on';
%                 % app.ConfigurationButton.Enable = 'on';
%         end
%     end
% end
% 
% 
% if MATDSS_StrComp(fieldnames(app),'MATDSSApplicationUIFigure') > 0
%     % app.DefaultButton.Enable = 'on';
%     app.BrowseButton.Enable = 'on';
%     app.OpenDSSFilesListBox.Enable = 'on';
%     app.SimulateButton.Enable = 'on';
%     app.P0SetFunctionButtonGroup.Enable = 'on';
%     app.TimeSettingsPanel.Enable = 'on';
%     app.SaveButton.Enable = 'on';
%     % app.ConfigurationButton.Enable = 'on';
% end
% 
% clearvars -except app