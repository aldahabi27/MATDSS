function MATDSSApp_Plot(app)
% MATDSSApp_Plot(app)
% This function updates the current plot in the MATDSS application based on
% the selected curves in the Plot Properties Panel (PPP) on the right.
%
% It identifies the currently visible plot and updates it with the selected
% measurement curves and simulation data. This allows users to visualize 
% different aspects of the simulation results interactively.
%
% Future enhancements may include the implementation of separate legends 
% for each plot in the Properties Panel, improving user experience and 
% data interpretation.
%
% Parameters:
%   - app: The MATDSS application instance, containing all relevant 
%          properties and methods for plot management.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com



% Retrieve necessary data from the application state
MATDSS = app.MyRun.MATDSS;  % Main data structure containing MATDSS info
DER = app.MyRun.DER;        % Distributed Energy Resources (DER) data

% Plot color settings
MyColors = app.MyRun.Plot.Colors; 

% Retrieve names of DERs and measurement nodes/branches
AllDERNames = {DER.DSSName};
k_v = MATDSS.Meas.k_v;  % Measurement index for nodes
AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

k_i = MATDSS.Meas.k_i;  % Measurement index for branches
AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);

% Switch case to determine which plot is currently selected in the PPP
switch app.PPPTabGroup.SelectedTab.Title
    case [char(916) 'P0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   DeltaP0 Plot                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Clear existing children from the P0Panel
        delete(app.P0Panel.Children)
        
        % Create a new axes object within the P0Panel
        P0ax = axes(app.P0Panel); % subplot(isubplot(iGainLoop))
        
        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        
        % Enable grid and hold for multiple plots on the same axes
        grid(P0ax,'on')
        hold(P0ax,'on')
        box(P0ax,'on')
        
        % Initialize legend cell array
        P0Legend = {};
        
        % Extract control element for plotting limits
        ContE = MATDSS.Cont.CA(1).E;
        
        % Plot measured P0 data if included in the list box
        for i = 1:length(app.P0MeasListBox.Value)
            MyNode = app.P0MeasListBox.Value{i};
            switch MyNode
                case [char(916) 'P0']
                    % Plot measured P0 data (gray 'x' markers)
                    plot(P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MyColors(8,:));
                    P0Legend = [P0Legend; '$\Delta P_{0,\rm{meas}}$'];
                    
                case [char(916) 'P0,set']
                    % Plot requested setpoints (dashed line)
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P0Legend = [P0Legend; '$\Delta P_{0,set-meas}$'];
                    
                case [char(916) 'P0,set - limits']
                    % Plot requested setpoints limits (application limits in green)
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    P0Legend = [P0Legend;'$\Delta P_{{0,set}_{ll}}$'];
                    P0Legend = [P0Legend;'$\Delta P_{{0,set}_{ul}}$'];
                    
                case [char(916) 'P0,set - E-limits']
                    % Plot controller limits (dark green)
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3, 'Color',MyColors(11,:));
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3, 'Color',MyColors(11,:));
                    P0Legend = [P0Legend;""]; % Empty legend entries
                    P0Legend = [P0Legend;""]; % Empty legend entries
                    
                otherwise
                    if strcmp(MyNode(1:2),'P_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            % Plot specific DER data if found
                            plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DERIndex).W(1:end-1,2)./1e3);
                            P0Legend = [P0Legend;['$P_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end
        
        % Plot simulation data if included in the list box
        for i = 1:length(app.P0SimListBox.Value)
            MyNode = app.P0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'P0'] % Delta P0
                    % Plot simulation Delta P0 (line width 2)
                    plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    P0Legend = [P0Legend;'$\Delta P_{0}$'];
                    
                case [char(916) 'P0,set'] % Delta P0,set
                    % Plot simulation Delta P0,set (dashed line)
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P0Legend = [P0Legend;'$\Delta P_{0,Set}$'];
            end
        end
        
        % Set axis labels and legend
        xlabel(P0ax,'$t \ (s)$');
        ylabel(P0ax,'$\Delta P \ (kW)$');
        legend(P0ax, P0Legend, 'interpreter', 'latex','Location','southeast');
        set(P0ax, 'FontSize', 14);
        title(P0ax,['Global G = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        
        % Set axis limits
        xlim(P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        
        % Pause briefly to ensure updates are visible
        pause(0.001);



    case ['CA-' char(916) 'P0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %              CA - DeltaP0 Plot                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear existing children from the CA_P0Panel
        delete(app.CA_P0Panel.Children)

        % Create a new axes object within the CA_P0Panel
        CA_P0ax = axes(app.CA_P0Panel); % subplot(isubplot(iGainLoop))

        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        % Enable grid and hold for multiple plots on the same axes
        grid(CA_P0ax,'on')
        hold(CA_P0ax,'on')
        box(CA_P0ax,'on')

        % Initialize legend cell array for CA plot
        CA_P0Legend = {};

        % Plot measured and requested P0 data from CA controllers if included in the list box
        for i = 1:length(app.CAP0MeasListBox.Value)
            MyNode = app.CAP0MeasListBox.Value{i};
            MyNodeText = MyNode;

            % Extract CA index from node text
            i_index = strfind(MyNodeText,'_');
            iCA = str2double(MyNode(3:i_index-1));

            switch MyNodeText(i_index+1:end)
                case [char(916) 'P0']
                    CA = MATDSS.Cont.CA(iCA);
                    % Plot measured P0 data (line width 2)
                    if length(CA.P0(1,:)) > length(MATDSS.Meas.at)
                        CA.P0 = CA.P0(:,2:end);
                    end
                    plot(CA_P0ax,MATDSS.Meas.at(1:end)-MATDSS.Time.Sim.ST,(sum(CA.P0(:,1:end),1) - sum(CA.P0(:,1),1))./1e3,'LineWidth',2);
                    CA_P0Legend = [CA_P0Legend; ['$CA' num2str(iCA) '-\Delta P_{0}$']];

                case [char(916) 'P0,set']
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.P0Set) > length(MATDSS.Meas.at)
                        CA.P0Set = CA.P0Set(2:end);
                    end
                    % Plot requested setpoints (dashed line)
                    plot(CA_P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1))./1e3,'--','LineWidth',2);
                    CA_P0Legend = [CA_P0Legend; ['$CA' num2str(iCA) '-\Delta P_{0,\rm{set}}$']];

                case [char(916) 'P0,set - E-limits']
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.P0Set) > length(MATDSS.Meas.at)
                        CA.P0Set = CA.P0Set(2:end);
                    end
                    % Plot controller limits (dark green)
                    plot(CA_P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1) - CA.E)./1e3, 'Color', MyColors(11,:));
                    plot(CA_P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1) + CA.E)./1e3, 'Color', MyColors(11,:));
                    CA_P0Legend = [CA_P0Legend;""]; % Empty legend entries
                    CA_P0Legend = [CA_P0Legend;""]; % Empty legend entries

                otherwise
                    if strcmp(MyNode(1:2),'P_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            % Plot specific DER data if found
                            plot(CA_P0ax,DER(DERIndex).W(1:end-1,1),DER(DERIndex).W(1:end-1,2)./1e3);
                            CA_P0Legend = [CA_P0Legend;['$P_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end

        % Plot simulation data if included in the list box
        for i = 1:length(app.CAP0SimListBox.Value)
            MyNode = app.CAP0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'P0'] % Delta P0
                    % Plot simulation Delta P0 (line width 2)
                    plot(CA_P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    CA_P0Legend = [CA_P0Legend;'$\Delta P_{0}$'];

                case [char(916) 'P0,set'] % Delta P0,set
                    % Plot simulation Delta P0,set (dashed line)
                    plot(CA_P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    CA_P0Legend = [CA_P0Legend;'$\Delta P_{0,\rm{set}}$'];
            end
        end

        % Set axis labels and legend
        xlabel(CA_P0ax,'$t \ (s)$');
        ylabel(CA_P0ax,'$\Delta P \ (kW)$');
        legend(CA_P0ax, CA_P0Legend, 'interpreter', 'latex','Location','southeast');
        set(CA_P0ax, 'FontSize', 14);
        title(CA_P0ax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')

        % Set axis limits
        xlim(CA_P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])

        % Pause briefly to ensure updates are visible
        pause(0.001);




    case ['P_DER']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    P_DER Plot                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear existing children from the P_DERPanel
        delete(app.P_DERPanel.Children)

        % Create a new axes object within the P_DERPanel
        P_DERax = axes(app.P_DERPanel);

        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        % Enable grid and hold for multiple plots on the same axes
        grid(P_DERax,'on')
        hold(P_DERax,'on')
        box(P_DERax,'on')

        % Initialize legend cell array for P_DER plot
        P_DERLegend = {};

        % Get controller error value
        ContE = MATDSS.Cont.CA(1).E;

        % Plot measured and requested P0 data from DER controllers if included in the list box
        for i = 1:length(app.P_DERMeasListBox.Value)
            MyNode = app.P_DERMeasListBox.Value{i};
            switch MyNode
                case [char(916) 'P0']
                    % Plot measured P0 data (x marker, gray color)
                    plot(P_DERax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MyColors(8,:));
                    P_DERLegend = [P_DERLegend; '$\Delta P_{0,\rm{meas}}$'];

                case [char(916) 'P0,set']
                    % Plot requested setpoints (dashed line)
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P_DERLegend = [P_DERLegend; '$\Delta P_{0,set-meas}$'];

                case [char(916) 'P0,set - limits']
                    % Plot setpoints limits (application limits, colors: red for setpoint, green for limits)
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    P_DERLegend = [P_DERLegend;'$\Delta P_{{0,set}_{ll}}$'];
                    P_DERLegend = [P_DERLegend;'$\Delta P_{{0,set}_{ul}}$'];

                case [char(916) 'P0,set - L2C E-limits']
                    % Plot L2C controller limits (dark green)
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3, 'Color',MyColors(11,:));
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3, 'Color',MyColors(11,:));
                    P_DERLegend = [P_DERLegend;""]; % Empty legend entries
                    P_DERLegend = [P_DERLegend;""]; % Empty legend entries

                otherwise
                    if strcmp(MyNode(1:2),'P_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            % Plot specific DER data if found
                            plot(P_DERax,DER(DERIndex).W(1:end-1,1),DER(DERIndex).W(1:end-1,2)./1e3);
                            P_DERLegend = [P_DERLegend;['$P_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end

        % Plot simulation data if included in the list box
        for i = 1:length(app.P_DERSimListBox.Value)
            MyNode = app.P_DERSimListBox.Value{i};
            switch MyNode
                case [char(916) 'P0'] % Delta P0
                    % Plot simulation Delta P0 (line width 2)
                    plot(P_DERax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    P_DERLegend = [P_DERLegend;'$\Delta P_{0}$'];

                case [char(916) 'P0,set'] % Delta P0,set
                    % Plot simulation Delta P0,set (dashed line)
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P_DERLegend = [P_DERLegend;'$\Delta P_{0,\rm{set}}$'];
            end
        end

        % Set axis labels and legend
        legend(P_DERax,P_DERLegend, 'interpreter', 'latex','Location','southeast');
        xlabel(P_DERax,'$t \ (s)$');
        ylabel(P_DERax,'$\Delta P_{\mathrm{DER}} \ (kW)$');
        set(P_DERax, 'FontSize', 14);
        title(P_DERax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')

        % Set axis limits
        xlim(P_DERax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])

        % Pause briefly to ensure updates are visible
        pause(0.001);


    case [char(916) 'Q0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                DeltaQ0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear existing children from the Q0Panel
        delete(app.Q0Panel.Children)

        % Create a new axes object within the Q0Panel
        Q0ax = axes(app.Q0Panel);

        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        % Enable grid and hold for multiple plots on the same axes
        grid(Q0ax,'on')
        hold(Q0ax,'on')
        box(Q0ax,'on')

        % Initialize legend cell array for Q0 plot
        Q0Legend = {};

        % Get controller error value
        ContE = MATDSS.Cont.CA(1).E;

        % Plot measured and requested Q0 data from DER controllers if included in the list box
        for i = 1:length(app.Q0MeasListBox.Value)
            MyNode = app.Q0MeasListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0']
                    % Plot measured Q0 data (x marker, gray color)
                    plot(Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.Q0,1) - MATDSS.ControlSignals.Q0Set(1))./1e3,'x','Color',MyColors(8,:));
                    Q0Legend = [Q0Legend; '$\Delta Q_{0,\rm{meas}}$'];

                case [char(916) 'Q0,set']
                    % Handle case for Q0 setpoints if needed (currently empty)

                case [char(916) 'Q0,set - limits']
                    % Handle case for Q0 setpoint limits if needed (currently empty)

                case [char(916) 'Q0,set - E-limits']
                    % Handle case for Q0 setpoint E-limits if needed (currently empty)

                otherwise
                    if strcmp(MyNode(1:2),'Q_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            % Plot specific DER Q0 data if found
                            plot(Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DERIndex).W(1:end-1,2)./1e3);
                            Q0Legend = [Q0Legend;['$Q_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end

        % Plot simulation data if included in the list box
        for i = 1:length(app.Q0SimListBox.Value)
            MyNode = app.Q0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0'] % Delta Q0
                    % Plot simulation Delta Q0 (line width 2)
                    plot(Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.Q0,1) - sum(MATDSS.Sim.Meas.Q0(:,1),1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    Q0Legend = [Q0Legend;'$\Delta Q_{0}$'];
            end
        end

        % Set axis labels and legend
        xlabel(Q0ax,'$t \ (s)$');
        ylabel(Q0ax,'$\Delta Q \ (kVar)$');
        legend(Q0ax, Q0Legend, 'interpreter', 'latex','Location','southeast');
        set(Q0ax, 'FontSize', 14);
        title(Q0ax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')

        % Set axis limits
        xlim(Q0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])

        % Pause briefly to ensure updates are visible
        pause(0.001);



    case ['CA-' char(916) 'Q0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %           CA - DeltaQ0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear existing children from the CA_Q0Panel
        delete(app.CA_Q0Panel.Children)

        % Create a new axes object within the CA_Q0Panel
        CA_Q0ax = axes(app.CA_Q0Panel); % Create axes for the plot

        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        % Enable grid and hold for multiple plots on the same axes
        grid(CA_Q0ax,'on')
        hold(CA_Q0ax,'on')
        box(CA_Q0ax,'on')

        % Initialize legend cell array for CA_Q0 plot
        CA_Q0Legend = {};

        % Plot measured and requested Q0 data from CA controllers if included in the list box
        for i = 1:length(app.CAQ0MeasListBox.Value)
            MyNode = app.CAQ0MeasListBox.Value{i};
            MyNodeText = MyNode;
            i_index = strfind(MyNodeText,'_');
            iCA = str2double(MyNode(3:i_index-1)); % Extract CA index
            switch MyNodeText(i_index+1:end)
                case [char(916) 'Q0']
                    CA = MATDSS.Cont.CA(iCA); % Get the controller CA
                    % Ensure CA.Q0 length matches the measurement length
                    if length(CA.Q0(1,:)) > length(MATDSS.Meas.at)
                        CA.Q0 = CA.Q0(:,2:end);
                    end
                    % Plot measured Q0 data (line width 2)
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(CA.Q0,1) - sum(CA.Q0(:,1),1))./1e3,'LineWidth',2);
                    CA_Q0Legend = [CA_Q0Legend; ['$CA' num2str(iCA) '-\Delta Q_{0}$']];

                case [char(916) 'Q0,set']
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.Q0Set) > length(MATDSS.Meas.at)
                        CA.Q0Set = CA.Q0Set(2:end);
                    end
                    % Plot Delta Q0 requested setpoints (dashed line)
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1))./1e3,'--','LineWidth',2);
                    CA_Q0Legend = [CA_Q0Legend; ['$CA' num2str(iCA) '-\Delta Q_{0,\rm{set}}$']];

                case [char(916) 'Q0,set - E-limits']
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.Q0Set) > length(MATDSS.Meas.at)
                        CA.Q0Set = CA.Q0Set(2:end);
                    end
                    % Plot Q0 setpoint limits (colored dark green)
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1) - CA.E)./1e3, 'Color', MyColors(11,:));
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1) + CA.E)./1e3, 'Color', MyColors(11,:));
                    CA_Q0Legend = [CA_Q0Legend;""]; % Adding empty entries for limits
                    CA_Q0Legend = [CA_Q0Legend;""];

                otherwise
                    if strcmp(MyNode(1:2),'Q_')
                        DERName = MyNode.Text(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            % Plot specific DER Q0 data if found
                            plot(CA_Q0ax,DER(DERIndex).Var(1:end-1,1),DER(DERIndex).Var(1:end-1,2)./1e3);
                            plot(CA_Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DERIndex).Var(1:end-1)./1e3);
                            CA_Q0Legend = [CA_Q0Legend;['$Q_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end

        % Plot simulation data if included in the list box
        for i = 1:length(app.CAQ0SimListBox.Value)
            MyNode = app.CAQ0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0'] % Delta Q0
                    % Plot simulation Delta Q0 (line width 2)
                    plot(CA_Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.Q0,1) - MATDSS.ControlSignals.Q0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    CA_Q0Legend = [CA_Q0Legend;'$\Delta Q_{0}$'];

                case [char(916) 'Q0,set'] % Delta Q0,set
                    % Plot simulation Delta Q0 requested setpoints (dashed line)
                    plot(CA_Q0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    CA_Q0Legend = [CA_Q0Legend;'$\Delta Q_{0,Set}$'];
            end
        end

        % Set axis labels and legend
        xlabel(CA_Q0ax,'$t \ (s)$');
        ylabel(CA_Q0ax,'$\Delta Q \ (kVar)$');
        legend(CA_Q0ax, CA_Q0Legend, 'interpreter', 'latex','Location','southeast');
        set(CA_Q0ax, 'FontSize', 14);
        title(CA_Q0ax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')

        % Set axis limits
        xlim(CA_Q0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])

        % Pause briefly to ensure updates are visible
        pause(0.001);




    case ['Q_DER']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    Q_DER Plot                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear existing children from the Q_DERPanel
        delete(app.Q_DERPanel.Children)

        % Create a new axes object within the Q_DERPanel
        Q_DERax = axes(app.Q_DERPanel);

        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        % Enable grid and hold for multiple plots on the same axes
        grid(Q_DERax,'on')
        hold(Q_DERax,'on')
        box(Q_DERax,'on')

        % Initialize legend cell array for Q_DER plot
        Q_DERLegend = {};
        ContE = MATDSS.Cont.CA(1).E; % Obtain the controller E value

        % Plot measured Q0 data and DER data if included in the list box
        for i = 1:length(app.Q_DERMeasListBox.Value)
            MyNode = app.Q_DERMeasListBox.Value{i};

            switch MyNode
                case [char(916) 'Q0']
                    % Plot measured Q0 data (line width 2, color gray)
                    plot(Q_DERax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.Q0,1)-sum(MATDSS.Meas.Q0(:,1),1))./1e3,'x','Color',MyColors(8,:));
                    Q_DERLegend = [Q_DERLegend; '$\Delta Q_{0,\rm{meas}}$'];

                otherwise
                    if strcmp(MyNode(1:2),'Q_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            % Plot specific DER Q0 data if found
                            plot(Q_DERax,DER(DERIndex).Var(1:end-1,1),DER(DERIndex).Var(1:end-1,2)./1e3);
                            Q_DERLegend = [Q_DERLegend;['$Q_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end

        % Plot simulation data if included in the list box
        for i = 1:length(app.Q_DERSimListBox.Value)
            MyNode = app.Q_DERSimListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0'] % Delta Q0
                    % Plot simulation Delta Q0 (line width 2)
                    plot(Q_DERax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.Q0,1) - sum(MATDSS.Sim.Meas.Q0(:,1),1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    Q_DERLegend = [Q_DERLegend;'$\Delta Q_{0}$'];
            end
        end

        % Set axis labels and legend
        xlabel(Q_DERax,'$t \ (s)$');
        ylabel(Q_DERax,'$\Delta Q_{\mathrm{DER}} \ (kVar)$');
        legend(Q_DERax,Q_DERLegend, 'interpreter', 'latex','Location','southeast');
        set(Q_DERax, 'FontSize', 14);
        title(Q_DERax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')

        % Set axis limits
        xlim(Q_DERax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])

        % Pause briefly to ensure updates are visible
        pause(0.001);



    case 'VP'

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  V_Profile Plot                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear existing children from the VoltagePanel
        delete(app.VoltagePanel.Children)

        % Create a new axes object within the VoltagePanel
        Vax = axes(app.VoltagePanel);

        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        % Enable grid and hold for multiple plots on the same axes
        grid(Vax,'on')
        hold(Vax,'on')
        box(Vax,'on')

        % Initialize legend cell array for voltage profile plot
        VLegend = {};

        % Plot measured voltage profiles if included in the list box
        for i = 1:length(app.VMeasListBox.Value)
            MyNode = app.VMeasListBox.Value{i};
            VNodeName = MyNode;
            VNodeIndex = MATDSS_StrComp(AllMeasNodeNames,VNodeName);
            if VNodeIndex > 0
                % Plot measured voltage magnitude profile
                plot(Vax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.VMagProfilePu(VNodeIndex,:));
                VLegend = [VLegend, VNodeName];
            end
        end

        % Plot simulated voltage profiles if included in the list box
        for i = 1:length(app.VSimListBox.Value)
            MyNode = app.VSimListBox.Value{i};
            VNodeName = MyNode;
            VNodeIndex = MATDSS_StrComp(AllNodeNames,VNodeName);
            if VNodeIndex > 0
                % Plot simulated voltage magnitude profile
                plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.VMagProfilePu(VNodeIndex,:));
                VLegend = [VLegend, VNodeName];
            end
        end

        % Set axis labels, legend, and other properties
        xlabel(Vax,'$t \ (s)$');
        ylabel(Vax,'$V \ (p.u.)$');
        legend(Vax,VLegend, 'interpreter', 'latex','Location','southeast');
        set(Vax, 'FontSize', 14);
        title(Vax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(Vax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])



    case 'IP'
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  I_Profile Plot                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Clear existing children from the CurrentPanel
        delete(app.CurrentPanel.Children)

        % Create a new axes object within the CurrentPanel
        Iax = axes(app.CurrentPanel);

        % Set default interpreters for LaTeX formatting
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');

        % Enable grid and hold for multiple plots on the same axes
        grid(Iax,'on')
        hold(Iax,'on')
        box(Iax,'on')

        % Initialize legend cell array for current profile plot
        ILegend = {};

        % Plot measured current profiles if included in the list box
        for i = 1:length(app.IMeasListBox.Value)
            MyNode = app.IMeasListBox.Value{i};
            IBranchName = MyNode;
            IBranchIndex = MATDSS_StrComp(AllMeasBranchNames,IBranchName);
            if IBranchIndex > 0
                % Plot measured current profile
                plot(Iax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.IProfile(IBranchIndex,:));
                ILegend = [ILegend, IBranchName];
            end
        end

        % Plot simulated current profiles if included in the list box
        for i = 1:length(app.ISimListBox.Value)
            MyNode = app.ISimListBox.Value{i};
            IBranchName = MyNode;
            IBranchIndex = MATDSS_StrComp(AllBranchNames,IBranchName);
            if IBranchIndex > 0
                % Plot simulated current profile
                plot(Iax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.IProfile(IBranchIndex,:));
                ILegend = [ILegend, IBranchName];
            end
        end

        % Set axis labels, legend, and other properties
        xlabel(Iax,'$t \ (s)$');
        ylabel(Iax,'$I \ (A)$');
        legend(Iax,ILegend, 'interpreter', 'latex','Location','southeast');
        set(Iax, 'FontSize', 14);
        title(Iax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(Iax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])

end


end %function


% Old Function Code (without comments)

%{ 

function MATDSSApp_Plot(app)
% MATDSSApp_Plot(app)
% This function will update the current plot following the selected curves
% in the PPP (Plot Properties Panel) on the right.
%
% The function will determine which 'plot' is currently visible, and update
% it accordingly. The code here could have some improvments in the future,
% and add features (like separate legend in the panel on the right).
%
%
% 
%
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



MyColors = app.MyRun.Plot.Colors;

AllDERNames = {DER.DSSName};
k_v = MATDSS.Meas.k_v;
AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

k_i = MATDSS.Meas.k_i;
AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);


switch app.PPPTabGroup.SelectedTab.Title
    case [char(916) 'P0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   DeltaP0 Plot                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(app.P0Panel.Children)
        P0ax = axes(app.P0Panel);%subplot(isubplot(iGainLoop))
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(P0ax,'on')
        hold(P0ax,'on')
        box(P0ax,'on')
        P0Legend = {};
        ContE = MATDSS.Cont.CA(1).E;
        for i = 1:length(app.P0MeasListBox.Value)
            MyNode = app.P0MeasListBox.Value{i};
            switch MyNode
                case [char(916) 'P0']
                    % Measured P0 that our controller sees (represented by
                    % 'x' in the plot and colored gray)
                    plot(P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MyColors(8,:));
                    P0Legend = [P0Legend; '$\Delta P_{0,\rm{meas}}$'];
                case [char(916) 'P0,set']
                    % Delta P0 requested setpoints
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P0Legend = [P0Legend; '$\Delta P_{0,set-meas}$'];
                case [char(916) 'P0,set - limits']
                    % Delta P0 requested setpoints limits (application
                    % limits) (colors-> setpoint is nice red, limits are
                    % Green)
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    P0Legend = [P0Legend;'$\Delta P_{{0,set}_{ll}}$'];
                    P0Legend = [P0Legend;'$\Delta P_{{0,set}_{ul}}$'];
                case [char(916) 'P0,set - E-limits']
                    % L2C Controller limits (colored dark green)
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3, 'Color',MyColors(11,:));
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3, 'Color',MyColors(11,:));
                    P0Legend = [P0Legend;""];
                    P0Legend = [P0Legend;""];
                otherwise
                    if strcmp(MyNode(1:2),'P_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DERIndex).W(1:end-1,2)./1e3);
                            P0Legend = [P0Legend;['$P_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end
        for i = 1:length(app.P0SimListBox.Value)
            MyNode = app.P0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'P0'] % Delta P0
                    plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    P0Legend = [P0Legend;'$\Delta P_{0}$'];
                case [char(916) 'P0,set'] % Delta P0,set
                    plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P0Legend = [P0Legend;'$\Delta P_{0,Set}$'];
            end

        end

        xlabel(P0ax,'$t \ (s)$');
        ylabel(P0ax,'$\Delta P \ (kW)$');
        legend(P0ax, P0Legend, 'interpreter', 'latex','Location','southeast');
        set(P0ax, 'FontSize', 14);
        title(P0ax,['Global G = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);

    case ['CA-' char(916) 'P0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %              CA - DeltaP0 Plot                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(app.CA_P0Panel.Children)

        CA_P0ax = axes(app.CA_P0Panel);%subplot(isubplot(iGainLoop))
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(CA_P0ax,'on')
        hold(CA_P0ax,'on')
        box(CA_P0ax,'on')
        CA_P0Legend = {};
        for i = 1:length(app.CAP0MeasListBox.Value)
            MyNode = app.CAP0MeasListBox.Value{i};
            MyNodeText = MyNode;
            i_index = strfind(MyNodeText,'_');
            iCA = str2double(MyNode(3:i_index-1));
            switch MyNodeText(i_index+1:end)
                case [char(916) 'P0']
                    CA = MATDSS.Cont.CA(iCA);
                    % Measured P0 that our controller sees (represented by
                    % 'x' in the plot and colored gray)
                    if length(CA.P0(1,:)) > length(MATDSS.Meas.at)
                        CA.P0 = CA.P0(:,2:end);
                    end
                    plot(CA_P0ax,MATDSS.Meas.at(1:end)-MATDSS.Time.Sim.ST,(sum(CA.P0(:,1:end),1) - sum(CA.P0(:,1),1))./1e3,'LineWidth',2);
                    CA_P0Legend = [CA_P0Legend; ['$CA' num2str(iCA) '-\Delta P_{0}$']];
                case [char(916) 'P0,set']
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.P0Set) > length(MATDSS.Meas.at)
                        CA.P0Set = CA.P0Set(2:end);
                    end
                    % Delta P0 requested setpoints
                    plot(CA_P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1))./1e3,'--','LineWidth',2);
                    CA_P0Legend = [CA_P0Legend; ['$CA' num2str(iCA) '-\Delta P_{0,\rm{set}}$']];
                case [char(916) 'P0,set - E-limits']
                    % L2C Controller limits (colored dark green)
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.P0Set) > length(MATDSS.Meas.at)
                        CA.P0Set = CA.P0Set(2:end);
                    end
                    plot(CA_P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1) - CA.E)./1e3, 'Color', MyColors(11,:));
                    plot(CA_P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.P0Set - CA.P0Set(1) + CA.E)./1e3, 'Color', MyColors(11,:));
                    CA_P0Legend = [CA_P0Legend;""];
                    CA_P0Legend = [CA_P0Legend;""];
                otherwise
                    if strcmp(MyNode(1:2),'P_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            plot(CA_P0ax,DER(DERIndex).W(1:end-1,1),DER(DERIndex).W(1:end-1,2)./1e3);
                            CA_P0Legend = [CA_P0Legend;['$P_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end
        for i = 1:length(app.CAP0SimListBox.Value)
            MyNode = app.CAP0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'P0'] % Delta P0
                    plot(CA_P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    CA_P0Legend = [CA_P0Legend;'$\Delta P_{0}$'];
                case [char(916) 'P0,set'] % Delta P0,set
                    plot(CA_P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    CA_P0Legend = [CA_P0Legend;'$\Delta P_{0,\rm{set}}$'];
            end
        end


        xlabel(CA_P0ax,'$t \ (s)$');
        ylabel(CA_P0ax,'$\Delta P \ (kW)$');
        legend(CA_P0ax, CA_P0Legend, 'interpreter', 'latex','Location','southeast');
        set(CA_P0ax, 'FontSize', 14);
        title(CA_P0ax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(CA_P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);


    case ['P_DER']

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    P_DER Plot                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(app.P_DERPanel.Children)
        P_DERax = axes(app.P_DERPanel);
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(P_DERax,'on')
        hold(P_DERax,'on')
        box(P_DERax,'on')
        P_DERLegend = {};
        ContE = MATDSS.Cont.CA(1).E;

        for i = 1:length(app.P_DERMeasListBox.Value)
            MyNode = app.P_DERMeasListBox.Value{i};
            switch MyNode
                case [char(916) 'P0']
                    % Measured P0 that our controller sees (represented by
                    % 'x' in the plot and colored gray)
                    plot(P_DERax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MyColors(8,:));
                    P_DERLegend = [P_DERLegend; '$\Delta P_{0,\rm{meas}}$'];
                case [char(916) 'P0,set']
                    % Delta P0 requested setpoints
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P_DERLegend = [P_DERLegend; '$\Delta P_{0,set-meas}$'];
                case [char(916) 'P0,set - limits']
                    % Delta P0 requested setpoints limits (application
                    % limits) (colors-> setpoint is nice red, limits are
                    % Green)
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + 10.*ContE)./1e3,'LineWidth', 2, 'Color',MyColors(5,:));
                    P_DERLegend = [P_DERLegend;'$\Delta P_{{0,set}_{ll}}$'];
                    P_DERLegend = [P_DERLegend;'$\Delta P_{{0,set}_{ul}}$'];
                case [char(916) 'P0,set - L2C E-limits']
                    % L2C Controller limits (colored dark green)
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3, 'Color',MyColors(11,:));
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3, 'Color',MyColors(11,:));
                    P_DERLegend = [P_DERLegend;""];
                    P_DERLegend = [P_DERLegend;""];
                otherwise
                    if strcmp(MyNode(1:2),'P_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            plot(P_DERax,DER(DERIndex).W(1:end-1,1),DER(DERIndex).W(1:end-1,2)./1e3);
                            P_DERLegend = [P_DERLegend;['$P_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end

        end

        for i = 1:length(app.P_DERSimListBox.Value)
            MyNode = app.P_DERSimListBox.Value{i};
            switch MyNode
                case [char(916) 'P0'] % Delta P0
                    plot(P_DERax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0,1) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    P_DERLegend = [P_DERLegend;'$\Delta P_{0}$'];
                case [char(916) 'P0,set'] % Delta P0,set
                    plot(P_DERax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    P_DERLegend = [P_DERLegend;'$\Delta P_{0,\rm{set}}$'];
            end

        end

        legend(P_DERax,P_DERLegend, 'interpreter', 'latex','Location','southeast');
        xlabel(P_DERax,'$t \ (s)$');
        ylabel(P_DERax,'$\Delta P_{\mathrm{DER}} \ (kW)$');
        set(P_DERax, 'FontSize', 14);
        title(P_DERax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(P_DERax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);




    case [char(916) 'Q0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                DeltaQ0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(app.Q0Panel.Children)
        Q0ax = axes(app.Q0Panel);
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(Q0ax,'on')
        hold(Q0ax,'on')
        box(Q0ax,'on')
        Q0Legend = {};
        ContE = MATDSS.Cont.CA(1).E;

        for i = 1:length(app.Q0MeasListBox.Value)

            MyNode = app.Q0MeasListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0']
                    % Measured Q0 that our controller sees (represented by
                    % 'x' in the plot and colored gray)
                    plot(Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.Q0,1) - MATDSS.ControlSignals.Q0Set(1))./1e3,'x','Color',MyColors(8,:));
                    Q0Legend = [Q0Legend; '$\Delta P_{0,\rm{meas}}$'];
                case [char(916) 'Q0,set']
                case [char(916) 'Q0,set - limits']
                case [char(916) 'Q0,set - E-limits']
                otherwise
                    if strcmp(MyNode(1:2),'Q_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            plot(Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DERIndex).W(1:end-1,2)./1e3);
                            Q0Legend = [Q0Legend;['$Q_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end
        for i = 1:length(app.Q0SimListBox.Value)
            MyNode = app.Q0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0'] % Delta Q0
                    plot(Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.Q0,1) - sum(MATDSS.Sim.Meas.Q0(:,1),1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    Q0Legend = [Q0Legend;'$\Delta Q_{0}$'];
            end
        end


        xlabel(Q0ax,'$t \ (s)$');
        ylabel(Q0ax,'$\Delta Q \ (kVar)$');
        legend(Q0ax, Q0Legend, 'interpreter', 'latex','Location','southeast');
        set(Q0ax, 'FontSize', 14);
        title(Q0ax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(Q0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);

    case ['CA-' char(916) 'Q0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %           CA - DeltaQ0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(app.CA_Q0Panel.Children)
        CA_Q0ax = axes(app.CA_Q0Panel);%subplot(isubplot(iGainLoop))
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(CA_Q0ax,'on')
        hold(CA_Q0ax,'on')
        box(CA_Q0ax,'on')
        CA_Q0Legend = {};

        for i = 1:length(app.CAQ0MeasListBox.Value)
            MyNode = app.CAQ0MeasListBox.Value{i};
            MyNodeText = MyNode;
            i_index = strfind(MyNodeText,'_');
            iCA = str2double(MyNode(3:i_index-1));
            switch MyNodeText(i_index+1:end)
                case [char(916) 'Q0']
                    CA = MATDSS.Cont.CA(iCA);
                    % Measured Q0 that our controller sees (represented by
                    % 'x' in the plot and colored gray)
                    if length(CA.Q0(1,:)) > length(MATDSS.Meas.at)
                        CA.Q0 = CA.Q0(:,2:end);
                    end
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(CA.Q0,1) - sum(CA.Q0(:,1),1))./1e3,'LineWidth',2);
                    CA_Q0Legend = [CA_Q0Legend; ['$CA' num2str(iCA) '-\Delta Q_{0}$']];
                case [char(916) 'Q0,set']
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.Q0Set) > length(MATDSS.Meas.at)
                        CA.Q0Set = CA.Q0Set(2:end);
                    end
                    % Delta Q0 requested setpoints
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1))./1e3,'--','LineWidth',2);
                    CA_Q0Legend = [CA_Q0Legend; ['$CA' num2str(iCA) '-\Delta Q_{0,\rm{set}}$']];
                case [char(916) 'Q0,set - E-limits']
                    CA = MATDSS.Cont.CA(iCA);
                    if length(CA.Q0Set) > length(MATDSS.Meas.at)
                        CA.Q0Set = CA.Q0Set(2:end);
                    end
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1) - CA.E)./1e3, 'Color', MyColors(11,:));
                    plot(CA_Q0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(CA.Q0Set - CA.Q0Set(1) + CA.E)./1e3, 'Color', MyColors(11,:));
                    CA_Q0Legend = [CA_Q0Legend;""];
                    CA_Q0Legend = [CA_Q0Legend;""];
                otherwise
                    if strcmp(MyNode(1:2),'Q_')
                        DERName = MyNode.Text(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            plot(CA_Q0ax,DER(DERIndex).Var(1:end-1,1),DER(DERIndex).Var(1:end-1,2)./1e3);
                            plot(CA_Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DERIndex).Var(1:end-1)./1e3);
                            CA_Q0Legend = [CA_Q0Legend;['$Q_{' DER(DERIndex).DSSName '}$']];
                        end
                    end
            end
        end
        for i = 1:length(app.CAQ0SimListBox.Value)
            MyNode = app.CAQ0SimListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0'] % Delta Q0
                    plot(CA_Q0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.Q0,1) - MATDSS.ControlSignals.Q0Set(1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    CA_Q0Legend = [CA_Q0Legend;'$\Delta Q_{0}$'];
                case [char(916) 'Q0,set'] % Delta Q0,set
                    plot(CA_Q0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.Q0Set - MATDSS.ControlSignals.Q0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(12,:));
                    CA_Q0Legend = [CA_Q0Legend;'$\Delta Q_{0,Set}$'];
            end
        end

        xlabel(CA_Q0ax,'$t \ (s)$');
        ylabel(CA_Q0ax,'$\Delta Q \ (kVar)$');
        legend(CA_Q0ax, CA_Q0Legend, 'interpreter', 'latex','Location','southeast');
        set(CA_Q0ax, 'FontSize', 14);
        title(CA_Q0ax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(CA_Q0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);




    case ['Q_DER']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    Q_DER Plot                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        delete(app.Q_DERPanel.Children)
        Q_DERax = axes(app.Q_DERPanel);
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(Q_DERax,'on')
        hold(Q_DERax,'on')
        box(Q_DERax,'on')
        Q_DERLegend = {};
        ContE = MATDSS.Cont.CA(1).E;

        for i = 1:length(app.Q_DERMeasListBox.Value)
            MyNode = app.Q_DERMeasListBox.Value{i};

            switch MyNode
                case [char(916) 'Q0']
                    % Measured Q0 that our controller sees (represented by
                    % 'x' in the plot and colored gray)
                    plot(Q_DERax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.Q0,1)-sum(MATDSS.Meas.Q0(:,1),1))./1e3,'x','Color',MyColors(8,:));
                    Q_DERLegend = [Q_DERLegend; '$\Delta Q_{0,\rm{meas}}$'];                   
                otherwise
                    if strcmp(MyNode(1:2),'Q_')
                        DERName = MyNode(3:end);
                        DERIndex = MATDSS_StrComp(AllDERNames,DERName);
                        if DERIndex > 0
                            plot(Q_DERax,DER(DERIndex).Var(1:end-1,1),DER(DERIndex).Var(1:end-1,2)./1e3);
                            Q_DERLegend = [Q_DERLegend;['$Q_{' DER(DERIndex).DSSName '}$']];
                        end
                    end

            end
        end
        for i = 1:length(app.Q_DERSimListBox.Value)
            MyNode = app.Q_DERSimListBox.Value{i};
            switch MyNode
                case [char(916) 'Q0'] % Delta Q0
                    plot(Q_DERax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.Q0,1) - sum(MATDSS.Sim.Meas.Q0(:,1),1))./1e3,'LineWidth',2,'Color', MyColors(1,:));
                    Q_DERLegend = [Q_DERLegend;'$\Delta Q_{0}$'];
            end
        end
        legend(Q_DERax,Q_DERLegend, 'interpreter', 'latex','Location','southeast');
        xlabel(Q_DERax,'$t \ (s)$');
        ylabel(Q_DERax,'$\Delta Q_{\mathrm{DER}} \ (kVar)$');
        set(Q_DERax, 'FontSize', 14);
        title(Q_DERax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(Q_DERax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);


    case 'VP'

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  V_Profile Plot                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        delete(app.VoltagePanel.Children)
        Vax = axes(app.VoltagePanel);
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(Vax,'on')
        hold(Vax,'on')
        box(Vax,'on')
        VLegend = {};

        for i = 1:length(app.VMeasListBox.Value)
            MyNode = app.VMeasListBox.Value{i};
            VNodeName = MyNode;
            VNodeIndex = MATDSS_StrComp(AllMeasNodeNames,VNodeName);
            if VNodeIndex > 0
                plot(Vax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.VMagProfilePu(VNodeIndex,:));
                VLegend = [VLegend, VNodeName];
            end
        end
        for i = 1:length(app.VSimListBox.Value)
            MyNode = app.VSimListBox.Value{i};
            VNodeName = MyNode;
            VNodeIndex = MATDSS_StrComp(AllNodeNames,VNodeName);
            if VNodeIndex > 0
                plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.VMagProfilePu(VNodeIndex,:));
                VLegend = [VLegend, VNodeName];
            end

        end
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(Vax,'on')
        box(Vax,'on')
        xlabel(Vax,'$t \ (s)$');
        ylabel(Vax,'$V \ (p.u.)$');
        legend(Vax,VLegend, 'interpreter', 'latex','Location','southeast');
        set(Vax, 'FontSize', 14);
        title(Vax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(Vax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])


    case 'IP'
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  I_Profile Plot                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(app.CurrentPanel.Children)
        Iax = axes(app.CurrentPanel);
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(Iax,'on')
        hold(Iax,'on')
        box(Iax,'on')
        ILegend = {};

        for i = 1:length(app.IMeasListBox.Value)
            MyNode = app.IMeasListBox.Value{i};
            IBranchName = MyNode;
            IBranchIndex = MATDSS_StrComp(AllMeasBranchNames,IBranchName);
            if IBranchIndex > 0
                plot(Iax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.IProfile(IBranchIndex,:));
                ILegend = [ILegend, IBranchName];
            end

        end


        for i = 1:length(app.ISimListBox.Value)
            MyNode = app.ISimListBox.Value{i};
            IBranchName = MyNode;
            IBranchIndex = MATDSS_StrComp(AllBranchNames,IBranchName);
            if IBranchIndex > 0
                plot(Iax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.IProfile(IBranchIndex,:));
                ILegend = [ILegend, IBranchName];
            end
        end


        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(Iax,'on')
        box(Iax,'on')
        xlabel(Iax,'$t \ (s)$');
        ylabel(Iax,'$I \ (A)$');
        legend(Iax,ILegend, 'interpreter', 'latex','Location','southeast');
        set(Iax, 'FontSize', 14);
        title(Iax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
        xlim(Iax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])




end
end %function


%}