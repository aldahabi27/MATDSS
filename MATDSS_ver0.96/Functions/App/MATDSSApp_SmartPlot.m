function PlotTiles = MATDSSApp_SmartPlot(app, MATDSS, DER, PlotTiles, Options)
% PlotHandle = MATDSSApp_SmartPlot(app, MATDSS, DER, PlotHandle)
% This function will update the plots based on the plot properties selected
% in the GUI, and will show the corresponding output on the plots. You may
% configure each plot differently (different curves for each one), or hit
% "Show all" to show all curves on the plot. This would include Simulation
% based valeus (measurements not visible to L2C and are assumed to be
% unmeasured/measurement not available for L2C/LLC operators in real-life).
%
%
% Options = Reset or Add;


if nargin < 5
    Options = 'reset';
end
Options = lower(Options);

switch Options
    case 'reset'
        MATDSS_index = 1:length(app.MyRun.MATDSS);
        app.MyRun.Plot.PlotHandle = [];
    case 'add'
        MATDSS.index = length(app.MyRun.MATDSS);
    otherwise
        MATDSSApp_Details(app,'Error in plotting. The index of MATDSS run is missing!',1);
        return;
end

MeasCheckedboxes = app.MeasTree.CheckedNodes;
SimCheckedboxes = app.SimTree.CheckedNodes;

switch PlotTiles.Tag
    case 'P0Tiles'
        P0ax = nexttile(PlotTiles);%subplot(isubplot(iGainLoop))
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(P0ax,'on')
        hold(P0ax,'on')
        box(P0ax,'on')
        for i = 1:size(MeasCheckedboxes,1)
            MeasCurve = MeasCheckedboxes(i);
            if MeasCurve.Tag ~= "Tree" % MeasCurve is a treenode - corresponding to a curve
                switch MeasCurve.Text
                    case [char(916) 'P0'] %Measured Delta P0 (x)
                        % Measured P0 that our controller sees (represented by
                        % 'x' in the plot and colored gray)
                        plot(P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MATDSS.Plot.Colors(8,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0,meas}$'];
                    case [char(916) 'P0,set'] %Measured Delta P0,set points (x)
                        plot(P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(sim_meas_at_index) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MATDSS.Plot.Colors(8,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0,set-meas}$'];
                    case [char(916) 'P0,set - limits'] %Delta P0,set limits
                        % Delta P0 requested setpoint limits
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*MATDSS.L2C.E)./1e3,'LineWidth', 2, 'Color',MATDSS.Plot.Colors(5,:));
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + 10.*MATDSS.L2C.E)./1e3,'LineWidth', 2, 'Color',MATDSS.Plot.Colors(5,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{{0,set}_{ll}}$'];
                        PlotLegend = [PlotLegend;'$\Delta P_{{0,set}_{ul}}$'];
                    case [char(916) 'P0,set - L2C E-limits'] % L2C - E limits
                        % L2C Controller limits (colored dark green)
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - MATDSS.L2C.E)./1e3, 'Color',MATDSS.Plot.Colors(11,:));
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + MATDSS.L2C.E)./1e3, 'Color',MATDSS.Plot.Colors(11,:));
                        PlotLegend = [PlotLegend;""];
                        PlotLegend = [PlotLegend;""];
                    otherwise %it is DER P curve
                        DERDSSName = MeasCurve.Text(3:end);
                        DER_index = MATDSS_StrComp(DER_names,DERDSSName);
                        plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DER_index).W(1:end-1)./1e3);
                        PlotLegend = [PlotLegend;['$P_{' DER(DER_index).DSSName '}$']];
                end
            end
        end
        for i = 1:size(SimCheckedboxes,1)
            SimCurve = SimCheckedboxes(i);
            if SimCurve.Tag ~= "Tree"
                switch SimCurve.Text
                    case [char(916) 'P0'] % Delta P0
                        plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MATDSS.Plot.Colors(1,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0}$']
                    case [char(916) 'P0,set'] % Delta P0,set
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MATDSS.Plot.Colors(12,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0,Set}$']
                end
            end
        end
        legend(P0ax, PlotLegend, 'interpreter', 'latex','Location','southeast');
        xlabel(P0ax,'$t \ (s)$');
        ylabel(P0ax,'$\Delta P \ (kW)$');
        set(P0ax, 'FontSize', 14);
        title(P0ax,['Gain = ' num2str(MATDSS.L2C.Gain), ', $\alpha$ = ' num2str(MATDSS.L2C.alpha)],'Interpreter','latex')
        xlim(P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);


end
%{

PlotHandle = app.PlotHandle.P0Tiles;

%}
PlotLegend = {};
nDER = MATDSS.Sim.nDER;
DER_names = {};
for i = 1:nDER
    DER_names = [DER_names;DER(i).DSSName];
end
[~,sim_meas_at_index] = ismember(MATDSS.Meas.at,MATDSS.Sim.Meas.at);



switch PlotTiles.Tag
    case 'P0Tiles'
        P0ax = nexttile(PlotTiles);%subplot(isubplot(iGainLoop))
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        grid(P0ax,'on')
        hold(P0ax,'on')
        box(P0ax,'on')
        for i = 1:size(MeasCheckedboxes,1)
            MeasCurve = MeasCheckedboxes(i);
            if MeasCurve.Tag ~= "Tree" % MeasCurve is a treenode - corresponding to a curve
                switch MeasCurve.Text
                    case [char(916) 'P0'] %Measured Delta P0 (x)
                        % Measured P0 that our controller sees (represented by
                        % 'x' in the plot and colored gray)
                        plot(P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MATDSS.Plot.Colors(8,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0,meas}$'];
                    case [char(916) 'P0,set'] %Measured Delta P0,set points (x)
                        plot(P0ax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(sim_meas_at_index) - MATDSS.ControlSignals.P0Set(1))./1e3,'x','Color',MATDSS.Plot.Colors(8,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0,set-meas}$'];
                    case [char(916) 'P0,set - limits'] %Delta P0,set limits
                        % Delta P0 requested setpoint limits
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - 10.*MATDSS.L2C.E)./1e3,'LineWidth', 2, 'Color',MATDSS.Plot.Colors(5,:));
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + 10.*MATDSS.L2C.E)./1e3,'LineWidth', 2, 'Color',MATDSS.Plot.Colors(5,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{{0,set}_{ll}}$'];
                        PlotLegend = [PlotLegend;'$\Delta P_{{0,set}_{ul}}$'];
                    case [char(916) 'P0,set - L2C E-limits'] % L2C - E limits
                        % L2C Controller limits (colored dark green)
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - MATDSS.L2C.E)./1e3, 'Color',MATDSS.Plot.Colors(11,:));
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + MATDSS.L2C.E)./1e3, 'Color',MATDSS.Plot.Colors(11,:));
                        PlotLegend = [PlotLegend;""];
                        PlotLegend = [PlotLegend;""];
                    otherwise %it is DER P curve
                        DERDSSName = MeasCurve.Text(3:end);
                        DER_index = MATDSS_StrComp(DER_names,DERDSSName);
                        plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,DER(DER_index).W(1:end-1)./1e3);
                        PlotLegend = [PlotLegend;['$P_{' DER(DER_index).DSSName '}$']];
                end
            end
        end
        for i = 1:size(SimCheckedboxes,1)
            SimCurve = SimCheckedboxes(i);
            if SimCurve.Tag ~= "Tree"
                switch SimCurve.Text
                    case [char(916) 'P0'] % Delta P0
                        plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MATDSS.Plot.Colors(1,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0}$']
                    case [char(916) 'P0,set'] % Delta P0,set
                        plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MATDSS.Plot.Colors(12,:));
                        PlotLegend = [PlotLegend;'$\Delta P_{0,Set}$']
                end
            end
        end
        legend(P0ax, PlotLegend, 'interpreter', 'latex','Location','southeast');
        xlabel(P0ax,'$t \ (s)$');
        ylabel(P0ax,'$\Delta P \ (kW)$');
        set(P0ax, 'FontSize', 14);
        title(P0ax,['Gain = ' num2str(MATDSS.L2C.Gain), ', $\alpha$ = ' num2str(MATDSS.L2C.alpha)],'Interpreter','latex')
        xlim(P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
        pause(0.001);

app.MyRun.Plot.PlotHandle.P0Tiles = PlotTiles;
        
    case 'P_DERTiles'
    case 'Q0Tiles'
    case 'Q_DERTiles'
    case 'VoltageTiles'
    case 'CurrentTiles'
    otherwise
        MATDSSApp_Details(app,{'Error encounterd in MATDSSApp_SmartPlot function.';' '; ['MATDSS Application cannot find ' PlotTiles.Tag 'tiled layout handle!'];'******************************************************';' '});

end
