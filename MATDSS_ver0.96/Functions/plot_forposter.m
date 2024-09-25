% Ilyas Code for plotting three subplots stacked - same x-axis
%
%
%
%

close all
ncase = 0;
h = figure(1);%% Plotting Current

%%

ILegend = {};
Iax = subplot(3,1,3);
hold on;
ILim = 125; %Current limit

MATDSS = MyPlotsData{1}.MATDSS;

plot(Iax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat(ILim,size(MATDSS.Sim.Meas.at)), '--', 'Color',[MyColors(8,:), 0.], 'linewidth', 2);
ILegend = [ILegend;""];

for i = 1:length(MyPlotsData) - ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    Ibranches = {'l3.1','l3.2','l3.3'};
    Legendflag = true;
    for j = 1:length(Ibranches)
        IBranchName = Ibranches{j};
        IBranchIndex = find(strcmp(AllMeasBranchNames,IBranchName));
        if IBranchIndex > 0
            plot(Iax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.IProfile(IBranchIndex,:),'-','LineWidth',1,'Color', MyColors(i,:));
            if Legendflag
                ILegend = [ILegend; ['$\textrm{Config. }' num2str(i) '$']];
                Legendflag = false;
            else
                ILegend = [ILegend;""];
            end
        end
    end
end
grid on;
grid(Iax,'on')
box(Iax,'on')
xlabel(Iax,'$\rm{t} \ (s)$');
ylabel(Iax,'$\rm{L3} \ I \ (A)$');
set(Iax, 'FontSize', 12);
hold off;

disp('done')


%% Plotting Voltage

VLegend = {};
Vax = subplot(3,1,2);
hold on
VLim = 0.05; %Voltage limits

MATDSS = MyPlotsData{1}.MATDSS;


plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 - VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(8,:), 'linewidth', 2);
plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 + VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(8,:), 'linewidth', 2);
VLegend = [VLegend;""];
VLegend = [VLegend;""];

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    Vnodes = {'n5.1','n5.2','n5.3'};
    Legendflag = true;
    for j = 1:length(Vnodes)
        VNodeName = Vnodes{j};
        VNodeIndex = find(strcmp(AllNodeNames,VNodeName));
        if VNodeIndex > 0
            plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.VMagProfilePu(VNodeIndex,:),'-','LineWidth',1,'Color', MyColors(i,:));
            if Legendflag
                VLegend = [VLegend; ['$\textrm{Config. }' num2str(i) '$']];
                Legendflag = false;
            else
                VLegend = [VLegend; ""];
            end
        end
    end
end
grid on;
grid(Vax,'on')
box(Vax,'on')
ylabel(Vax,'$\rm{n}5 \ V \ (p.u.)$');

set(Vax, 'FontSize', 12);
xlim(Vax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
hold off;


%% Plotting Delta P0
P0ax = subplot(3,1,1);
hold on
P0Legend = {};
ContE = 5000; %5kW limits

MATDSS = MyPlotsData{1}.MATDSS;
DER = MyPlotsData{1}.DER;
AllDERNames = {DER.DSSName};
k_v = MATDSS.Meas.k_v;
AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

k_i = MATDSS.Meas.k_i;
AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);

% Plotting 1-time curves

% Delta P0 set-points
plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));
P0Legend = [P0Legend; '$\Delta P_{0,set-meas}$'];

% Delta P0 Goal Limits
plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3, '--', 'LineWidth', 2, 'Color',MyColors(8,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3, '--', 'LineWidth', 2, 'Color',MyColors(8,:));
P0Legend = [P0Legend;""];
P0Legend = [P0Legend;""];

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
    P0Legend = [P0Legend;['$\Delta P_{0,\ \textrm{Config. }' num2str(i) '}$']];
end
grid on;
box(P0ax,'on')
ylabel(P0ax,'$\rm{\Delta P}_{0} \ (kW)$');
set(P0ax, 'FontSize', 12);
xlim(P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])


legend(P0ax,'','','','1CA','$\rm{2CA}_{config. 1}$', '$\rm{2CA}_{config. 2}$','interpreter','latex','NumColumns',4, 'Position', [0.347214096628135 0.868641111770643 0.548428949600873 0.0506666682107106]);



%%

P0ax_position = get(P0ax,'Position');
xlabel(P0ax,"")
box(P0ax,"off")
xticklabels(P0ax, {})
% xticks(P0ax,{})
% xtickformat(P0ax,)
set(P0ax.XAxis, "TickDir", "both");


Vax_position = P0ax_position + [0, -P0ax_position(4), 0, 0];
set(Vax,'Position', Vax_position)
set(Vax.XAxis, "TickDir", "both");
xticklabels(Vax,{})
xlabel(Vax,"")
box(Vax,"off")


Iax_position = Vax_position + [0, -Vax_position(4), 0, 0];
set(Iax,'Position', Iax_position)
% xticklabels(Iax,{})
% xlabel(Vax,"")
box(Iax,"off")
% set(Vax, 'Position', [0.08 0.5 0.9 0.15]);
% set(P0ax, 'Position', [0.08 0.7 0.9 0.15]);


myxlim = [0 60];
myxlim = [0 10];
% myxlim = [0 5];
set(P0ax, "XLim", myxlim)
set(Iax, "XLim", myxlim)
set(Vax, "XLim", myxlim)

set(P0ax, "XTickMode", "auto")
set(Vax, "XTickMode", "auto")
set(Iax, "XTickMode", "auto")

P0axylim = [0 120];
P0axylim = [-650 0];
P0axylim = [-260 0];
Vaxylim = [0.946 0.965];
Iaxylim = [112 135];
set(P0ax, "YLim", P0axylim)
set(Vax, "YLim", Vaxylim)
set(Iax, "YLim", Iaxylim)

% legend(P0ax, P0Legend, 'interpreter', 'latex','Location','northeast','NumColumns',2);




MyPositionPlot = [0 0 12.80/1.5 7.20/1.5];

%%
set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)
print (h,'-dpng', '-r300', [CaseFolder '\' CaseNumber '_Combined2.png'])
% print(h, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
saveas(h, [CaseFolder '\' CaseNumber '_Combined2.fig'])
saveas(h, [CaseFolder '\' CaseNumber '_Combined2.eps'])
orient(h,'landscape');
saveas(h, [CaseFolder '\' CaseNumber '_Combined2.pdf'])
% stackedplot(copyobj(figure,P0ax))




%%

izoom = 1;
xzoom = [14 30];
set(P0ax, "XLim", xzoom)
set(Iax, "XLim", xzoom)
set(Vax, "XLim", xzoom)


p0axzoom = [-290 -80];
vaxzoom = [0.945 0.957];
iaxzoom = [108 136];
set(P0ax, "YLim", p0axzoom)
set(Vax, "YLim", vaxzoom)
set(Iax, "YLim", iaxzoom)


%%

% set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)
print (h,'-dpng', '-r300', [CaseFolder '\' CaseNumber '_zoomed' num2str(izoom) '.png'])
% print(h, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
saveas(h, [CaseFolder '\' CaseNumber '_zoomed' num2str(izoom) '.fig'])
saveas(h, [CaseFolder '\' CaseNumber '_zoomed' num2str(izoom) '.eps'], 'epsc')
orient(h,'landscape');
saveas(h, [CaseFolder '\' CaseNumber '_zoomed' num2str(izoom) '.pdf'])
% stackedplot(copyobj(figure,P0ax))


