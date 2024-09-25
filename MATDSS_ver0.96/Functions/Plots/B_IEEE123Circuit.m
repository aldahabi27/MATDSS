clear;
clc;
close all;

currentFilePath = mfilename('fullpath');
[currentDirectory, ~, ~] = fileparts(currentFilePath);

i_slash = strfind(currentDirectory,'\');
% Plotting Step test results
%%
close all
ncase = 0;
h = figure(1);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (h, 'Units', 'normalized', 'Position', [0.25,0.25,0.4,0.3],'Visible', 'on');

CaseName = 'Step';
CaseFolder = [currentDirectory(1:i_slash(end-1)) 'exports\IEEE 123\' CaseName];     % initial subdirectory
% CaseFolder = uigetdir("C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\05 May\MATDSS_ver0.84\Exports");     % initial subdirectory


% Load all mat files
[~,CaseFiles] = dos(['dir /s /b ' '"' fullfile(CaseFolder,'*.mat') '"' ]);  % OS dir command
CaseFiles=textscan(CaseFiles,'%s','delimiter','\n'); CaseFiles=CaseFiles{:};

CaseNumber = strfind(CaseFolder,"\");
CaseNumber = CaseFolder(CaseNumber(end)+1);

CaseNumber(strfind(CaseNumber," ")) = "_";
MyPlotsData =[];
for i = 1:length(CaseFiles)
    MyPlotsData = [MyPlotsData;{load(CaseFiles{i})}];
end



% General parameters to control plots
MyColors = MyPlotsData{1}.Plot.Colors;
MyColors = [MyColors; '#C9C9C9';'#B9B9B9'];
MyPositionPlot = [0 0 12.80 7.20];
ContE = [5000, 10000]; %5kW & 10kW limits
MATDSS = MyPlotsData{1}.MATDSS;
DER = MyPlotsData{1}.DER;
AllDERNames = {DER.DSSName};
k_v = MATDSS.Meas.k_v;
AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

k_i = MATDSS.Meas.k_i;
AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);



% Plotting Step response
%{
%%% I_Plot
ILegend = {};
Iax = subplot(3,1,3);
hold on;
ILim = 125; %Current limit

MATDSS = MyPlotsData{1}.MATDSS;

plot(Iax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat(ILim,size(MATDSS.Sim.Meas.at)), '--', 'Color',[MyColors(end,:)], 'linewidth', 2);
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

% disp('done')


%%% Plotting Voltage

VLegend = {};
Vax = subplot(3,1,2);
hold on
VLim = 0.05; %Voltage limits

MATDSS = MyPlotsData{1}.MATDSS;


plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 - VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);
plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 + VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);
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


%%% Plotting Delta P0
P0ax = subplot(3,1,1);
hold on
P0Legend = {};
% ContE = 5000; %5kW limits

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
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/12)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/12) - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/12)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/12) - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/12:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/12:end) - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/12:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/12:end) - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
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


    
leg.ItemTokenSize = [10,15];

% Stack the plots

for i = 1:2

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

pause(1)
end


%}
P0ax = axes(h);
hold on
P0Legend = {};
% ContE = 5000; %5kW limits

MATDSS = MyPlotsData{1}.MATDSS;
DER = MyPlotsData{1}.DER;
AllDERNames = {DER.DSSName};
k_v = MATDSS.Meas.k_v;
AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

k_i = MATDSS.Meas.k_i;
AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);

% Delta P0 set-points
plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));
P0Legend = [P0Legend; '$\Delta P_{0,set-meas}$'];

% Delta P0 Goal Limits
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/6)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/6) - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/6)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/6) - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/6:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/6:end) - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/6:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/6:end) - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
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


leg = legend(P0ax,'','','','','','1CA','$\rm{6CA}_{RoU = 1}$', '$\rm{6CA}_{RoU = 5}$','interpreter','latex','NumColumns',1, 'Position', [0.794 0.805 0.1 0.05]);
    
leg.ItemTokenSize = [10,15];
%% Customize the plot

% General Results
myxlim = [0 25];
set(P0ax, "XLim", myxlim)
% set(Iax, "XLim", myxlim)
% set(Vax, "XLim", myxlim)

set(P0ax, "XTickMode", "auto")
% set(Vax, "XTickMode", "auto")
% set(Iax, "XTickMode", "auto")

P0axylim = [-275 0];
% Vaxylim = [0.946 0.965];
% Iaxylim = [112 135];
set(P0ax, "YLim", P0axylim)
% set(Vax, "YLim", Vaxylim)
% set(Iax, "YLim", Iaxylim)


% VaxYTicks = [0.95, 0.96];
% set(Vax, 'YTick', VaxYTicks)%, 'YTickLabelRotation', 80)

% IaxYTicks = [115:15:135];
% set(Iax, 'YTick', IaxYTicks)




%% MyPositionPlot = [0 0 12.80/1.5 7.20/1.5];
% set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)

figname = 'IEEE123StepTrackingWithDisturbance';

% Use this to export eps and force vector format
exportgraphics(h, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (h,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(h, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(h,'landscape');
set(h, 'Renderer', 'painters');
exportgraphics(h, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')


% set(h, 'color', 'none');
% print (h,'-depsc', '-r300', [CaseFolder '\' figname 'testingtest.eps'], '-tight')


% stackedplot(copyobj(figure,P0ax))




print(h, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');




%{
%% Prepare to Show current control on the right
set (h, 'Units', 'normalized', 'Position', [0.25,0.25,0.05,0.4],'Visible', 'on');

set(P0ax, "YLabel", [])
set(Iax, "YLabel", [])
set(Vax, "YLabel", [])

myxlim = [10 30];
set(P0ax, "XLim", myxlim)
set(Iax, "XLim", myxlim)
set(Vax, "XLim", myxlim)

set(P0ax, "XTickMode", "auto")
set(Vax, "XTickMode", "auto")
set(Iax, "XTickMode", "auto")

P0axylim = [-215 -185];
Vaxylim = [0.946 0.965];
Iaxylim = [123 126];
set(P0ax, "YLim", P0axylim)
set(Vax, "YLim", Vaxylim)
set(Iax, "YLim", Iaxylim)


VaxYTicks = [0.95, 0.96];
set(Vax, 'YTick', VaxYTicks)%, 'YTickLabelRotation', 80)

% IaxYTicks = [115:15:135];
set(Iax, 'YTickMode', 'auto')

legend(P0ax,'off')

% MyPositionPlot = [0 0 12.80/1.5 7.20/1.5];

% set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)



for i = 1:2

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

pause(1)
end

%
figname = '5BusStepTrackingWithDisturbance_CurrentControl';

% Use this to export eps and force vector format
exportgraphics(h, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (h,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(h, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(h,'landscape');
exportgraphics(h, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(h, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');

% stackedplot(copyobj(figure,P0ax))


% Generate inset figures (zoomed-in inset figures to show the curves)

% Current zoomed in figure
close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
Izax = axes(g);
hold on;
ILim = 125; %Current limit


plot(Izax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat(ILim,size(MATDSS.Sim.Meas.at)), '--', 'Color',[MyColors(end,:)], 'linewidth', 2);

for i = 1:length(MyPlotsData) - ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    Ibranches = {'l3.1','l3.2','l3.3'};
    for j = 1:length(Ibranches)
        IBranchName = Ibranches{j};
        IBranchIndex = find(strcmp(AllMeasBranchNames,IBranchName));
        if IBranchIndex > 0
            plot(Izax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.IProfile(IBranchIndex,:),'-','LineWidth',1,'Color', MyColors(i,:));
        end
    end
end
grid(Izax,'on')
box(Izax,'on')
set(Izax, 'FontSize', 18);
hold off;


set(Izax, "XLim", [0.7 2.5], "YLim", [114 118])


figname = '5BusStepTrackingWithDisturbance_Current1';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');

% stackedplot(copyobj(figure,P0ax))



% Voltage zoomed in figure

close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
Vzax = axes(g);
hold on
VLim = 0.05; %Voltage limits

plot(Vzax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 - VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);
plot(Vzax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 + VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    Vnodes = {'n5.1','n5.2','n5.3'};
    Legendflag = true;
    for j = 1:length(Vnodes)
        VNodeName = Vnodes{j};
        VNodeIndex = find(strcmp(AllNodeNames,VNodeName));
        if VNodeIndex > 0
            plot(Vzax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.VMagProfilePu(VNodeIndex,:),'-','LineWidth',1,'Color', MyColors(i,:));
        end
    end
end
grid(Vzax,'on')
box(Vzax,'on')

set(Vzax, 'FontSize', 18);
set(Vzax, "XLim", [0.6 2.5], "YLim", [0.95 0.9535])
hold off;


figname = '5BusStepTrackingWithDisturbance_Voltage1';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');

% stackedplot(copyobj(figure,P0ax))

%


%}
% Power zoomed in figure1
close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
P0zax = axes(g);
hold on;



% Delta P0 set-points
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));

% Delta P0 Goal Limits
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    plot(P0zax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
end
grid on;
box(P0zax,'on')
set(P0zax, 'FontSize', 18);
set(P0zax, "XLim", [0.6 2.5], "YLim", [-240 -160])
hold off;

figname = 'IEEE123StepTrackingWithDisturbance_Power1';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');


% Power zoomed in figure1
close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
P0zax = axes(g);
hold on;



% Delta P0 set-points
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));

% Delta P0 Goal Limits
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    plot(P0zax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
end
grid on;
box(P0zax,'on')
set(P0zax, 'FontSize', 18);
set(P0zax, "XLim", [10.6 12.5], "YLim", [-230 -155])
hold off;

figname = 'IEEE123StepTrackingWithDisturbance_Power2';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');




%% Plotting Ramp response
close all
clear h g Vax Iax P0ax Vzax Izax P0zax
close all
ncase = 0;
h = figure(1);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (h, 'Units', 'normalized', 'Position', [0.25,0.25,0.4,0.3],'Visible', 'on');

CaseName = 'Ramp';
CaseFolder = [currentDirectory(1:i_slash(end-1)) 'exports\IEEE 123\' CaseName];     % initial subdirectory
% CaseFolder = uigetdir("C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\05 May\MATDSS_ver0.84\Exports");     % initial subdirectory


% Load all mat files
[~,CaseFiles] = dos(['dir /s /b ' '"' fullfile(CaseFolder,'*.mat') '"' ]);  % OS dir command
CaseFiles=textscan(CaseFiles,'%s','delimiter','\n'); CaseFiles=CaseFiles{:};

CaseNumber = strfind(CaseFolder,"\");
CaseNumber = CaseFolder(CaseNumber(end)+1);

CaseNumber(strfind(CaseNumber," ")) = "_";
MyPlotsData =[];
for i = 1:length(CaseFiles)
    MyPlotsData = [MyPlotsData;{load(CaseFiles{i})}];
end



% General parameters to control plots
MyColors = MyPlotsData{1}.Plot.Colors;
MyColors = [MyColors; '#C9C9C9';'#B9B9B9'];
MyPositionPlot = [0 0 12.80 7.20];
ContE = [5000, 10000]; %5kW & 10kW limits
MATDSS = MyPlotsData{1}.MATDSS;
DER = MyPlotsData{1}.DER;
AllDERNames = {DER.DSSName};
k_v = MATDSS.Meas.k_v;
AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

k_i = MATDSS.Meas.k_i;
AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);



% Plotting Step response
%{
%%% I_Plot
ILegend = {};
Iax = subplot(3,1,3);
hold on;
ILim = 125; %Current limit

MATDSS = MyPlotsData{1}.MATDSS;

plot(Iax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat(ILim,size(MATDSS.Sim.Meas.at)), '--', 'Color',[MyColors(end,:)], 'linewidth', 2);
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

% disp('done')


%%% Plotting Voltage

VLegend = {};
Vax = subplot(3,1,2);
hold on
VLim = 0.05; %Voltage limits

MATDSS = MyPlotsData{1}.MATDSS;


plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 - VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);
plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 + VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);
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


%%% Plotting Delta P0
P0ax = subplot(3,1,1);
hold on
P0Legend = {};
% ContE = 5000; %5kW limits

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
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/12)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/12) - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/12)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/12) - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/12:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/12:end) - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/12:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/12:end) - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
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


    
leg.ItemTokenSize = [10,15];

% Stack the plots

for i = 1:2

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

pause(1)
end


%}
P0ax = axes(h);
hold on
P0Legend = {};
% ContE = 5000; %5kW limits

MATDSS = MyPlotsData{1}.MATDSS;
DER = MyPlotsData{1}.DER;
AllDERNames = {DER.DSSName};
k_v = MATDSS.Meas.k_v;
AllMeasNodeNames = MATDSS.Sim.Meas.AllNodesNames(k_v);
AllNodeNames = MATDSS.Sim.Meas.AllNodesNames;

k_i = MATDSS.Meas.k_i;
AllMeasBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames(k_i));
AllBranchNames = cellstr(MATDSS.Sim.Meas.AllBranchesPhasesNames);

% Delta P0 set-points
plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));
P0Legend = [P0Legend; '$\Delta P_{0,set-meas}$'];

% Delta P0 Goal Limits
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/6)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/6) - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(1:length(MATDSS.Time.Sim.TimeSpan)/6)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(1:length(MATDSS.ControlSignals.P0Set)/6) - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/6:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/6:end) - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan(length(MATDSS.Time.Sim.TimeSpan)/6:end)-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set(length(MATDSS.ControlSignals.P0Set)/6:end) - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
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


leg = legend(P0ax,'','','','','','1CA','$\rm{6CA}_{RoU = 1}$', '$\rm{6CA}_{RoU = 5}$','interpreter','latex','NumColumns',1, 'Position', [0.468 0.805 0.1 0.05]);
    
leg.ItemTokenSize = [10,15];
%% Customize the plot

% General Results
% myxlim = [0 25];
% set(P0ax, "XLim", myxlim)
% set(Iax, "XLim", myxlim)
% set(Vax, "XLim", myxlim)

% set(P0ax, "XTickMode", "auto")
% set(Vax, "XTickMode", "auto")
% set(Iax, "XTickMode", "auto")

% P0axylim = [-275 0];
% Vaxylim = [0.946 0.965];
% Iaxylim = [112 135];
% set(P0ax, "YLim", P0axylim)
% set(Vax, "YLim", Vaxylim)
% set(Iax, "YLim", Iaxylim)


% VaxYTicks = [0.95, 0.96];
% set(Vax, 'YTick', VaxYTicks)%, 'YTickLabelRotation', 80)

% IaxYTicks = [115:15:135];
% set(Iax, 'YTick', IaxYTicks)




%% MyPositionPlot = [0 0 12.80/1.5 7.20/1.5];
% set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)

figname = 'IEEE123RampTrackingWithDisturbance';

% Use this to export eps and force vector format
exportgraphics(h, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (h,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(h, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(h,'landscape');
set(h, 'Renderer', 'painters');
exportgraphics(h, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')


% set(h, 'color', 'none');
% print (h,'-depsc', '-r300', [CaseFolder '\' figname 'testingtest.eps'], '-tight')


% stackedplot(copyobj(figure,P0ax))




print(h, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');




%{
%% Prepare to Show current control on the right
set (h, 'Units', 'normalized', 'Position', [0.25,0.25,0.05,0.4],'Visible', 'on');

set(P0ax, "YLabel", [])
set(Iax, "YLabel", [])
set(Vax, "YLabel", [])

myxlim = [10 30];
set(P0ax, "XLim", myxlim)
set(Iax, "XLim", myxlim)
set(Vax, "XLim", myxlim)

set(P0ax, "XTickMode", "auto")
set(Vax, "XTickMode", "auto")
set(Iax, "XTickMode", "auto")

P0axylim = [-215 -185];
Vaxylim = [0.946 0.965];
Iaxylim = [123 126];
set(P0ax, "YLim", P0axylim)
set(Vax, "YLim", Vaxylim)
set(Iax, "YLim", Iaxylim)


VaxYTicks = [0.95, 0.96];
set(Vax, 'YTick', VaxYTicks)%, 'YTickLabelRotation', 80)

% IaxYTicks = [115:15:135];
set(Iax, 'YTickMode', 'auto')

legend(P0ax,'off')

% MyPositionPlot = [0 0 12.80/1.5 7.20/1.5];

% set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)



for i = 1:2

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

pause(1)
end

%
figname = '5BusStepTrackingWithDisturbance_CurrentControl';

% Use this to export eps and force vector format
exportgraphics(h, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (h,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(h, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(h,'landscape');
exportgraphics(h, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(h, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');

% stackedplot(copyobj(figure,P0ax))


% Generate inset figures (zoomed-in inset figures to show the curves)

% Current zoomed in figure
close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
Izax = axes(g);
hold on;
ILim = 125; %Current limit


plot(Izax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat(ILim,size(MATDSS.Sim.Meas.at)), '--', 'Color',[MyColors(end,:)], 'linewidth', 2);

for i = 1:length(MyPlotsData) - ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    Ibranches = {'l3.1','l3.2','l3.3'};
    for j = 1:length(Ibranches)
        IBranchName = Ibranches{j};
        IBranchIndex = find(strcmp(AllMeasBranchNames,IBranchName));
        if IBranchIndex > 0
            plot(Izax,MATDSS.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Meas.IProfile(IBranchIndex,:),'-','LineWidth',1,'Color', MyColors(i,:));
        end
    end
end
grid(Izax,'on')
box(Izax,'on')
set(Izax, 'FontSize', 18);
hold off;


set(Izax, "XLim", [0.7 2.5], "YLim", [114 118])


figname = '5BusStepTrackingWithDisturbance_Current1';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');

% stackedplot(copyobj(figure,P0ax))



% Voltage zoomed in figure

close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
Vzax = axes(g);
hold on
VLim = 0.05; %Voltage limits

plot(Vzax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 - VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);
plot(Vzax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 + VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(end,:), 'linewidth', 2);

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    Vnodes = {'n5.1','n5.2','n5.3'};
    Legendflag = true;
    for j = 1:length(Vnodes)
        VNodeName = Vnodes{j};
        VNodeIndex = find(strcmp(AllNodeNames,VNodeName));
        if VNodeIndex > 0
            plot(Vzax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,MATDSS.Sim.Meas.VMagProfilePu(VNodeIndex,:),'-','LineWidth',1,'Color', MyColors(i,:));
        end
    end
end
grid(Vzax,'on')
box(Vzax,'on')

set(Vzax, 'FontSize', 18);
set(Vzax, "XLim", [0.6 2.5], "YLim", [0.95 0.9535])
hold off;


figname = '5BusStepTrackingWithDisturbance_Voltage1';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');

% stackedplot(copyobj(figure,P0ax))

%


%}
% Power zoomed in figure1
close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
P0zax = axes(g);
hold on;



% Delta P0 set-points
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));

% Delta P0 Goal Limits
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    plot(P0zax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
end
grid on;
box(P0zax,'on')
set(P0zax, 'FontSize', 18);
set(P0zax, "XLim", [13 15], "YLim", [-320 -220])
hold off;
%%
figname = 'IEEE123RampTrackingWithDisturbance_Power1';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');


% Power zoomed in figure1
close(figure(2))
g = figure(2);

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (g, 'Units', 'normalized', 'Position', [0.75,0.25,0.1,0.1],'Visible', 'on');


%%% I_Plot
P0zax = axes(g);
hold on;



% Delta P0 set-points
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));

% Delta P0 Goal Limits
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
% plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(1))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end-1,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));
plot(P0zax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE(2))./1e3, '--', 'LineWidth', 1, 'Color',MyColors(end,:));

for i = 1:length(MyPlotsData)-ncase
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    plot(P0zax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
end
grid on;
box(P0zax,'on')
set(P0zax, 'FontSize', 18);
set(P0zax, "XLim", [43 45], "YLim", [-360 -260])
hold off;

figname = 'IEEE123RampTrackingWithDisturbance_Power2';

% Use this to export eps and force vector format
exportgraphics(g, [CaseFolder '\' figname '.eps'], 'BackgroundColor','none', 'ContentType','vector')

% Png
print (g,'-dpng', '-r300', [CaseFolder '\' figname '.png'])

saveas(g, [CaseFolder '\' figname '.fig'])

% Pdf vector format
orient(g,'landscape');
exportgraphics(g, [CaseFolder '\' figname '.pdf'], 'BackgroundColor','none', 'ContentType','vector')
print(g, [CaseFolder '\' figname 'vecText.pdf'], '-dpdf', '-painters', '-r600');

