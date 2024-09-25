clear;
clc;
close all;


CaseFolder = uigetdir("D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Exports");     % initial subdirectory
% CaseFolder = uigetdir("C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.84\Exports");     % initial subdirectory
[~,CaseFiles] = dos(['dir /s /b ' '"' fullfile(CaseFolder,'*.mat') '"' ]);  % OS dir command
CaseFiles=textscan(CaseFiles,'%s','delimiter','\n'); CaseFiles=CaseFiles{:};
% addpath {'D:\OneDrive - University of Waterloo\PhD UoW\2023\03 March\MATDSS_ver0.83\Functions\distinguishable_colors.m'}
% CaseFolder = "D:\OneDrive - University of Waterloo\PhD UoW\2023\03 March\MATDSS_ver0.83\Exports\MATDSS_Simple2Area_test_0404232312";

CaseNumber = strfind(CaseFolder,"\");
CaseNumber = CaseFolder(CaseNumber(end)+1);

CaseNumber(strfind(CaseNumber," ")) = "_";
MyPlotsData =[];
for i = 1:length(CaseFiles)
    MyPlotsData = [MyPlotsData;{load(CaseFiles{i})}];
end

MyColors = MyPlotsData{1}.Plot.Colors;

MyPositionPlot = [0 0 12.80 7.20];
VoltageZoomRange = [0.943 0.97];
CurrentRange = [120 145];
%% Plotting Delta P0
close(figure(1))
h = figure(1);
% If you want to make the figure unvisible to enahnce the experience, use the next line
% Note: This might lead to .fig file being invisible when you try to open it in MATLAB after that,
% so make sure you use the function below near .fig export line to reset the visibility when you reopen the file.
% ExpFig = figure('Name', 'ExpFig','Visible', 'off');

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (h, 'Units', 'normalized', 'Position', [0.25,0.25,0.5,0.5],'Visible', 'on');


axes('FontSize', 16);
hold on
% color_array = distinguishable_colors(NCurves);
P0Legend = {};
P0ax = h.CurrentAxes;
ContE = [5000, 10000]; %5kW limits

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
plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) - ContE)./1e3, '--', 'LineWidth', 2, 'Color',MyColors(11,:));
plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1) + ContE)./1e3, '--', 'LineWidth', 2, 'Color',MyColors(11,:));
P0Legend = [P0Legend;"";""];
P0Legend = [P0Legend;"";""];

for i = 1:length(MyPlotsData)
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
    P0Legend = [P0Legend;['$\Delta P_{0,\ \textrm{Config. }' num2str(i) '}$']];
end
grid on;
box(P0ax,'on')
xlabel(P0ax,'$t \ (s)$');
ylabel(P0ax,'$\Delta P_{0} \ (kW)$');
legend(P0ax, P0Legend, 'interpreter', 'latex','Location','northeast','NumColumns',2);
set(P0ax, 'FontSize', 14);
% title(P0ax,['Global G = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
xlim(P0ax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
% hold off;


%{


xlim(P0ax, [0,20])
xticks(P0ax, 0:2:20)
% xtickangle(P0ax, 00)
legend(P0ax, '$P_{0,\sf{set}}$','','','1CA', '2CA', 'No-Control', 'interpreter', 'latex','Location','northeast','NumColumns',2)




% P0ax2 = axes();
% plot(P0ax2,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
% semilogy(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'-','LineWidth',1,'Color', MyColors(i,:));
% ylim(P0ax2, [0 120])
% linkaxes([P0ax, P0ax2],'x') 
ylim(P0ax, [-270 0])
legend(P0ax, '$P_{0,\sf{set}}$','','','1CA', '2CA', 'No-Control', 'interpreter', 'latex','Location','northeast','NumColumns',2)
% breakyaxis([0, 90], 0.01)
xlim(P0ax, [0,10])
% xticks(P0ax, 0:2:20)
% xtickangle(P0ax, 00)

breakyaxis([200, 800]);

%}
set(findall(h,'-property','FontSize'),'FontSize', 16);

% ylim(P0ax,[-35 110])

set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)
%{
print (h,'-dpng', '-r300', [CaseFolder '\' CaseNumber '_' char(916) 'P0.png'])
% print(h, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
saveas(h, [CaseFolder '\' CaseNumber '_'  char(916) 'P0.fig'])
saveas(h, [CaseFolder '\' CaseNumber '_'  char(916) 'P0.eps'],'epsc')
saveas(h, [CaseFolder '\' CaseNumber '_'  'DeltaP0.eps'],'epsc')
orient(h,'landscape');
saveas(h, [CaseFolder '\' CaseNumber '_'  char(916) 'P0.pdf'])
%}
% clear h
% clear P0Legend legend
% clear P0ax
% close all
%% Plotting Voltage
clear VLegend
close(figure(2))
h = figure(2);
% If you want to make the figure unvisible to enahnce the experience, use the next line
% Note: This might lead to .fig file being invisible when you try to open it in MATLAB after that,
% so make sure you use the function below near .fig export line to reset the visibility when you reopen the file.
% ExpFig = figure('Name', 'ExpFig','Visible', 'off');

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (h, 'Units', 'normalized', 'Position', [0.25,0.25,0.5,0.5],'Visible', 'on');


axes('FontSize', 16);
hold on
% color_array = distinguishable_colors(NCurves);
VLegend = {};
Vax = h.CurrentAxes;
VLim = 0.05; %Voltage limits

MATDSS = MyPlotsData{1}.MATDSS;
% Plotting 1-time curves
% 
% % Delta P0 set-points
% plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));
% P0Legend = [P0Legend; '$\Delta P_{0,set-meas}$'];

% Delta P0 Goal Limits

plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 - VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(11,:), 'linewidth', 2);
plot(Vax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat((1 + VLim),size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(11,:), 'linewidth', 2);
VLegend = [VLegend;""];
VLegend = [VLegend;""];

for i = 1:length(MyPlotsData)
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    Vnodes = {'n5.1','n5.2','n5.3'};
    Legendflag = true
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
xlabel(Vax,'$t \ (s)$');
ylabel(Vax,'$\textrm{n}5 \ V \ (p.u.)$');
legend(Vax,VLegend, 'interpreter', 'latex','Location','northeast','NumColumns',2);
set(Vax, 'FontSize', 14);
% title(Vax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
xlim(Vax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
% ylim
hold off;

set(findall(h,'-property','FontSize'),'FontSize', 16);



set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)

%{

xlim(Vax, [0,20])
ylim(Vax, VoltageZoomRange)

xticks(Vax, 0:2:20)
% xtickangle(P0ax, 00)
legend(Vax, '','$\underbar{V}$','1CA','', '', '2CA', '', '', 'No-Control', 'interpreter', 'latex','Location','northeast','NumColumns',2)

%}


%{

print (h,'-dpng', '-r300', [CaseFolder '\' CaseNumber '_Voltage.png'])
% print(h, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
saveas(h, [CaseFolder '\' CaseNumber '_Voltage.fig'])
saveas(h, [CaseFolder '\' CaseNumber '_Voltage.eps'],'epsc')
orient(h,'landscape');
saveas(h, [CaseFolder '\' CaseNumber '_Voltage.pdf'])

legend(Vax,VLegend, 'interpreter', 'latex','Location','northeast','NumColumns',2);
ylim(Vax, VoltageZoomRange)

set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)
print (h,'-dpng', '-r300', [CaseFolder '\' CaseNumber '_Voltage_zoom.png'])
% print(h, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
saveas(h, [CaseFolder '\' CaseNumber '_Voltage_zoom.fig'])
saveas(h, [CaseFolder '\' CaseNumber '_Voltage_zoom.eps'],'epsc')
orient(h,'landscape');
saveas(h, [CaseFolder '\' CaseNumber '_Voltage_zoom.pdf'])

clear h
clear Vax
clear VLegend
%}
%% Plotting Current
clear ILegend
close(figure(3))
h = figure(3);
% If you want to make the figure unvisible to enahnce the experience, use the next line
% Note: This might lead to .fig file being invisible when you try to open it in MATLAB after that,
% so make sure you use the function below near .fig export line to reset the visibility when you reopen the file.
% ExpFig = figure('Name', 'ExpFig','Visible', 'off');

% h = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Figure_ID = Figure_ID + 1;
set (h, 'Units', 'normalized', 'Position', [0.25,0.25,0.5,0.5],'Visible', 'on');


axes('FontSize', 16);
hold on
% color_array = distinguishable_colors(NCurves);
ILegend = {};
Iax = h.CurrentAxes;
ILim = 125; %Current limit

MATDSS = MyPlotsData{1}.MATDSS;
% Plotting 1-time curves
% 
% % Delta P0 set-points
% plot(P0ax,MATDSS.Time.Sim.TimeSpan-MATDSS.Time.Sim.ST,(MATDSS.ControlSignals.P0Set - MATDSS.ControlSignals.P0Set(1))./1e3,'--','LineWidth',2,'Color',MyColors(10,:));
% P0Legend = [P0Legend; '$\Delta P_{0,set-meas}$'];

% Delta P0 Goal Limits

plot(Iax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,repmat(ILim,size(MATDSS.Sim.Meas.at)), '--', 'Color',MyColors(11,:), 'linewidth', 2);
% plot(Iax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(1 + ContE), '--', 'Color',MyColors(11,:));
ILegend = [ILegend;""];
% ILegend = [ILegend;""];

for i = 1:length(MyPlotsData)
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
xlabel(Iax,'$t \ (s)$');
ylabel(Iax,'$\textrm{L}3 \ I \ (A)$');
legend(Iax,ILegend, 'interpreter', 'latex','Location','southeast','NumColumns',2);
set(Iax, 'FontSize', 14);
% title(Iax,['Global Gain = ' num2str(MATDSS.Cont.Gain)],'Interpreter','latex')
xlim(Iax,[0, MATDSS.Time.Sim.TimeSpan(end)-MATDSS.Time.Sim.ST])
% ylim
hold off;

set(findall(h,'-property','FontSize'),'FontSize', 16);



set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)

ylim(Iax, CurrentRange)


%{

xlim(Iax, [0,20])
% ylim(Iax, VoltageZoomRange)

xticks(Iax, 0:2:20)
% xtickangle(P0ax, 00)
legend(Iax, '$\overline{\rm{i}}$','1CA','','','2CA','','','No-Control', 'interpreter', 'latex','Location','northeast','NumColumns',2)

%}


%{


print (h,'-dpng', '-r300', [CaseFolder '\' CaseNumber '_Current.png'])
% print(h, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
saveas(h, [CaseFolder '\' CaseNumber '_Current.fig'])
saveas(h, [CaseFolder '\' CaseNumber '_Current.eps'],'epsc')
orient(h,'landscape');
saveas(h, [CaseFolder '\' CaseNumber '_Current.pdf'])


% ylim(Iax, [0.945 0.965])
% 
% set(h,'PaperUnits','inches','PaperPosition',MyPositionPlot)
% print (h,'-dpng', '-r300', [CaseFolder '\' CaseNumber '_Current_zoom.png'])
% % print(h, '-dpng', '-r300', [OutputFileDir '\' app.PlotTabGroup.Children(i).Title '.png'])
% saveas(h, [CaseFolder '\' CaseNumber '_Current_zoom.fig'])
% saveas(h, [CaseFolder '\' CaseNumber '_Current_zoom.eps'])
% orient(h,'landscape');
% saveas(h, [CaseFolder '\' CaseNumber '_Current_zoom.pdf'])

clear h
clear Iax
clear ILegend
%}
disp('done')

%% Switching time results

switchingtimes = {};
switchingtimes10 = {};
switchingtimes5 = {};

Limit = ContE(1)./1e3;
YData_Flag_10_All = {};
YData_Flag_5_All = {};
for i = 1:length(MyPlotsData)
    MATDSS = MyPlotsData{i}.MATDSS;
    DER = MyPlotsData{i}.DER;
    % plot(P0ax,MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST,(sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(1))./1e3,'LineWidth',2,'Color', MyColors(i,:));
    t_index = -1;
    YData = (sum(MATDSS.Sim.Meas.P0) - MATDSS.ControlSignals.P0Set(11))./1e3;
    XData = MATDSS.Sim.Meas.at-MATDSS.Time.Sim.ST;

    YData_Flag = find(abs(YData) >= Limit);
    YData_Flag_10 = find(abs(YData) >= 10);
    YData_Flag_5 = find(abs(YData) >= 5);

    YData_Flag_5 = XData(YData_Flag_5+1);
    YData_Flag_10 = XData(YData_Flag_10+1);
    YData_Flag_10_All = [YData_Flag_10_All, YData_Flag_10'];
    YData_Flag_5_All = [YData_Flag_5_All, YData_Flag_5'];
    
    if YData_Flag(end) + 1 > length(XData)
        YData_Flag(end) = length(XData) - 1;
    end


    switchingtimes = [switchingtimes;...
                     {['case ' num2str(i)]}, {XData(YData_Flag(end)+1)}];

    YData_Flag = find(abs(YData) >= 10);
    YData_Flag(YData_Flag > 500) = [];
    switchingtimes10 = [switchingtimes10;...
            {['case ' num2str(i)]}, {XData(YData_Flag(end)+1)}];

    YData_Flag = find(abs(YData) >= 5);
    YData_Flag(YData_Flag > 2000) = [];
    switchingtimes5 = [switchingtimes5;...
            {['case ' num2str(i)]}, {XData(YData_Flag(end)+1)}];
    % P0Legend = [P0Legend;['$\Delta P_{0,\ \textrm{Config. }' num2str(i) '}$']];
end


cell2table([switchingtimes5, switchingtimes10, switchingtimes])



