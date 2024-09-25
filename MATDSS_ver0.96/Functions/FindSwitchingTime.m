% Open the figure with Delta P0 plot

% This code will find the switching time and report it back to you below.
% Default limit is +- 5kW
Limit = 5;


h = gca;

MyLines = findobj(h, 'Type', 'Line');

MyLinesNames = {MyLines.DisplayName}';

MyLineIndex = find(strcmp(MyLinesNames,'$\Delta P_{0}$' ));

TargetLine = MyLines(MyLineIndex);
XData = TargetLine.XData;
YData = TargetLine.YData;
t_index = -1;

YData_Flag = find(abs(YData) >= Limit);

display(['Settling Time for ' num2str(Limit) 'kW = ' num2str(XData(YData_Flag(end)+1)-1) 's'])

close all