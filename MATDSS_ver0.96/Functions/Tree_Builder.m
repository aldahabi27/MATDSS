%% This code is developed to build a tree of the network/feeder - IEEE 8500 in particular now

clear;
clc;
close all;

addpath Functions
% This code intiate the link with OpenDSS engine to run the specified
% OpenDSS file
MATDSS.Sim.DSSCompile = ['Compile "C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\DSS_Files\8500-Node\MATDSS_Master_1CA.dss"'];
% MATDSS.Sim.DSSCompile = ['Compile "D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\DSS_Files\8500-Node\MATDSS_Master_1CA.dss"'];

% Initiate DSS link
MATDSS.Sim.DSSObj = actxserver('OpenDSSEngine.DSS');

% Flag if link is established
MATDSS.Sim.DSSStartOk = MATDSS.Sim.DSSObj.Start(0);

if MATDSS.Sim.DSSStartOk
    % Define the Text interface (to run commands from MATLAB)
    MATDSS.Sim.DSSText = MATDSS.Sim.DSSObj.Text;

    % Compile the DSS file first to load the circuit
    MATDSS.Sim.DSSText.command = MATDSS.Sim.DSSCompile;

    % Set up the interface variables DSSCircuit and DSSSolution
    MATDSS.Sim.DSSCircuit = MATDSS.Sim.DSSObj.ActiveCircuit;
    MATDSS.Sim.DSSSolution = MATDSS.Sim.DSSCircuit.Solution;
    %     MATDSS.Sim.DSSLines
else
    msgbox('DSSEngine Failed to start!');
    disp({'!****************************************!';...
        '!* Could not initiate link with OpenDSS *!';...
        '!****************************************!';...
        ' ';'Simulation Aborted!';''});
    return
end


for i = 1:20
    MATDSS.Sim.DSSText.Command = 'solve';
end
%% Get organized list of buses with parent and children lists

% First we generate a list of all buses and their parent and children buses
% This won't be ordered yet. We build the proper tree structure next when
% we work on organizing the tree.


% Get all circuit elements and buses names
AllCktElements = MATDSS.Sim.DSSCircuit.AllElementNames;
AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames;

% This is my list variable
bus_connections = cell(length(AllBusNames), 5);


% we loop over the buses and check connected elements to each bus and what
% are the parent and children buses that connected to it.
for i = 1:length(AllBusNames)
    if mod(round(100*i/length(AllBusNames)),10) == 0 && round(100*i/length(AllBusNames)) ~= round(100*(i-1)/length(AllBusNames))
        disp([num2str(round(100*i/length(AllBusNames))), '%'])
    end

    % Select bus i
    myBus = MATDSS.Sim.DSSCircuit.SetActiveBus(AllBusNames{i});
    myBus = MATDSS.Sim.DSSCircuit.ActiveBus;

    % Get list of elements that are connected to bus i (AllPCE are all
    % power convergence elements such as loads...etc. and AllPDE are power
    % delievery elements such as lines and transformers
    ListOfElements = [myBus.AllPCEatBus; myBus.AllPDEatBus];
    % Clear any empty strings returned from OpenDSS
    ListOfElements(strcmpi(ListOfElements,'')) = [];
    % remove any "none" entries too
    ListOfElements(strcmpi(ListOfElements,'none')) = [];

    % Here we prepare the variables that we will use to organize the
    % neighbours of bus i (devices + buses)
    LeftBuses = {}; %higher order
    RightBuses = {};% lower order

    % Power flows from left to right



    % elements connected to left and right buses
    LeftElements = {};
    RightElements = {};


    MyOtherBuses = {}; % This is used to avoid having any repeated buses

    %clear any repeated element as OpenDSS might report multiple instances
    %of elements!
    L = length(ListOfElements);
    j = 1;
    while j < L - 1
        temp = ListOfElements{end-j+1};
        if sum(strcmpi(temp, ListOfElements(1:end-j)))
            ListOfElements(end-j+1) = [];
            L = L - 1;
        end
        j = j + 1;
    end

    % Check the unique list of elements one-by-one
    for j = 1:length(ListOfElements)

        % Select element j
        myElement = MATDSS.Sim.DSSCircuit.SetActiveElement(ListOfElements{j});
        myElement = MATDSS.Sim.DSSCircuit.ActiveElement;


        % myElement.Properties(1)


        % get list of buses connected to element j (usually this is ordered
        % in the power flow direction) and we rely on this to categorize
        % the buses
        myElementBuses = myElement.BusNames;

        % check reported buses
        for k = 1:length(myElementBuses)
            % Since some elements are single-phase, we care about the full
            % bus name without phases, therefore, we omit any phases in our
            % naming here
            testbus = myElementBuses{k};
            idot = strfind(testbus, '.');
            if ~isempty(idot)
                testbus = testbus(1:idot-1);
            end


            % check whether it is right or left bus (by checking if it has
            % k = 1 or k > 1
            if ~(strcmpi(testbus, AllBusNames{i}) || sum(strcmpi(MyOtherBuses,testbus))) && k==1
                LeftBuses = [LeftBuses;testbus];
                LeftElements = [LeftElements; ListOfElements{j}];
                MyOtherBuses = [MyOtherBuses, testbus];
            elseif ~(strcmpi(testbus, AllBusNames{i}) || sum(strcmpi(MyOtherBuses,testbus))) && k~=1
                RightBuses = [RightBuses;testbus];
                RightElements = [RightElements; ListOfElements{j}];
                MyOtherBuses = [MyOtherBuses, testbus];
            end


        end
    end


    bus_connections(i,:) = [{LeftBuses, LeftElements, AllBusNames{i}, RightElements ,RightBuses}];

end

bus_connections(1,1:2) = [{''}, {''}];


% For special cases, where multiple lines (swtich + line) are connected to
% a bus, we double that entry and make sure we have single parent bus per
% row in the list

disp('Building the extended list')
extended_bus_connections = {};
i = 1;
while true
    if mod(round(100*i/size(bus_connections,1)),10) == 0 && round(100*i/size(bus_connections,1)) ~= round(100*(i-1)/size(bus_connections,1))
        disp([num2str(round(100*i/size(bus_connections,1))), '%'])
    end
    % if the row has multiple parents
    if size(bus_connections{i,1},1) > 1
        val1 = bus_connections{i,1}; %get the list of parent buses
        val2 = bus_connections{i,2}; % get hte list of connections to parent buses
        for j = 1:length(val1)  % save the new entries
            extended_bus_connections = [extended_bus_connections; [val1(j), val2(j), bus_connections(i,3:end)]];
        end
    else %this is normal row, just copy it as is
        extended_bus_connections = [extended_bus_connections; bus_connections(i,:)];
    end
    i = i + 1; %keep looping
    if i > size(bus_connections,1) % if we go beyond the size of original list, exit
        break;
    end


end


%% Building the tree

% We now start from the source bus, and then move down the feeder to
% organize the tree.

disp('building the tree');
specialnodes = {};
AdjMatrix = zeros(length(AllBusNames));
AllBusesOrignal = string(AllBusNames);
Parents = string(extended_bus_connections(:,1));
AllBusesExtended = string(extended_bus_connections(:,3));
for i = 1:length(AllBusNames)
    if mod(round(100*i/length(AllBusNames)),10) == 0 && round(100*i/length(AllBusNames)) ~= round(100*(i-1)/length(AllBusNames))
        disp([num2str(round(100*i/length(AllBusNames))), '%'])
    end

    CurrentBusIndex = find(strcmpi(AllBusesExtended,{AllBusNames{i}}));
    CurrentBusParents = string({extended_bus_connections{CurrentBusIndex,1}});

    for j = 1:length(CurrentBusParents)
        ParentIndex = find(strcmpi(AllBusesOrignal,CurrentBusParents{j}));
        % AdjMatrix(i,ParentIndex) = 1;
        AdjMatrix(ParentIndex,i) = 1;
    end
end
%%
close all
G = digraph(AdjMatrix,AllBusNames);

h = figure(1);
subplot(1,2,1);
p1 = plot(G,'NodeLabel',G.Nodes.Name);
hold on;
highlight(p1,1,'nodecolor','r','MarkerSize',10);
subplot(1,2,2);
p2 = plot(G,'NodeLabel',1:length(AllBusNames));
hold on;
highlight(p2,1,'nodecolor','r','MarkerSize',10);

h2 = figure(2);
p3 = plot(G);
hold on;
highlight(p3,1,'NodeColor','r','MarkerSize',10);


for i = 1:size(bus_connections,1)
    if size(bus_connections{i,1},1) > 1
        node_index = find(strcmpi(AllBusesOrignal,bus_connections(i,3)));
        highlight(p3,node_index, 'NodeColor','g','MarkerSize',10)
        specialnodes = [specialnodes;{bus_connections(i,3)}];
    end
end
%%

descendants = [];
node_name = 'e182744';

node_index = find(strcmpi(AllBusesOrignal,node_name));
descendants = collectDescendants(descendants, node_index,G);

descendants_names = AllBusesOrignal(descendants);

myDescendantstable = table(descendants_names, descendants)


%% Plot circuit in Matlab
close(figure(3))
h3 = figure(3);
hold on
A = readlines("C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.84\DSS_Files\8500-Node\Buscoords.dss");
% A = readlines("D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.84\DSS_Files\8500-Node\Buscoords.dss");

bus_coordinates = {};
for i = 4:size(A,1)-1
    busline = convertStringsToChars(A(i));
    icomma = strfind(busline,',');
    if ~isempty(icomma)
        busname = busline(1:icomma(1)-1);
        busx = str2double(busline(icomma(1)+1:icomma(2)-1));
        busy = str2double(busline(icomma(2)+1:end));
        bus_coordinates = [bus_coordinates; {busname, busx, busy}];
    end
end

plot([bus_coordinates{:,2}], [bus_coordinates{:,3}], 'o')
offset = 100;
text([bus_coordinates{:,2}]+offset, [bus_coordinates{:,3}]+offset,bus_coordinates(:,1),'fontsize',8)
specialnodes_coord = {};
for i = 1:length(specialnodes)
    node_index = find(strcmpi(bus_coordinates(:,1),specialnodes{i}));
    plot([bus_coordinates{node_index,2}], [bus_coordinates{node_index,3}], 'og', 'MarkerFaceColor','g')
    specialnodes_coord = [specialnodes_coord; {bus_coordinates{node_index,2}, bus_coordinates{node_index,3}}];
end


%%
findbus = 'l3081380';
find(strcmpi(bus_connections(:,3),findbus))


%%
% 
% bus_connections(1714,:) = {bus_connections{1714,1}(1), bus_connections{1714,2}(1), bus_connections{1714,3:end}}; % clean 1714
% bus_connections(823,:) = {bus_connections{823,1}(1), bus_connections{823,2}(1), bus_connections{823,3:end}}; % clean 823
% bus_connections(2298,:) = {bus_connections{2298,1}(1), bus_connections{2298,2}(1), bus_connections{2298,3:end}}; % clean 2298
% bus_connections(1563,:) = {bus_connections{1563,1}(1), bus_connections{1563,2}(1), bus_connections{1563,3:end}}; % clean 1563
% bus_connections(1774,:) = {bus_connections{1774,1}(2), bus_connections{1774,2}(2), bus_connections{1774,3:end}}; % clean 1774

%% Building the tree

old_extended_bus_connections  = extended_bus_connections;

extended_bus_connections = bus_connections;
% We now start from the source bus, and then move down the feeder to
% organize the tree.

disp('building the tree');
specialnodes = {};
AdjMatrix = zeros(length(AllBusNames));
AllBusesOrignal = string(AllBusNames);
Parents = string(extended_bus_connections(:,1));
AllBusesExtended = string(extended_bus_connections(:,3));
for i = 1:length(AllBusNames)
    if mod(round(100*i/length(AllBusNames)),10) == 0 && round(100*i/length(AllBusNames)) ~= round(100*(i-1)/length(AllBusNames))
        disp([num2str(round(100*i/length(AllBusNames))), '%'])
    end

    CurrentBusIndex = find(strcmpi(AllBusesExtended,{AllBusNames{i}}));
    CurrentBusParents = string({extended_bus_connections{CurrentBusIndex,1}});

    for j = 1:length(CurrentBusParents)
        ParentIndex = find(strcmpi(AllBusesOrignal,CurrentBusParents{j}));
        % AdjMatrix(i,ParentIndex) = 1;
        AdjMatrix(ParentIndex,i) = 1;
    end
end
%%
close all
G = digraph(AdjMatrix,AllBusNames);

h = figure(1);
subplot(1,2,1);
p1 = plot(G,'NodeLabel',G.Nodes.Name);
hold on;
highlight(p1,1,'nodecolor','r','MarkerSize',10);
subplot(1,2,2);
p2 = plot(G,'NodeLabel',1:length(AllBusNames));
hold on;
highlight(p2,1,'nodecolor','r','MarkerSize',10);

h2 = figure(2);
p3 = plot(G);
hold on;
highlight(p3,1,'NodeColor','r','MarkerSize',10);


for i = 1:size(bus_connections,1)
    if size(bus_connections{i,1},1) > 1
        node_index = find(strcmpi(AllBusesOrignal,bus_connections(i,3)));
        highlight(p3,node_index, 'NodeColor','g','MarkerSize',10)
        specialnodes = [specialnodes;{bus_connections(i,3)}];
    end
end

%% Area 2
% descendants = [];
% node_name = 'l3081380';
% 
% node_index = find(strcmpi(AllBusesOrignal,node_name));
% descendants = collectDescendants(descendants, node_index,G);
% 
% descendants_names = AllBusesOrignal(descendants);
% 
% myDescendantstable = table(descendants_names, descendants);
% 
% Area2 = descendants_names;

%% Area 3
% descendants = [];
% node_name = 'l2955077';
% 
% node_index = find(strcmpi(AllBusesOrignal,node_name));
% descendants = collectDescendants(descendants, node_index,G);
% 
% descendants_names = AllBusesOrignal(descendants);
% 
% myDescendantstable = table(descendants_names, descendants);
% 
% Area3 = descendants_names;

%% Area 1

% descendants = [];
% node_name = 'sourcebus';
% 
% node_index = find(strcmpi(AllBusesOrignal,node_name));
% descendants = collectDescendants(descendants, node_index,G);
% 
% descendants_names = AllBusesOrignal(descendants);
% 
% myDescendantstable = table(descendants_names, descendants)

%
% Area1 = setdiff(AllBusesOrignal,Area2);
% Area1 = setdiff(Area1, Area3);
% 

%% Save areas in excel
% AreasData = repmat("",max([length(Area1), length(Area2), length(Area3)]),1);
% 
% AreasData(1:length(Area1),1) = Area1;
% AreasData(1:length(Area2),2) = Area2;
% AreasData(1:length(Area3),3) = Area3;
% 
% AreasTable = cell2table(cellstr(AreasData));
% AreasTable.Properties.VariableNames = {'Area1', 'Area2', 'Area3'};
% 
% 
% writetable(AreasTable, 'IEEE8500AreasTable.xlsx');
%% Plot circuit in Matlab
close(figure(3))
h3 = figure(3);
hold on
A = readlines("C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\DSS_Files\8500-Node\Buscoords.dss");
% A = readlines("D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\DSS_Files\8500-Node\Buscoords.dss");

bus_coordinates = {};
for i = 4:size(A,1)-1
    busline = convertStringsToChars(A(i));
    icomma = strfind(busline,',');
    if ~isempty(icomma)
        busname = busline(1:icomma(1)-1);
        busx = str2double(busline(icomma(1)+1:icomma(2)-1));
        busy = str2double(busline(icomma(2)+1:end));
        bus_coordinates = [bus_coordinates; {busname, busx, busy}];
    end
end

plot([bus_coordinates{:,2}], [bus_coordinates{:,3}], 'o')
offset = 10;
text([bus_coordinates{:,2}]+offset, [bus_coordinates{:,3}]+offset,bus_coordinates(:,1),'fontsize',8)
specialnodes_coord = {};
specialnodes = {'m1142843', 'l3081380', 'l2955077', 'm1166374'};
for i = 1:length(specialnodes)
    node_index = find(strcmpi(bus_coordinates(:,1),specialnodes{i}));
    plot([bus_coordinates{node_index,2}], [bus_coordinates{node_index,3}], 'og', 'MarkerFaceColor','g')
    specialnodes_coord = [specialnodes_coord; {bus_coordinates{node_index,2}, bus_coordinates{node_index,3}}];
end

%% Work on June 20 2023

% Defining the headbuses for our control areas. The idea is to identify the
% head buses, and then we let the code generate the list of buses for each
% area by looking at the children of those head buses. Then, it shall
% remove all buses of all other areas from each area to get teh list of
% unique buses that belongs to the area alone. This is not optimium
% approach, but it shall do the job we want.

addpath FuzzySearch 

InterfaceBuses = {"sourcebus";...
                  "m1209791";...
                  "l3160863";...
                  "l3253570";... % if this does not work properly, change it to "I2937757"
                  "p901951";...
                  "m3032977";...
                  "m1166377";...
                  "l3081380";...
                  "l2955077";...
                  "p829965";...
                  "l2861218";...
                  "n1230121";...
                  "l2933164";...
                  "p827521";...
                  "p827520";...
                  "m1069184";...
                  "p827508";...
                  "p827533";...
                  "l3104134";...
                  "e182748";...
                  "p827531";...
                  "l2973162";...
                  "l3067478";...
                  "e183472";...
                  "e183473";...
                  "p901946";...
                  "p829796";...
                  "m1108535";...
                  "n1134479";...
                  "p827548";...
                  "n1136663";...
                  "e206209";...
                  "m1047484";...
                  "p827563";...
                  "n1140519";...
                  "n1138594";...
                  "n1140819";...
                  "n1136666";...
                  "n1136670";...
                  "m1047766";...
                  "p827509";...
                  "p901894";...
                  "p827536";...
                  "m1026905";...
                  "m1026829";...
                  "m1009820";...
                  "m1026984";...
                  "n1134470";...
                  "e182727";...
                  };



MyAreas = struct;
MyAreas.Descendants = {};
MyAreas.Buses = {};
MyAreas.Index_Buses = {};
MyAreas.InterfaceBus = {};
MyAreas.Index_InterfaceBus = {};
MyAreas.nBuses = {};
MyAreas.Nodes = {};
MyAreas.Index_Nodes = {};
MyAreas.Buses_nPhases = [];

for i = 1:length(InterfaceBuses)
    disp(['Getting Descendants of Area ' num2str(i) ' with interface bus "' convertStringsToChars(InterfaceBuses{i}) '", (area ' num2str(i) ' of ' num2str(length(InterfaceBuses)) ')']);
    MyAreas(i).InterfaceBus = InterfaceBuses{i};
    interfacebus_index = find(strcmpi(AllBusesOrignal, MyAreas(i).InterfaceBus));
    MyAreas(i).Descendants = collectDescendants(MyAreas(i).Descendants, interfacebus_index,G);
    MyAreas(i).Index_InterfaceBus = interfacebus_index;

    MyAreas(i).Descendants_Flattended = flattenCellArray(MyAreas(i).Descendants);
end

% to onbtain the unique set of buses, I need to check first each bus if it
% is part of teh descendants of teh area, and then check if the area's
% interface bus is the closest head-parent bus

for i = 1:length(AllBusNames)
    % disp(['Bus ' num2str(i) ' of ' num2str(length(AllBusNames))])
    mybus = AllBusNames(i); % for bus i
    bus_areas = []; % areas where this bus is a descendants
    for j = 1:length(InterfaceBuses)
        if ismember(i, [MyAreas(j).Index_InterfaceBus;   MyAreas(j).Descendants_Flattended])
            bus_areas = [bus_areas; j];
        end
    end

    % now check which is the closest parent bus for this bus at hand.
    % note that parent bus for closes children area won't be considered as
    % this bus won't belong to it!
    ShortDist = inf;
    ShortDistArea = -1;
    for j = 1:length(bus_areas)
        [~, j_area_dist] = shortestpath(G, MyAreas(bus_areas(j)).InterfaceBus,mybus);
        if j_area_dist < ShortDist
            ShortDist = j_area_dist;
            ShortDistArea = bus_areas(j);
        end
    end

    MyAreas(ShortDistArea).Buses = [MyAreas(ShortDistArea).Buses; mybus];
    MyAreas(ShortDistArea).Index_Buses = [MyAreas(ShortDistArea).Index_Buses; i];    
end



%%

% checking that we got this right!
Areas_All_Buses = {};
for i = 1:length(MyAreas)
    % Check that all areas have only unique set of buses (should be, but
    % check anyway)

    if ~(isempty(setdiff(unique(MyAreas(i).Buses), MyAreas(i).Buses)))
        disp(['Area ' num2str(i) ' does not have unique set of buses']);
    end
    Areas_All_Buses = [Areas_All_Buses; MyAreas(i).Buses];
end

if isempty(setdiff(Areas_All_Buses, AllBusNames))
    disp('Looks like the work has been done correctly!')
    AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames;
    for i = 1:length(MyAreas)
        MyAreas(i).nBuses = length(MyAreas(i).Buses);
        for j = 1:MyAreas(i).nBuses
            mybus = MATDSS.Sim.DSSCircuit.Buses(MyAreas(i).Buses{j});
            mybus_nphases = mybus.NumNodes;
            mybus_iphases = mybus.Nodes;
            mybus_Nodes = {};
            mybus_iNodes = [];
            while ~isempty(mybus_iphases)
                NewNode = [mybus.Name '.' num2str(mybus_iphases(1))];
                mybus_Nodes = [mybus_Nodes; NewNode];
                mybus_iphases(1) = [];
                iNewNode = find(strcmpi(AllNodesNames, NewNode));
                mybus_iNodes = [mybus_iNodes; iNewNode];
            end
            
            MyAreas(i).Buses_nPhases = [MyAreas(i).Buses_nPhases; mybus_nphases];
            MyAreas(i).Nodes = [MyAreas(i).Nodes; mybus_Nodes];
            MyAreas(i).Index_Nodes = [MyAreas(i).Index_Nodes; mybus_iNodes];
        end
    end
else
    disp('There are buses that do not belong to any area!')
end


%% Color Coding the Control Areas on The tree and Coordination plots

% Plot circuit in Matlab
close(figure(3))
h3 = figure(3);
h3axes = axes(h3);
hold on
A = readlines("C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\DSS_Files\8500-Node\Buscoords.dss");
% A = readlines("D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\DSS_Files\8500-Node\Buscoords.dss");

bus_coordinates = {};
for i = 4:size(A,1)-1
    busline = convertStringsToChars(A(i));
    icomma = strfind(busline,',');
    if ~isempty(icomma)
        busname = busline(1:icomma(1)-1);
        busx = str2double(busline(icomma(1)+1:icomma(2)-1));
        busy = str2double(busline(icomma(2)+1:end));
        bus_coordinates = [bus_coordinates; {busname, busx, busy}];
    end
end

plot(h3axes, [bus_coordinates{:,2}], [bus_coordinates{:,3}], 'o')
offset = 10;
% text([bus_coordinates{:,2}]+offset, [bus_coordinates{:,3}]+offset,bus_coordinates(:,1),'fontsize',8)
specialnodes_coord = {};
MyColors = distinguishable_colors(length(MyAreas));
for j = 1:length(MyAreas)
    % specialnodes = {'m1142843', 'l3081380', 'l2955077', 'm1166374'};
    specialnodes = MyAreas(j).Buses;
    for i = 1:length(specialnodes)
        node_index = find(strcmpi(bus_coordinates(:,1),specialnodes{i}));
        plot(h3axes, [bus_coordinates{node_index,2}], [bus_coordinates{node_index,3}], 'o', 'MarkerFaceColor',MyColors(j,:), 'MarkerEdgeColor', MyColors(j,:));
        % specialnodes_coord = [specialnodes_coord; {bus_coordinates{node_index,2}, bus_coordinates{node_index,3}}];
    end
end


% Graph Plot - Tree plot
% close all
% close(figure(1))
close(figure(2))
% G = digraph(AdjMatrix,AllBusNames);

% h = figure(1);
% subplot(1,2,1);
% p1 = plot(G,'NodeLabel',G.Nodes.Name);
% hold on;
% highlight(p1,1,'nodecolor','r','MarkerSize',10);
% subplot(1,2,2);
% p2 = plot(G,'NodeLabel',1:length(AllBusNames));
% hold on;
% highlight(p2,1,'nodecolor','r','MarkerSize',10);

h2 = figure(2);
p3 = plot(G);
hold on;
highlight(p3,1,'NodeColor','r','MarkerSize',10);


for j = 1:length(MyAreas)
    % specialnodes = {'m1142843', 'l3081380', 'l2955077', 'm1166374'};
    specialnodes = MyAreas(j).Buses;
    for i = 1:length(specialnodes)
        node_index = find(strcmpi(AllBusNames,specialnodes{i}));
        highlight(p3, node_index, 'NodeColor',MyColors(j,:))
    end

    node_index = find(strcmpi(AllBusNames, MyAreas(j).InterfaceBus));
    highlight(p3,node_index, 'NodeColor', MyColors(j,:), 'MarkerSize', 10);
end




%%

%{
MyAreas = struct;
MyAreas.Descendants = {};
MyAreas.Buses = {};
MyAreas.Index_Buses = {};
MyAreas.InterfaceBus = {};
MyAreas.Index_InterfaceBus = {};
MyAreas.nBuses = {};
MyAreas.Nodes = {};
MyAreas.Index_Nodes = {};
MyAreas.Buses_nPhases = {};

% descendants = [];
% node_name = 'sourcebus';
% 
% node_index = find(strcmpi(AllBusesOrignal,node_name));
% descendants = collectDescendants(descendants, node_index,G);
% 
% descendants_names = AllBusesOrignal(descendants);
% 
% myDescendantstable = table(descendants_names, descendants)

%
% Area1 = setdiff(AllBusesOrignal,Area2);
% Area1 = setdiff(Area1, Area3);
%}

%% Generate my DERs on three phase buses only

myDERs_buses = {};

for i = 1:length(MyAreas)
    for j = 1:MyAreas(i).nBuses
        myDERs_buses = [myDERs_buses; {i, MyAreas(i).Buses{j}}];
    end
end
DERs_Unique_Areas = myDERs_buses(:,1);
DERs_Unique_Areas =  unique([DERs_Unique_Areas{:}]);


%% Consider only m buses for DERs locations
myDERs_buses_m = {};
myDERs_buses_l = {};

for i = 1:size(myDERs_buses,1)
    iDERBus = myDERs_buses(i,2);
    if strcmpi(iDERBus{1}(1), 'm')
        myDERs_buses_m = [myDERs_buses_m; myDERs_buses(i,:)];
    end
    if strcmpi(iDERBus{1}(1), 'l')
        myDERs_buses_l = [myDERs_buses_l; myDERs_buses(i,:)];
    end
end

m_DERs_Unique_Areas = myDERs_buses_m(:,1);
m_DERs_Unique_Areas =  unique([m_DERs_Unique_Areas{:}]);



l_DERs_Unique_Areas = myDERs_buses_l(:,1);
l_DERs_Unique_Areas =  unique([l_DERs_Unique_Areas{:}]);

%%
for i = 1:size(myDERs_buses_l,1)
    node_index = find(strcmpi(bus_coordinates(:,1),myDERs_buses_l{i,2}));
    plot(h3axes, [bus_coordinates{node_index,2}], [bus_coordinates{node_index,3}], '>', 'MarkerFaceColor',[0,0,0] , 'MarkerEdgeColor', [0.2, 0.8, 0.4]);

    node_index = find(strcmpi(AllBusNames,myDERs_buses_l{i,2}));
    highlight(p3, node_index, 'Marker', '>')
end

for i = 1:size(myDERs_buses_m,1)
    node_index = find(strcmpi(bus_coordinates(:,1),myDERs_buses_m{i,2}));
    plot(h3axes, [bus_coordinates{node_index,2}], [bus_coordinates{node_index,3}], '*', 'MarkerFaceColor',[0,0,0] , 'MarkerEdgeColor', [0.2, 0.8, 0.4]);

    node_index = find(strcmpi(AllBusNames,myDERs_buses_m{i,2}));
    highlight(p3, node_index, 'Marker', '*')
end



%% Fomrulatijng my DERS Table and checking number of phases for each bus!
DERsTableData = {};
for i = 1:size(myDERs_buses_m)
    DSSBus = MATDSS.Sim.DSSCircuit.SetActiveBus(myDERs_buses_m{i,2});
    DSSBus = MATDSS.Sim.DSSCircuit.ActiveBus;
    nphase = DSSBus.NumNodes;
    phases = DSSBus.Nodes;
    DERsTableData = [DERsTableData;myDERs_buses_m(i,:), nphase, phases];
end

for i = 1:size(myDERs_buses_l)
    DSSBus = MATDSS.Sim.DSSCircuit.SetActiveBus(myDERs_buses_l{i,2});
    DSSBus = MATDSS.Sim.DSSCircuit.ActiveBus;
    nphase = DSSBus.NumNodes;
    phases = DSSBus.Nodes;
    DERsTableData = [DERsTableData;myDERs_buses_l(i,:), nphase, phases];
end


%% Checking interface buses and their parent bus!

InterfaceBusesTableData = {};
for i = 1:length(InterfaceBuses)
    ibusname = InterfaceBuses{i};
    ibusindex = find(strcmpi(AllBusNames,ibusname));
    iparentindex = getParentNode(G, ibusindex);

    if strcmpi(ibusname, 'sourcebus')
        iparentbus = 'sourcebus';
    else
        iparentbus = AllBusNames{iparentindex};
    end

    InterfaceBusesTableData = [InterfaceBusesTableData; {ibusindex, ibusname, iparentindex, iparentbus}];

end

%% Export Results to Excel
AreasData = repmat("",max([MyAreas(:).nBuses]),1);

AreasName = [];
for i = 1:length(MyAreas)
    AreasData(1:MyAreas(i).nBuses,i) = MyAreas(i).Buses;
    AreasName = [AreasName; convertCharsToStrings(['area ' num2str(i)])];
end

AreasTable = cell2table(cellstr(AreasData));
AreasTable.Properties.VariableNames = AreasName;


% writetable(AreasTable, 'C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500AreasTable.xlsx');
writetable(AreasTable, 'D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500AreasTable.xlsx');




DERsData = myDERs_buses;
DERsTable = cell2table(DERsData);
DERsTable.Properties.VariableNames = {'Area Number', 'm & l buses'};
% writetable(DERsTable, 'C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500DERsTable.xlsx');
writetable(DERsTable, 'D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500DERsTable.xlsx');


DERsTable_ml = cell2table(DERsTableData);
DERsTable_ml.Properties.VariableNames = {'Area number', 'm & l buses', 'nPhase', 'Phases'};
% writetable(DERsTable_ml, 'C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500DERsTable_ml.xlsx')
writetable(DERsTable_ml, 'D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500DERsTable_ml.xlsx')

%% Export Interface buses detials to excel

InterfaceBusesTable = cell2table(InterfaceBusesTableData);
InterfaceBusesTable.Properties.VariableNames = {'Interface Bus Index', 'Interface Bus Name', 'Parent Bus Index', 'Parent Bus Name'};
writetable(InterfaceBusesTable, 'C:\Users\aldah\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500InterfaceBuses.xlsx')
% writetable(InterfaceBusesTable, 'D:\OneDrive - University of Waterloo\PhD UoW\2023\T & D Paper\Simulation Results\MATDSS_ver0.87\Functions\IEEE8500InterfaceBuses.xlsx')

%%
function descendants = collectDescendants(descendants, currentNodeIndex,G)
    % Get the direct successors (children) of the current node.
    children = successors(G, currentNodeIndex);
    
    % Append the children to the descendants list.
    descendants = [descendants; children];
    
    % Recursively call the function for each child.
    for i = 1:numel(children)
        descendants = collectDescendants(descendants, children(i),G);
    end
end


function flattenedArray = flattenCellArray(cellArray)
    flattenedArray = [];
    stack = cellArray; % Initialize a stack with the original cell array
    
    if iscell(stack)
        while ~isempty(stack)
            if size(stack,1) < 2
                disp('hellow');
            end
            element = stack{end};
            stack = stack(1:end-1);

            if length(element) > 1
                flattenedArray= [flattenedArray, element'];
            else
                flattenedArray = [flattenedArray, element];
            end
        end

    else
        flattenedArray = cellArray;
    end


    if size(flattenedArray,1) < size(flattenedArray,2)
        flattenedArray = flattenedArray';
    end
end



function parentNodeIndex = getParentNode(graph, nodeIndex)
    % Get the adjacency matrix of the graph
    adjacencyMatrix = adjacency(graph);
    
    % Check if the node has any parents
    if sum(adjacencyMatrix(:, nodeIndex)) == 0
        parentNodeIndex = []; % No parent node exists
    else
        % Find the index of the parent node
        parentNodeIndex = find(adjacencyMatrix(:, nodeIndex));
    end
end
