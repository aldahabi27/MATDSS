function [MATDSS, DER] = MATDSS_Controllers(app,MATDSS,DER)
% function MATDSS = MATDSS_Controllers(app,MATDSS,DER)
% This function will define the controllers needed following the provided
% information in the Configuration files. Use the configuration option in
% MATDSS Application to define the control areas and assign the buses to
% them. Also specify main controller parameters from within the application
% too.
%
% Make sure you specify the area to include all buses within it. Otherwise,
% the controller will be confused when defining the branches. This
% functionality will be autmated later to ensure the areas are being
% defined as they should be.
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     Parameters needed to define the conrollers    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CA = [];
AllBusNodesNames = {};
AllNodesNames = MATDSS.Sim.Meas.AllNodesNames;
AllBusNames = MATDSS.Sim.Meas.AllBusNames;
for i = 1:length(AllNodesNames)
    iNodeName = AllNodesNames{i};
    iDot = strfind(iNodeName,'.'); % find the location of "dots" in the node name to find the busname
    iBusname = iNodeName(1:iDot(1)-1);
    AllBusNodesNames = [AllBusNodesNames;iBusname];
end
SourceBusName = AllBusNames{1};
% Get all lines/branches names
AllLinesNames = MATDSS.Sim.Meas.AllBranchesPhasesNames;
% Get all lines/branches names
MyLines = MATDSS.Sim.DSSCircuit.Lines;
MyLinesAllNames = MyLines.AllNames;


NumberOfCA = length(unique(str2double(MATDSS.TableData.CASettings(:,3))));


MyPer = 0.55; MyPerEnd = 0.75;
MyPerSec = (MyPerEnd - MyPer)/4;
MATDSSApp_Status(app,MyPer);

% Extracting Control Areas settings from the tables.
for i = 1:NumberOfCA
    Controller = struct;

    [...
        Controller.Area,...
        Controller.Type,...
        Controller.alpha,...
        Controller.rp,...
        Controller.rbard,...
        Controller.E,...
        Controller.vul,...
        Controller.vll,...
        Controller.iul,...
        Controller.arho,...
        Controller.asigma,...
        Controller.alambda,...
        Controller.amu,...
        Controller.aeta,...
        Controller.apsi,...
        Controller.agamma,...
        Controller.anu,...
        Controller.azeta,...
        Controller.crho,...
        Controller.csigma,...
        Controller.clambda,...
        Controller.cmu,...
        Controller.ceta,...
        Controller.cpsi,...
        Controller.cgamma,...
        Controller.cnu,...
        Controller.czeta] = deal(MATDSS.TableData.ControllersSettings{i,:});

    ControllerFields = fieldnames(Controller);
    % Changing strings to numbers in the corresponding locations.
    for j = 1:length(ControllerFields)
        if j ~= [2,3,9]
            Controller.(ControllerFields{j}) = str2double(Controller.(ControllerFields{j}));
        end

        if j == 3 && ~strcmp(Controller.(ControllerFields{j}),'auto')
            Controller.(ControllerFields{j}) = str2double(Controller.(ControllerFields{j}));
        end

    end

    Controller.Buses = cell(0,4); % Check which buses belongs to this controller
    % Adding list of buses that belong to the controller
    for j = 1:size(MATDSS.TableData.CASettings,1)
        if str2double(MATDSS.TableData.CASettings{j,3}) == i
            Controller.Buses = [Controller.Buses; MATDSS.TableData.CASettings(j,1),str2double(MATDSS.TableData.CASettings(j,2)),...
                str2double(MATDSS.TableData.CASettings(j,3)),MATDSS.TableData.CASettings(j,4),...
                MATDSS.TableData.CASettings(j,5)];
        end
    end


    % Check if this area is not L2C, then find the interface bus, and add
    % it to the list of buses.
    if MATDSS_StrComp(Controller.Buses(:,1),SourceBusName) <= 0
        CABus0 = Controller.Buses(find(strcmp('t',Controller.Buses(:,4))),5);
        CASettingsIndex = find(strcmp(MATDSS.TableData.CASettings(:,1),CABus0));
        Controller.Buses = [MATDSS.TableData.CASettings(CASettingsIndex,1),str2double(MATDSS.TableData.CASettings(CASettingsIndex,2)),...
            str2double(MATDSS.TableData.CASettings(CASettingsIndex,3)),MATDSS.TableData.CASettings(CASettingsIndex,4),...
            MATDSS.TableData.CASettings(CASettingsIndex,5);...
            Controller.Buses];
    end


    % Find all indecies of all phases in this area
    CABusesIndices = [];
    ControllerBusesIndices = [];
    for j = 2:size(Controller.Buses,1)
        CABusesIndices = [CABusesIndices; find(strcmp(AllBusNodesNames,Controller.Buses(j,1)))];
        if Controller.Buses{j,3} == i
            ControllerBusesIndices = [ControllerBusesIndices;find(strcmp(AllBusNodesNames,Controller.Buses(j,1)))];
        end
    end

    
    % Save these variables in the controller struct, then save it in CA
    % variable.
    Controller.NodesIndices = CABusesIndices;
    Controller.nBuses = size(Controller.Buses,1);
    Controller.ControllerBusesIndices = ControllerBusesIndices;


    CA = [CA; Controller];


    MATDSSApp_Status(app, MyPer + MyPerSec * i / NumberOfCA);

    clear Controller

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Get All Interface Buses and Their Corresponding Information %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate informationm about virtual DERs and their location (to be used
% mainly in ABMH function). We do that by looking at the branches, where
% the buses belongs to different control areas.

CAInterfaceBuses = {};
% CA#, Bus 1, Bus1CA, Bus 2, Bus2CA, nPhases, Nodes
% 1 CA#:is the Control area where the virtual DER is needed
% 2 Bus1: is where the virtual DER is located.
% 3 Bus1CA: is the control area where Bus1 is located (this is same as CA#)
% 4 Bus2: is the bus connected to Bus1
% 5 Bus2CA: is the Control area connected to the higher level control area. So
% it is the one required to provide the information about P,Q and their
% limits as a VDER.
% 6 nPhases: number of phases at the interface bus. In some cases, this can
% be less than 3 phases interface.
% 7 Nodes: The nodes connected at Bus2 (to Bus1).
%


CABusesTIndex = find(strcmp(MATDSS.TableData.CASettings(:,4),'t'));
CABus0 = {};

% Loop over all interface buses/lines/connections
for i = 1:length(CABusesTIndex)
    % Take bus 1 from teh table
    Bus1 = MATDSS.TableData.CASettings{CABusesTIndex(i),5};
    % Check to which area it belongs
    for j = 1:NumberOfCA
        if MATDSS_StrComp(CA(j).Buses(:,1),Bus1) > 0
            Bus1CA = j; % Get CA number of Bus1
            CANumber = j; % This is the area for which we are finding the connection line.
            break;
        end
    end

    % Get bus 2 from the table, also get its area too
    Bus2 = MATDSS.TableData.CASettings{CABusesTIndex(i),1};
    Bus2CA = MATDSS.TableData.CASettings{CABusesTIndex(i),3};
    if strcmpi(Bus2,SourceBusName) % If here we have it with the source bus, this means we are in L2C, simple!
        nPhases = 3; % Currently we are considering 3phase interface with transmission network. We can later add code to automate this with all possible cases.
        CABus0Phases = {};
        for j = 1:nPhases
            CABus0Phases = [CABus0Phases; [Bus2, '.', num2str(j)]];
        end
        CABus0 = [CABus0; {Bus2CA, MATDSS_StrComp(AllNodesNames,CABus0Phases)}]; % save it then in CABus0
    else % if this is not L2C (we are in lower levels)

        % In this case, we need to do more work. Let's find the line
        % connecting the bus0 and then get all phases information/order.
        for k = 1:length(MyLinesAllNames)
            % Loop over all lines in the circuit, to find the branch that
            % connects our interface bus of LLC to L2C or other LLC.
            MyLines.Name = MyLinesAllNames{k};
            LineBus1 = MyLines.Bus1; % Get bus 1 name
            iDotBus1 = strfind(LineBus1,'.'); % find the location of "dots" in the node name to find the busname
            if isempty(iDotBus1) % If the bus name we got was only the name without nodes, add the nodes manually and then remove them. This is important so our code can work with all cases. 
                LineBus1 = [LineBus1 '.1.2.3'];
                iDotBus1 = strfind(LineBus1,'.'); % find the location of "dots" in the node name to find the busname
            end
            LineBus1FullName = LineBus1; % Keep a copy of the full name for later
            LineBus1 = LineBus1(1:iDotBus1(1)-1); % here is the short name without phases


            % Repeat what we did to Bus1 for Bus2.
            LineBus2 = MyLines.Bus2;
            iDotBus2 = strfind(LineBus2,'.'); % find the location of "dots" in the node name to find the busname
            if isempty(iDotBus2)
                LineBus2 = [LineBus2 '.1.2.3'];
                iDotBus2 = strfind(LineBus2,'.'); % find the location of "dots" in the node name to find the busname
            end
            LineBus2FullName = LineBus2;
            LineBus2 = LineBus2(1:iDotBus2(1)-1);


            % Now check, if one of those buses that we have on this branch
            % is the interface bus that we are interested in.
            if (strcmp(Bus1,LineBus1) && strcmp(Bus2,LineBus2)) || ((strcmp(Bus2,LineBus1) && strcmp(Bus1,LineBus2)))
                nPhases = MyLines.Phases; % Get the number of phases of this CA
                
                % Check if Bus1 is the interface bus we want, else flip
                % them
                if ~strcmp(Bus1,LineBus1)
                    LineBus1FullName = LineBus2FullName;
                    iDotBus1 = iDotBus2;
                end

                CABus0Phases = {};
                % Get phases order
                for j = 1:length(iDotBus1)
                    CABus0Phases = [CABus0Phases; [LineBus1FullName(1:iDotBus1(1)) LineBus1FullName(iDotBus1(j)+1)]];
                end

                % Save Interface bus information in CABus0
                CABus0 = [CABus0; {Bus2CA, MATDSS_StrComp(AllNodesNames,CABus0Phases)}];
                break;
            end
        end % Done with detecting buses and lines that connects interface bus of CA we are considering here

        % generating list of nodes to build our CAInterfaceBuses variable.
        bus1nodes = '[';
        for j = 1:length(iDotBus1)
            bus1nodes = [bus1nodes, LineBus1FullName(iDotBus1(j)+1),','];
        end
        bus1nodes = [bus1nodes(1:end-1),']'];

        CAInterfaceBuses = [CAInterfaceBuses; {CANumber, Bus1, Bus1CA, Bus2, str2num(Bus2CA), nPhases, bus1nodes}];
    end
    MATDSSApp_Status(app, MyPer + MyPerSec*(1 + i / length(CABusesTIndex)));

end % End of our for loop over all interface buses


% Save CA and CAInterfaceBuses variable in MATDSS.Cont
MATDSS.Cont.CA = CA;
MATDSS.Cont.CAInterfaceBuses = CAInterfaceBuses;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                Defing Virtual DERs                %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% So now, we are ready to define our VDERs!
for i = 1:size(CAInterfaceBuses,1)

    VDERsTable = app.VDERsTable.Data.Variables;
    VDERIndexTable = find(CAInterfaceBuses{i,5}==str2double(VDERsTable(:,1)));
    VDERInfo = VDERsTable(VDERIndexTable,[3,4,5,6,7,8,9]);


    [MATDSS, DER] = MATDSS_DERNew(MATDSS,DER,...
        num2str(size(DER,2)+1),...                                                                       DER Index
        ['VDER_CA' num2str(CAInterfaceBuses{i,1}) '_For_CA' num2str(CAInterfaceBuses{i,5})],... DER Name
        strcat("'", CAInterfaceBuses{i,2}, "'"),...                                             DER Bus Location
        VDERInfo(1),...                                                                              Tau
        "VDER",...                                                                              DER Type --> VDER
        CAInterfaceBuses{i,7},...                                                               Nodes for the VDER
        VDERInfo(2),...                                                                               Type of connection (for now will go with Wye, but we might need to have a smarter way of selecting why or delta)
        num2str(CAInterfaceBuses{i,6}),...                                                      Number of phases
        VDERInfo(3),...                                                                          DER Mode
        VDERInfo(4),...                                                                          P(x)
        VDERInfo(5),...                                                                          Q(x)
        '-1e99','1e99','-1e99','1e99',...                                                       Pmin,Pmax,Qmin,Qmax
        VDERInfo(6),VDERInfo(7));                                                                              %ax and cx
    
        MATDSSApp_Status(app, MyPer + MyPerSec*(2 + i / size(CAInterfaceBuses,1)));
end

% Now, we are ready to re-define our control areas (as isolated networks).
% We define our CA-DER list for each CA. We map the location of the DERs
% (DERs and VDERs) to the bus numbers inside the redefined control areas.
%
% Then we will be able to call our ABMH matrix. We define a new measurement
% function that would map the voltages, currents and powers to the new
% isolated networks.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%         GET CA Properties and Information         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:NumberOfCA
    

    % Initialize some properties
    CAnDER = 0; %nDER of CA
    CADERIndex = []; % Corresponding indecies
    DERBuses = [DER(:).BusNum]'; % DER Buses
    %Pmax, Pmin, Qmax, Qmin of CA VDER. Currenlty we sum up all available
    %resources within CA to get an estimate of available energy.
    %
    % This section needs to be rewritten to start from lower areas up to
    % L2C
    %
    % If DER is in the area, add up its limits to our VDER
    Pmax = 0;
    Pmin = 0;
    Qmax = 0;
    Qmin = 0;
    for j = 1:length(DERBuses)
        IsDERInCA = find(DERBuses(j) == [MATDSS.Cont.CA(i).Buses{:,2}]');
        if ~isempty(IsDERInCA)
            if MATDSS.Cont.CA(i).Buses{IsDERInCA,3} == i
                CAnDER = CAnDER + 1;
                CADERIndex = [CADERIndex; j];
                DER(j).CAIndex = i;
                Pmin = Pmin - DER(j).Pmax;
                Pmax = Pmax - DER(j).Pmin;
                Qmin = Qmin - DER(j).Qmax;
                Qmax = Qmax - DER(j).Qmin;
            end
        end
    end
    
    MATDSS.Cont.CA(i).nDER = CAnDER; %Save the number of DERs in CA
    MATDSS.Cont.CA(i).DERIndex = CADERIndex; % Save DERs indecies in CA
    area_index = find (i==str2double([CABus0(:,1)])); % CA index in the list of CABus0
    MATDSS.Cont.CA(i).CABus0 = CABus0{area_index,2}; % Get corresponding Bus0 information
    MATDSS.Cont.CA(i).NodesIndices = [MATDSS.Cont.CA(i).CABus0;MATDSS.Cont.CA(i).NodesIndices]; % All nodes indicies including Bus0
    MATDSS.Cont.CA(i).nNodes = length(MATDSS.Cont.CA(i).NodesIndices); %nNodes = length of all nodes in CA
    MATDSS.Cont.CA(i).nBuses = size(MATDSS.Cont.CA(i).Buses,1); % nBuses = length of Buses in CA
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                    YBus of CA                     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We use Yreduced method to shrink Ybus of the whole network to
    % Yreduced of CA
    FullYbus = MATDSS.Sim.Ybus;
    ControllerYbus = [];
    ControllerYbus = [FullYbus(:,MATDSS.Cont.CA(i).NodesIndices),FullYbus(:,setdiff(1:length(AllBusNodesNames),MATDSS.Cont.CA(i).NodesIndices))];
    ControllerYbus = [ControllerYbus(MATDSS.Cont.CA(i).NodesIndices,:);ControllerYbus(setdiff(1:length(AllBusNodesNames),MATDSS.Cont.CA(i).NodesIndices),:)];

    nphase = length(AllBusNodesNames);
    nCAphase = length(MATDSS.Cont.CA(i).NodesIndices);
    Ynn = ControllerYbus(1:nCAphase,1:nCAphase);
    Yrr = ControllerYbus(nCAphase+1:nphase,nCAphase+1:nphase);
    Ynr = ControllerYbus(1:nCAphase,nCAphase+1:nphase);
    Yrn = ControllerYbus(nCAphase+1:nphase,1:nCAphase);
    ControllerYred = Ynn-Ynr*(Yrr^-1)*Yrn;
    MATDSS.Cont.CA(i).Ybus = ControllerYred;

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      Find P0, Q0, P0Set, Q0Set of CA and VDER     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % If this is L2C (highest control area), then no VDER is needed.
    if strcmp(MATDSS.Cont.CA(i).Buses(1,1),SourceBusName) 
        if ~isfield(MATDSS.Cont,'at')
            MATDSS.Cont.at = MATDSS.Time.Sim.TimeSpan(find(mod(MATDSS.Time.Sim.TimeSpan,MATDSS.Time.Cont.TimeStep)==0));
        end

        MATDSS.Cont.CA(i).P0Set = zeros(1,length(MATDSS.Cont.at));
        MATDSS.Cont.CA(i).Q0Set = zeros(1,length(MATDSS.Cont.at));
        temp = MATDSS.Cont.CA(i).CABus0;
        while ~isempty(temp)
            MATDSS.Cont.CA(i).ControllerBusesIndices(MATDSS.Cont.CA(i).ControllerBusesIndices == temp(1)) = [];
            temp(1) = [];
        end
        MATDSS.Cont.CA(i).VDERIndex = -1; % Indicating that we should use The P0 and Q0 from MATDSS.ControlSignals
    else % Get the VDER index for the CA
        for j = 1:size(CAInterfaceBuses,1) % Record the corresponding VDER index for the CA to retrieve the P and Q Setpoints
            if CAInterfaceBuses{j,5} == MATDSS.Cont.CA(i).Area
                nNotVDER = MATDSS.Sim.nDER - size(CAInterfaceBuses,1);
                MATDSS.Cont.CA(i).VDERIndex = nNotVDER + j;

                if MATDSS.Cont.CA(i).VDERIndex > 0 % This if statement is not needed
                    % Foimd 
                    MATDSS.Sim.DSSCircuit.SetActiveBus(CAInterfaceBuses{j,2});
                    MyBus = MATDSS.Sim.DSSCircuit.ActiveBus;
                    MyBusV = MyBus.Voltages;
                    MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end);
                    MyBusV = MyBusV.';%.*1e3;

                    MyLines = MATDSS.Sim.DSSCircuit.Lines;
                    MyLinesAllNames = MyLines.AllNames;

                    for k = 1:length(MyLinesAllNames)
                        success = MATDSS.Sim.DSSCircuit.SetActiveElement(['Line.' MyLinesAllNames{k}]); % set active element to line k
                        if ~success
                            error('Error in setting active element in Branch current measurements in MATDSS_IPQ.m');
                        end

                        MyLine = MATDSS.Sim.DSSCircuit.ActiveElement; % handle of the current line i
                        MyLineBuses = MyLine.BusNames;
                        
                        iDotBus1 = strfind(MyLineBuses{1},'.'); % find the location of "dots" in the node name to find the busname
                        if ~isempty(iDotBus1) % If the bus name we got was only the name without nodes, add the nodes manually and then remove them. This is important so our code can work with all cases.
                            MyBus1 = MyLineBuses{1};
                            MyLineBuses(1) = {MyBus1(1:iDotBus1-1)};
                        end
                        iDotBus2 = strfind(MyLineBuses{2},'.'); % find the location of "dots" in the node name to find the busname
                        if ~isempty(iDotBus2) % If the bus name we got was only the name without nodes, add the nodes manually and then remove them. This is important so our code can work with all cases.
                            MyBus2 = MyLineBuses{2};
                            MyLineBuses(2) = {MyBus2(1:iDotBus2-1)};
                        end


                        if sum(strcmpi(MyLineBuses,CAInterfaceBuses{j,2})) && sum(strcmpi(MyLineBuses,CAInterfaceBuses{j,4}))
                            MyLineI = MyLine.Currents;
                            MyLineI = MyLineI(1:2:end) + 1i*MyLineI(2:2:end);
                            MyBusV = MyLine.Voltages;
                            MyBusV = MyBusV(1:2:end) + 1i*MyBusV(2:2:end);
                            if strcmpi(MyLineBuses{1},CAInterfaceBuses{j,2})
                                MyLineI = MyLineI(1:length(MyLineI)/2).';
                                MyBusV = MyBusV(1:length(MyBusV)/2).';
                                break;

                            else
                                MyLineI = MyLineI(length(MyLineI)/2 + 1:end).';
                                MyBusV = MyBusV(length(MyBusV)/2 + 1:end).';
                                break;
                            end
%                             if strcmpi(MyLineBuses{1},CAInterfaceBuses{j,2})
%                                 MyLineI = MyLineI(1:length(MyLineI)/2).';
%                                 break;
%                             else
%                                 MyLineI = MyLineI(length(MyLineI)/2 + 1:end).';
%                                 break;
%                             end
                        end
                    end

                    MyCAInterfaceS0 = MyBusV.*conj(MyLineI);
                    MyCAInterfaceP0 = real(MyCAInterfaceS0);
                    MyCAInterfaceQ0 = imag(MyCAInterfaceS0);
                    i_WVar = DER(MATDSS.Cont.CA(i).VDERIndex).WVarIndex;
                    DER(MATDSS.Cont.CA(i).VDERIndex).W = ones(1, length(MATDSS.Cont.at)).*(-sum(MyCAInterfaceP0));
                    DER(MATDSS.Cont.CA(i).VDERIndex).Var  = ones(1, length(MATDSS.Cont.at)).*(-sum(MyCAInterfaceQ0));
                    DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint(:,2) = [ones(length(MATDSS.Cont.at),1).*(-sum(MyCAInterfaceP0))];
                    DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint(:,2) = [ones(length(MATDSS.Cont.at),1).*(-sum(MyCAInterfaceQ0))];
                    MATDSS.Cont.CA(i).P0Set = ones(1,length(MATDSS.Cont.at)).*(-DER(MATDSS.Cont.CA(i).VDERIndex).WSetpoint(1,2));
                    MATDSS.Cont.CA(i).Q0Set = ones(1,length(MATDSS.Cont.at)).*(-DER(MATDSS.Cont.CA(i).VDERIndex).VarSetpoint(1,2));
                    DER(MATDSS.Cont.CA(i).VDERIndex).Pmax = DER(MATDSS.Cont.CA(i).VDERIndex).W(1) + Pmax;
                    DER(MATDSS.Cont.CA(i).VDERIndex).Pmin = DER(MATDSS.Cont.CA(i).VDERIndex).W(1) + Pmin;
                    DER(MATDSS.Cont.CA(i).VDERIndex).Qmax = DER(MATDSS.Cont.CA(i).VDERIndex).Var(1) + Qmax;
                    DER(MATDSS.Cont.CA(i).VDERIndex).Qmin = DER(MATDSS.Cont.CA(i).VDERIndex).Var(1) + Qmin;

%                     DER(MATDSS.Cont.CA(i).VDERIndex).Pmax = DER(MATDSS.Cont.CA(i).VDERIndex).W + 10e3;
%                     DER(MATDSS.Cont.CA(i).VDERIndex).Pmin = DER(MATDSS.Cont.CA(i).VDERIndex).W + -10e3;
%                     DER(MATDSS.Cont.CA(i).VDERIndex).Qmax = DER(MATDSS.Cont.CA(i).VDERIndex).Var + 10e3;
%                     DER(MATDSS.Cont.CA(i).VDERIndex).Qmin = DER(MATDSS.Cont.CA(i).VDERIndex).Var + -10e3;
                end

            end
        end

        %         MATDSS.Cont.CA(i).CABus0 = setdiff(MATDSS.Cont.CA(i).NodesIndices,MATDSS.Cont.CA(i).ControllerBusesIndices);
    end


    %     MATDSS.Cont.CA(i).CABus0 = CABus0;
    [MATDSS.Cont.CA(i).k_v, MATDSS.Cont.CA(i).k_vIndex,~] = intersect(MATDSS.Meas.k_v,MATDSS.Cont.CA(i).ControllerBusesIndices);


    % Find k_i
    k_i = [];
    if MATDSS.Meas.k_i > 0
        k_iLines = AllLinesNames(MATDSS.Meas.k_i);
        for j = 1:length(k_iLines)
            iDotline = strfind(k_iLines{j},'.');
            LineName = k_iLines{j};
            MyLines.Name = LineName(1:iDotline(1)-1); % Get the line name first and handle that line in mylines variable
            LineBus1 = MyLines.Bus1; % Get bus 1 name with all phases
            LineBus2 = MyLines.Bus2; % Get bus 2 name with all phases

            iDotBus1 = strfind(LineBus1,'.');
            iDotBus2 = strfind(LineBus2,'.');
            if ~isempty(iDotBus1)
                LineBus1 = LineBus1(1:iDotBus1(1)-1);
            end
            if ~isempty(iDotBus2)
                LineBus2 = LineBus2(1:iDotBus2(1)-1);
            end
            iBus1 = find(strcmp(MATDSS.Cont.CA(i).Buses(:,1),LineBus1));
            iBus2 = find(strcmp(MATDSS.Cont.CA(i).Buses(:,1),LineBus2));

            if ~isempty(iBus1) && ~isempty(iBus2)
                k_i = [k_i; MATDSS.Meas.k_i(j)];
            end

        end
    end
    MATDSS.Cont.CA(i).k_i = k_i;
    [~,MATDSS.Cont.CA(i).k_iIndex,~] = intersect(MATDSS.Meas.k_i,MATDSS.Cont.CA(i).k_i);


    [A, B, M, H, Mv, MvIndex, Mi,MiBuses] = MATDSS_ABMH(app,MATDSS,DER,MATDSS.Cont.CA(i));
    MATDSS.Cont.CA(i).ABMH.A = A;
    MATDSS.Cont.CA(i).ABMH.B = B;
    MATDSS.Cont.CA(i).ABMH.M = M;
    MATDSS.Cont.CA(i).ABMH.H = H;
    % Initialize the controller dual variables
    MATDSS.Cont.CA(i).Duals.rho = zeros(size(M,1),length(MATDSS.Cont.at)+1);
    MATDSS.Cont.CA(i).Duals.sigma = zeros(size(H,1),length(MATDSS.Cont.at)+1);
    MATDSS.Cont.CA(i).Duals.lambda = zeros(1,length(MATDSS.Cont.at)+1);
    MATDSS.Cont.CA(i).Duals.mu = zeros(1,length(MATDSS.Cont.at)+1);
    MATDSS.Cont.CA(i).Duals.eta = zeros(1,length(MATDSS.Cont.at)+1);
    MATDSS.Cont.CA(i).Duals.psi = zeros(1,length(MATDSS.Cont.at)+1);
    %     if MvIndex == 0
    %         MATDSS.Cont.CA(i).Duals.gamma(:,1) = 0;
    %         MATDSS.Cont.CA(i).Duals.nu(:,1) = 0;
    %     else
    MATDSS.Cont.CA(i).Duals.gamma = zeros(max(size(Mv)),length(MATDSS.Cont.at)+1);
    MATDSS.Cont.CA(i).Duals.nu = zeros(max(size(Mv)),length(MATDSS.Cont.at)+1);
    %     end
    % if isempty(k_i)
        % MATDSS.Cont.CA(i).Duals.zeta(:,1) = 0;
    % else
    MATDSS.Cont.CA(i).Duals.zeta = zeros(max(size(Mi)),length(MATDSS.Cont.at)+1);
    % end


    if strcmpi(MATDSS.Cont.CA(i).alpha,'auto')
        MATDSS.Cont.CA(i).alpha = MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
    else
        MATDSS.Cont.CA(i).alpha = MATDSS.Cont.CA(i).alpha.*MATDSS.Cont.Gain.*MATDSS.Time.Cont.TimeStep;
    end
    MATDSS.Cont.CA(i).ABMH = [];
    %     MATDSS.Cont.CA(i).P0Set = [];
    %     MATDSS.Cont.CA(i).Q0Set = [];

    MATDSS.Cont.CA(i).P0 = zeros(length(MATDSS.Cont.CA(i).CABus0),length(MATDSS.Cont.at));
    MATDSS.Cont.CA(i).Q0 = zeros(length(MATDSS.Cont.CA(i).CABus0),length(MATDSS.Cont.at));
    MATDSSApp_Status(app, MyPer + MyPerSec*(3 + i / NumberOfCA));
end
end
