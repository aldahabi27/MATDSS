function [A, B, M, H, Mv, MvIndex, Mi,MiBuses] = MATDSS_ABMH(app,MATDSS,DER,CA)
%MATDSS_ABMH returns A, B, M and H matrices for MATDSS and connected DERs
% [A, B, M, H] = MATDSS_ABMH(MATDSS,DER)
%
% The function will obtain A, B, M and H matrices following the derivations
% in section 1.4 and 1.50.
% The code here requires the following information to be defined in MATDSS
% and DER structs
%
% > A valid circuit loaded in OpenDSS
% > At least 1 DER device
% > Define Mv (currently it is being set manually)
% > Define Mi (Currently it is being set manually)
%
%
% note that VI (VoltageIndex) means the voltage Indices of all nodes
% except the interface bus (so it starts from node 4 --> 4:end)
%
% In addition, this function will return Mv set and corresponding Indices
% in list of Nodes. Similarly, it returns Mi set and corresponding Indices
% in Mi.


% Defining easy callbacks variables
nNodes = CA.nNodes;
nBuses = CA.nBuses;
nDERs = CA.nDER;
AllNodesNames = MATDSS.Sim.Meas.AllNodesNames;
AllBusNames = MATDSS.Sim.Meas.AllBusNames;
AllLinesNames = MATDSS.Sim.Meas.AllBranchesPhasesNames;
SourceBusName = AllBusNames{1};
CAAllNodesIndices = CA.NodesIndices;
CAAllNodesNames = AllNodesNames(CAAllNodesIndices);
CANodesIndices = CA.ControllerBusesIndices;
CAnNodes = length(CANodesIndices);
YBus = CA.Ybus;
VV = 1:length(CANodesIndices);% VV is voltage index of all nodes except interface bus
VVs = length(VV); % VVs is voltage index size length(VV)
t = MATDSS.Sim.t;
CANodesNames = AllNodesNames(CANodesIndices);
CAAllBusNames = CA.Buses(:,1);
k_v = CA.k_v;
k_i = CA.k_i;
%% QLv (unitless Matrix)
% QLv is the matrix we use to select/order the voltages of the phases
% (nodes) as in Mv array. This is used mainly when we refer to the
% linearized voltage v of the phaes/buses that we want to ensure voltage is
% within the allowed limits (for which our DERs would respond to make them
% within the limits if they are not).



%% QLv (unitless Matrix)
% QLv is the matrix we use to select/order the voltages of the phases
% (nodes) as in Mv array. This is used mainly when we refer to the
% linearized voltage v of the phaes/buses that we want to ensure voltage is
% within the allowed limits (for which our DERs would respond to make them
% within the limits if they are not).


if k_v >= 0
    Mv = AllNodesNames(k_v); % Vector of boolean flags for which nodes to extract from all voltages (nodes not buses)
    % To get Mv Indices in all nodes names
    MvIndex = MATDSS_StrComp(CANodesNames,Mv);
    
    
    
    QLv = zeros(size(Mv,1),VVs); %VVs is size of VV (nNodes-3)
    
    % the way to build QLv is by looping at each Mv element, and then
    % accordinlgy, map its bus and phase to set 1 to the corresponding node in
    % the nodeslist. So if Mv(1) is for bus 6 phase 2 (let's say this is the
    % element #18 in the nodes list), then the QLv would have 1 at the 18th
    % element and then 0 elsewhere in that row.
    
    
    if sum(MvIndex < 1)
        error('Error in QLv calculations in MATDSS_ABMH');
    else
        for i = 1:length(MvIndex)
            QLv(i,MvIndex(i)) = 1;
        end
    end
else
    Mv = 0;
    MvIndex = 0;
    QLv = 0;
end
% With this, we have QLv ready

% debugging:
% If your buses/phases are inorder as in "AllNodesNames", then your QLv
% should have a diagonal "1" entries for the corresponding buses/phases.
% Always check this back when you select all phases to be controlled, and
% then look at the corresponding QLv matrix. (it should be the identity of
% the size = length(AllNodesNames)-3.



%% QR (Unitless)
% QR is unitless because it maps the P and Q of the DER to the
% corresponding location in the buses/phases order. The units if MWye,
% MDelta, KWye, KDelta, GWye, GDelta, JWye and JDelta would map the P and Q
% correctly to V, I and VA. QR just maps the location of P and Q to the
% correpsonding bus/phase rows.
%
%


% We first initiallize QR and set its size!
if nDERs >= 1
    QR = zeros(4*(VVs),2*nDERs);
    for i = 1:nDERs % we will loop over each DER device to set its corresponding values (columns) in QR
        busname = DER(CA.DERIndex(i)).BusName; % check to which bus is DER i connected
        busnum = MATDSS.Sim.DSSCircuit.SetActiveBus(busname); %get bus number
        mybus = MATDSS.Sim.DSSCircuit.Buses(busname); %handle the bus with mybus variable
        if MATDSS_StrComp(DER(CA.DERIndex(i)).ConnType, 'Wye') %check the connection, if wye
            for j = 1:DER(CA.DERIndex(i)).NPhases %loop over the phases now (nodes)
                nodename = [busname '.' num2str(DER(CA.DERIndex(i)).Nodes(j))]; % check which phase we have included
                k = MATDSS_StrComp(CANodesNames,nodename); % get the corresponding index in AllNodesNames
                if sum( k < 1) %confirm that we don't have any errors in search
                    error('Error in QR calculations in MATDSS_ABMH')
                end
                QR(k,(i-1)*2 + 1) = DER(CA.DERIndex(i)).MultiPhaseRatio; %set the correspoinding ratio/value at that node for that DER (P)
                QR(k + CAnNodes,(i-1)*2 + 2) = DER(CA.DERIndex(i)).MultiPhaseRatio; %set the correspoinding ratio/value at that node for that DER (Q)
            end
        else % if delta connection
            if nchoosek(length(DER(CA.DERIndex(i)).Nodes),2) >= DER(CA.DERIndex(i)).NPhases % Confirm we have complete set of phases defined correclty for that DER
                for j = 1:DER(CA.DERIndex(i)).NPhases % loop over the phases
                    nodename = [busname '.' num2str(DER(CA.DERIndex(i)).Nodes(j))];
                    k = MATDSS_StrComp(CANodesNames,nodename); % Get the corresponding phase index in AllNodesNames
                    if sum( k < 1 )
                        error('Error in QR calculations in MATDSS_ABMH')
                    end
                    QR(k + 2*(CAnNodes),(i-1)*2 + 1) = DER(CA.DERIndex(i)).MultiPhaseRatio; %set the correspoinding ratio/value at that node for that DER (P)
                    QR(k + 3*(CAnNodes),(i-1)*2 + 2) = DER(CA.DERIndex(i)).MultiPhaseRatio; %set the correspoinding ratio/value at that node for that DER (Q)
                end
            else
                msgbox('Error in QR in MATDSS_ABMH')
            end
        end
    end

else
    QR = 0;
end


% QR = -QR;
%% Obtaining MWye and MDelta
%{
to get MWye and MDelta, we need the following:
1. YLL: Ybus matrix where the connection to the inerface bus is ignored. So,
we need to remove the first three rows and cols from Ybus matrix to get
this YLL

2. vhat: This is the load flow results (voltage profile). We just need to
also remove hte interface bus from the vector.

3. G matrix: This is the matrix of Gammas. We need to pay attention to the
number of phases connected at each bus. Also, pay attention to the phase
orders in the vhat vector!
%}


%% YBus -> Y00, Y0L, YL0 and YLL (unit is siemens)
Y00 = YBus(1:length(CA.CABus0),1:length(CA.CABus0));
YL0 = YBus(length(CA.CABus0)+1:end,1:length(CA.CABus0));
Y0L = YBus(1:length(CA.CABus0),length(CA.CABus0)+1:end);
YLL = YBus(length(CA.CABus0)+1:end,length(CA.CABus0)+1:end);


%% vhat --> This is phase to ground voltage (Unit in Volts)
% MATDSS.Sim.DSS_Solution.Solve;
% MATDSS = MATDSS_VProfile(MATDSS);
vhat = MATDSS.Sim.Meas.VProfile(CAAllNodesIndices,end); % Get all voltages in vhat
v0 = vhat(1:length(CA.CABus0)); % save v0 from vhat
vhat(1:length(CA.CABus0)) = []; % remove interface bus volts from vhat


%% G Matrix (This matrix is unitless, used to convert line to line to line to ground)
gamma = [1, -1, 0; 0, 1, -1; -1, 0, 1];

% to get G, we first define it as zeros, then will loop through each bus to
% obtain the gamma equivalent components

G = zeros(VVs);
BusNodes = cell(nBuses,2);
for i = 1:nNodes
    iNodeName = CAAllNodesNames{i};
    iDot = strfind(iNodeName,'.'); % find the location of "dots" in the node name to find the busname
    iBusname = iNodeName(1:iDot(1)-1);
    iBusnum = MATDSS_StrComp(CAAllBusNames,iBusname); % Get the bus number in the bus names list
    BusNodes{iBusnum,2} = [BusNodes{iBusnum,2}, str2double(iNodeName(iDot+1:end))];
    BusNodes{iBusnum,1} = iBusnum;
end


G_matrix = [];
% Generating Gamma for each bus excluding interface bus
for i = 2:size(BusNodes,1)
    igamma = gamma(BusNodes{i,2},BusNodes{i,2});
    G_matrix = blkdiag(G_matrix,igamma);
end


%% now we are ready to compute MWye and MDelta   (Unit is 1/(S*V) = 1/A ---> 1/Current)
% The Mwye and MDelta links our DER P and Q to phases voltages. Check
% equation 1.5 in T & D Coordination Report (Iteration #3)
%
% v = MWye xWye + MDelta xDelta + a
    

MWye = [YLL\(diag(conj(vhat))\eye(size(diag(conj(vhat))))), -1i*(YLL\eye(size(YLL)))*(diag(conj(vhat))\eye(size(diag(conj(vhat)))))];
MDelta = [YLL\ctranspose(G_matrix)*(diag(G_matrix*conj(vhat))\eye(size(diag(G_matrix*conj(vhat))))), -1i*(YLL\eye(size(YLL)))*ctranspose(G_matrix)*(diag(G_matrix*conj(vhat))\eye(size(diag(G_matrix*conj(vhat)))))];


% Below is temp, delete later;

% MWye = [inv(YLL)*inv(diag(conj(vhat))), -1i*inv(YLL)*inv(diag(conj(vhat)))];
% MDelta = [inv(YLL)*conj(G_matrix)*inv(diag(G_matrix*conj(vhat))), -1i*inv(YLL)*conj(G_matrix)*inv(diag(G_matrix*conj(vhat)))];
%
%% to get w, we run CalcVoltageBases --> This is phase to ground voltage (Unit in Volts)
% Note we have "w" which is the vector and we have "W" which is diag(w)
MATDSS.Sim.DSS_Text.Command = 'CalcVoltageBases';
w = MATDSS.Sim.DSSCircuit.AllBusVolts; %Phase to ground
w = w(1:2:end) + 1i*w(2:2:end);
w = w(CAAllNodesIndices);
w = reshape(w,nNodes,[]);
w = w(length(CA.CABus0)+1:end); %exclude interface bus voltage


%% Now, we have MWye, MDelta and w, Let's get KWye and KDelta (Same units as MWye and MDelta, 1/S*V = 1/A)
W = diag(w);
KWye = abs(W)*real(W\MWye);
KDelta = abs(W)*real(W\MDelta);



%% Now, we can get GWye and GDelta for interface bus equations (M and H) (Unitless, V*S/A = V*S/V*S = 1)
GWye = diag(v0)*conj(Y0L)*conj(MWye);
GDelta = diag(v0)*conj(Y0L)*conj(MDelta);

GBar = [GWye GDelta]*QR;

%% Lastly, we need to find B matrix.
% JWye and JDelta units is S/Q1V*S = 1/V.
%

% We need to figure out the lines connections in the circuit.

% set of branch currents to watch, we define it by an array of bus
% Indices. A for loop after that will confirm that these lines do exist
% The code developed for obtaining B matrix is written differently as the
% current is concerned with two buses, and the equation includes more
% variables (E, Z, Y) and JWye and JDelta. Therefore, look at this code
% with the derivations at hand to match the varibales and understand the
% logic here.

% This is the main handle that will be used to extract branch details
myLines = MATDSS.Sim.DSSCircuit.Lines;

if sum(k_i <= 0) || isempty(k_i)
    Mi = [];
    MiBuses = 0;
    B = [zeros(1,size(QR,2))];
else

    Mi = AllLinesNames(k_i);
    % to build the index for Mi
    if min(size(Mi)) < 1 % No current tracking
        B = [zeros(1,size(QR,2))];
    else
        SourceBusIndex = []; %Indices of lines connected to the source bus.
        % for i = 1:size(Mi,1)
        %     iDotline = strfind(Mi{i},'.');
        %     LineName = Mi{i};
        %     myLines.Name = LineName(1:iDotline(1)-1); % Get the line name first and handle that line in mylines variable
        %     b1 = myLines.Bus1;
        %     b2 = myLines.Bus2;
        %     if MATDSS_StrComp(lower(b1),'sourcebus') || MATDSS_StrComp(lower(b2),'sourcebus')
        %         SourceBusIndex = [SourceBusIndex; i];
        %         % break;
        %     end
        % end
        % 
    
        Mi(SourceBusIndex) = [];

        MiBuses = cell(length(Mi),3); % Define a general cell matrix to save all variables in, to be handy for debugging.
        Bbar = {}; % This is Bbar which we will have for each branch
        B = []; % This is the B matrix thaty we will build.
    
        for i = 1:size(Mi,1)
            iDotline = strfind(Mi{i},'.');
            LineName = Mi{i};
            myLines.Name = LineName(1:iDotline(1)-1); % Get the line name first and handle that line in mylines variable
            %     {myLines.Bus1, myLines.Bus2}
            MiBuses{i,1} = Mi{i}; % save the name in MiBuses
            MiBuses{i,2} = myLines.Bus1; % Get bus 1 name with all phases
            MiBuses{i,3} = myLines.Bus2; % Get bus 2 name with all phases
    
    
            iDotBus1 = strfind(myLines.Bus1,'.'); % find the location of "dots" in the node name to find the busname
            % if iDotBus1 is empty, this means the connection is three phase!
            % the code below handles the "mentioned" phases correctly. We will
            % develop a special code here to handle the case where phases are not
            % present. We will simply add them to the bus name! and then continue
            % using the same code
            if isempty(iDotBus1)
                MiBuses{i,2} = [MiBuses{i,2}, '.1.2.3'];
                iBus1name = MiBuses{i,2};
                iDotBus1 = strfind(iBus1name,'.'); % find the location of "dots" in the node name to find the busname
                iBus1name = iBus1name(1:iDotBus1(1)-1); %get bus 1 name alone to handle it
            else
                iBus1name = myLines.Bus1(1:iDotBus1(1)-1); %get 2 us 2 name
            end
            iBus1num = find(contains(AllBusNames,iBus1name)); % Get the bus number in the bus names list
            MiBuses{i,4} = iBus1num; %save bus 1 num
    
    
    
            iDotBus2 = strfind(myLines.Bus2,'.'); % find the location of "dots" in the node name to find the busname
            % WE repreat the same process here to fix the bus name
            if isempty(iDotBus2)
                MiBuses{i,3} = [MiBuses{i,3}, '.1.2.3'];
                iBus2name = MiBuses{i,3};
                iDotBus2 = strfind(iBus2name,'.'); % find the location of "dots" in the node name to find the busname
                iBus2name = iBus2name(1:iDotBus2(1)-1); %get bus 1 name alone to handle it
            else
                iBus2name = myLines.Bus2(1:iDotBus2(1)-1); %get 2 us 2 name
            end
            iBus2num = find(contains(AllBusNames,iBus2name)); % Get the bus number in the bus names list
            MiBuses{i,5} = iBus2num; % save bus 2 num
    
            %calculating zij
            %%%%%%%% consider revisiting - check yprim in cktelement %%%%%%%%
            n = myLines.Phases; % number of phases in this line
            r = myLines.Rmatrix; r = reshape(r,[],n).'; % in Ohms/unitlength             %%%%%%%%%%%%%%%%%%%%%% Multiply with the line length (r, x and y)
            x = myLines.Xmatrix; x = reshape(x,[],n).'; % in Ohms/unitlength
            zij = r + 1i.*x; %zij matrix
            MiBuses{i,6} = zij; %save it
            %     MiIndex = [MiIndex; [MiIndexRef(i).*ones(n,1) + [0;1;2].*[ones(n,1);zeros(3-n,1)]]];
    
            % calculating yij
            omega = 2*pi*60; %f = 60 Hz
            c = myLines.Cmatrix; % in nF/unitlength
            c = c*1e-9; % F/unitlength
            c = reshape(c,[],n); %get it in 3x3 shape (or 2x2 or 1x1)
            ysij = 1i.*omega.*c; %S/unitlength
            MiBuses{i,7} = ysij; % get yc
    
    
            %     if n < 3
            %         disp('hi');
            %     end
            % calculating Ei and Ej
            Ei = zeros(n,VVs); %initiallize Ei
            for j = 1:n % loop over the phases
                Ei_eye_index = MATDSS_StrComp(CAAllNodesNames,[MiBuses{i,2}(1:iDotBus1(1)), MiBuses{i,2}(iDotBus1(j)+1)]); % Get the bus number in the bus names list
                if Ei_eye_index > 3 % Ei_eye_index is not a source bus node
                    %         iphase = str2double(MiBuses{i,2}(iDotBus1(j)+1)); %index of the corresponding phase
                    Ei(j,Ei_eye_index-3) = 1; % set Ei to eye(n) at the correct bus
                end
            end
            %     Ei = [zeros(n,n*(iBus1num-1)),eye(n),zeros(n,n*(nBuses-iBus1num-1))];
            %     Ej = [zeros(n,n*(iBusnum2-1)),eye(n),zeros(n,n*(nBuses-iBusnum2-1))];
    
    
            Ej = zeros(n,VVs); %initiallize Ej
            for j = 1:n % loop over the phases
                Ej_eye_index = MATDSS_StrComp(CAAllNodesNames,[MiBuses{i,3}(1:iDotBus2(1)), MiBuses{i,3}(iDotBus2(j)+1)]); % Get the bus number in the bus names list
    
                if Ej_eye_index > 3 % Ej_eye_index is not a source bus node
                    %         iphase = str2double(MiBuses{i,3}(iDotBus2(j)+1)); %index of the corresponding phase
                    Ej(j,Ej_eye_index-3) = 1; % set Ei to eye(n) at the correct bus
                end
            end
    
            %     Ej = zeros(n,VVs); %repeat Ei process for Ej
            %     Ej_eye_index = find(contains(AllNodesNames,myLines.Bus2(1:iDotBus2(2)-1))); % Get the bus number in the bus names list
            %     Ej(:,Ej_eye_index-3) = eye(n,length(Ej_eye_index));
    
            MiBuses{i,8} = Ei; % now save both Ei and Ej
            MiBuses{i,9} = Ej;
            JWyeij = ((ysij + zij\eye(size(zij)))*Ei - (zij\eye(size(zij)))*Ej)*MWye; %calculate JWyeij
            JDeltaij = ((ysij + zij\eye(size(zij)))*Ei - (zij\eye(size(zij)))*Ej)*MDelta; %calculate JDeltaij
            JWyeij = abs(JWyeij);
            JDeltaij = abs(JDeltaij);
            MiBuses{i,10} = JWyeij; %save both
            MiBuses{i,11} = JDeltaij;
            BbarLine = [JWyeij JDeltaij]*QR; %Get Bbar_i
            Bbar{i} = BbarLine(str2double(LineName(iDotline+1)),:);
            B = [B;Bbar{i}]; %Save it directly in B also!
            
            %     clear Ei Ej c x r n zij ysij
        end
    end
end
% With this, we have the B matrix ready



%% Obtaining A, B, M and H now

if QLv == 0
    A = [zeros(1,size(QR,2))];
else
    A = QLv*[KWye KDelta]*QR;  % Unit is [unitless] * [1/A] * [unitless]
end
% B = B;  % Unit is [1/V] * [unitless] = 1/V
M = real(GBar);  % [unitless] * [unitless] = unitless
H = imag(GBar);  % [unitless] * [unitless] = unitless
