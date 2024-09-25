function MATDSS = MATDSS_VProfile(MATDSS)
% MATDSS_VProfile(MATDSS)
% This function evaluates and updates the voltage profile of the distribution
% system being simulated in MATDSS. It extracts the voltage magnitudes, per-unit
% voltages, and complex voltages (magnitude and angle) from OpenDSS for each bus
% and node.
%
% The voltage profiles are stored in both MATDSS.Meas and MATDSS.Sim.Meas, where:
%   - MATDSS.Meas stores the profiles at time steps when measurements need to be
%     updated for controllers and other modules that rely on specific sampling times.
%   - MATDSS.Sim.Meas stores all voltage profiles at every time step, which is
%     helpful for debugging and visualizing the simulation's progression.
%
% At the first time step (t=1), the function saves the bus names and voltage
% bases. It also calculates voltage phase indices and ensures that the voltage
% profiles align with the node list. The function supports both real-time
% simulations and loaded simulations (where pre-saved data is used).
%
% Parameters:
%   - MATDSS: The MATDSS application instance containing properties related to
%             the simulation and measurement configurations.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Define t and at for easier reference to current simulation time and absolute time
t = MATDSS.Sim.t;
at = MATDSS.Sim.at;

% At the first time step, save bus and node names from OpenDSS
if t == 1
    AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames;   % Get all bus names
    AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames; % Get all node names
else
    % Use previously stored bus and node names in subsequent steps
    AllBusNames = MATDSS.Sim.Meas.AllBusNames;
    AllNodesNames = MATDSS.Sim.Meas.AllNodesNames;
end

% Get voltage measurements from OpenDSS
AllBusVmag = MATDSS.Sim.DSSCircuit.AllBusVmag.';   % Voltage magnitudes of all buses
AllBusVmagPu = MATDSS.Sim.DSSCircuit.AllBusVmagPu.'; % Per-unit voltage magnitudes
v = MATDSS.Sim.DSSCircuit.AllBusVolts;             % Complex voltages (real + imag)
v = v(1:2:end) + 1i * v(2:2:end);                  % Reshape to get complex voltage
v = reshape(v, MATDSS.Sim.DSSCircuit.NumNodes, []); % Reshape voltage array to node count

% Check if it's time to update measurements based on the defined time step for measurements
if t == 1 || mod(at, MATDSS.Time.Meas.TimeStep) == 0
    % Update measurements at the current absolute time 'at'

    % Calculate voltage bases at the first time step if not using preloaded simulation data
    if t == 1
        if ~MATDSS.UseSimDataFlag
            myBusVBases = [];                      % Initialize array for bus voltage bases
            myVBases = zeros(size(AllNodesNames)); % Initialize zero array for node voltage bases

            % Loop through each bus to get its voltage base (in kV)
            for i = 1:length(AllBusNames)
                MATDSS.Sim.DSSCircuit.SetActiveBus(AllBusNames{i});
                mybus = MATDSS.Sim.DSSCircuit.ActiveBus;
                myBusVBases = [myBusVBases; mybus.kVBase];
            end

            % Map bus voltage bases to the corresponding nodes
            AllBusNodesNames = cellfun(@(x) x(1:end-2), AllNodesNames, 'UniformOutput', false);
            myIndex = MATDSS_StrComp(AllBusNames, AllBusNodesNames);
            myVBases = myBusVBases(myIndex);
            MATDSS.Meas.VBases = myVBases;

            % Determine voltage phase indices for tracking, using 'MATDSS_VI.xlsx'
            v_size = -1;
            for i = 1:size(MATDSS.TableData.VI, 1)
                if MATDSS.TableData.VI(i, 1) == "" || ismissing(MATDSS.TableData.VI(i, 1))
                    break;
                end
                v_size = i;
            end

            if v_size <= 0
                k_v = 0; % No valid phases found
            else
                A = AllNodesNames;                     % All node names from OpenDSS
                B = cellstr(MATDSS.TableData.VI(1:v_size, 1)); % Valid phases from 'MATDSS_VI.xlsx'
                [C, iA, iB] = intersect(A, B, 'stable');  % Find matching node indices
                k_v = iA;
            end
            MATDSS.Meas.k_v = k_v; % Store the phase indices

            % Store bus and node names in the simulation measurements for future reference
            MATDSS.Sim.Meas.AllNodesNames = AllNodesNames;
            MATDSS.Sim.Meas.AllBusNames = AllBusNames;
        else
            % Use preloaded simulation data for voltage base and phase index
            MATDSS.Meas.k_v = MATDSS.LoadedSimData.MATDSS.Meas.k_v;
            MATDSS.Meas.VBases = MATDSS.LoadedSimData.MATDSS.Meas.VBases;
            MATDSS.Sim.Meas.AllNodesNames = MATDSS.LoadedSimData.MATDSS.Sim.Meas.AllNodesNames;
            MATDSS.Sim.Meas.AllBusNames = MATDSS.LoadedSimData.MATDSS.Sim.Meas.AllBusNames;
        end

        % Initialize matrices to store voltage profiles for later time steps
        k_v = MATDSS.Meas.k_v;
        MATDSS.Meas.VProfile = zeros(length(v(k_v)), length(MATDSS.Meas.at));
        MATDSS.Meas.VProfileMag = zeros(length(abs(v(k_v))), length(MATDSS.Meas.at));
        MATDSS.Meas.VProfileAng = zeros(length(angle(v(k_v))), length(MATDSS.Meas.at));

        MATDSS.Meas.VMagProfile = zeros(length(AllBusVmag(k_v)), length(MATDSS.Meas.at));
        MATDSS.Meas.VMagProfilePu = zeros(length(AllBusVmagPu(k_v)), length(MATDSS.Meas.at));

        % Initialize simulation measurements for all time steps
        MATDSS.Sim.Meas.VMagProfile = zeros(length(AllBusVmag), length(MATDSS.Sim.Meas.at));
        MATDSS.Sim.Meas.VMagProfilePu = zeros(length(AllBusVmagPu), length(MATDSS.Sim.Meas.at));

        MATDSS.Sim.Meas.VProfile = zeros(length(v), length(MATDSS.Sim.Meas.at));
        MATDSS.Sim.Meas.VProfileMag = zeros(length(abs(v)), length(MATDSS.Sim.Meas.at));
        MATDSS.Sim.Meas.VProfileAng = zeros(length(angle(v) .* 180 / pi), length(MATDSS.Sim.Meas.at));
    end

    % Update the voltage profile at the current absolute time index 'at_index'
    at_index = find(MATDSS.Meas.at == at);

    % Save only required voltage-phases (from 'MATDSS_VI.xlsx')
    k_v = MATDSS.Meas.k_v;
    if k_v > 0
        % Store voltage profile, magnitude, and angle
        MATDSS.Meas.VProfile(:, at_index) = [v(k_v)];
        MATDSS.Meas.VProfileMag(:, at_index) = [abs(v(k_v))];
        MATDSS.Meas.VProfileAng(:, at_index) = [angle(v(k_v)) .* 180 / pi];

        % Store OpenDSS voltage magnitudes and per-unit values
        MATDSS.Meas.VMagProfile(:, at_index) = [AllBusVmag(k_v)];
        MATDSS.Meas.VMagProfilePu(:, at_index) = [AllBusVmagPu(k_v)];
    end
end

% Always save voltage measurements for debugging and plotting (MATDSS.Sim.Meas)
MATDSS.Sim.Meas.VMagProfile(:, t) = AllBusVmag;
MATDSS.Sim.Meas.VMagProfilePu(:, t) = AllBusVmagPu;

MATDSS.Sim.Meas.VProfile(:, t) = v;
MATDSS.Sim.Meas.VProfileMag(:, t) = abs(v);
MATDSS.Sim.Meas.VProfileAng(:, t) = angle(v) .* 180 / pi;

end


%% Old function code

%{

function MATDSS = MATDSS_VProfile(MATDSS)
% MATDSS = MATDSS_VProfile(MATDSS)
% This function will evaluate the voltage profile of the circuit and store
% the values in MATDSS.Meas and MATDSS.Sim.Meas.
%
% The function will update the measurements in MATDSS.Meas if the current
% time step is a point where measurements needs to be updated following the
% time settings provided for measurement updates.
%
% The function will always update the meausrements in MATDSS.Sim.Meas,
% which would be handy for debugging and plotting later to show the actual
% steps of values, while observing the response of the controller(s).
%
%
%


% define t and at to easily refer to them below
t = MATDSS.Sim.t;
at = MATDSS.Sim.at;


if t == 1 % at t = 1, save multiple values (bus names...etc)
    AllBusNames = MATDSS.Sim.DSSCircuit.AllBusNames; % get all bus names
    AllNodesNames = MATDSS.Sim.DSSCircuit.AllNodeNames; % get all nodes names
else
    AllBusNames = MATDSS.Sim.Meas.AllBusNames;
    AllNodesNames = MATDSS.Sim.Meas.AllNodesNames;
end


% Get new Voltage Profile measurements
AllBusVmag = MATDSS.Sim.DSSCircuit.AllBusVmag.';
AllBusVmagPu = MATDSS.Sim.DSSCircuit.AllBusVmagPu.';
v = MATDSS.Sim.DSSCircuit.AllBusVolts;
v = v(1:2:end) + 1i*v(2:2:end);
v = reshape(v,MATDSS.Sim.DSSCircuit.NumNodes,[]);


% Check if it is time to update the measurements for our controllers.
if t == 1 || mod(at,MATDSS.Time.Meas.TimeStep) == 0  % time to update the measurements
    % MATDSS.Meas.at = [MATDSS.Meas.at, at];

    % Calculate voltage bases (only when t = 1
    if t == 1
        if ~MATDSS.UseSimDataFlag

            % Old Vbase calculations, replaced with the one below
            %         MATDSS.Meas.VBases = MATDSS.Sim.DSS_Circuit.AllBusVmagPu.*MATDSS.Sim.DSS_Circuit.AllBusVmag; % Get all buses bases and save them for later.

            myBusVBases = []; % build an array of all bus voltage bases
            myVBases = zeros(size(AllNodesNames)); % build a zero array of all bases that will be used later
            % for voltage calculations. We always want to have a mapping
            % with the list of nodes not buses.

            % we build our bus Vbases array first by looping over the buses
            for i = 1:length(AllBusNames)
                MATDSS.Sim.DSSCircuit.SetActiveBus(AllBusNames{i});
                mybus = MATDSS.Sim.DSSCircuit.ActiveBus;
                myBusVBases = [myBusVBases;mybus.kVBase];
            end

            % We get the list of buses from the nodes list (bus where each node
            % belongs to!
            AllBusNodesNames = cellfun(@(x) x(1:end-2),AllNodesNames,'UniformOutput',false);
            % Get the index of the corresponding bus in the bus names list
            myIndex = MATDSS_StrComp(AllBusNames,AllBusNodesNames);
            % get the array of Vbases mapping to each node in our list.
            myVBases = myBusVBases(myIndex);
            MATDSS.Meas.VBases = myVBases;


            % Getting Voltage phases indecices for tracking
            v_size = -1;
            for i = 1:size(MATDSS.TableData.VI,1)
                if MATDSS.TableData.VI(i,1) == "" || ismissing(MATDSS.TableData.VI(i,1))
                    break;
                end
                v_size = i;
            end
            if v_size <= 0
                k_v = 0;
            else
                A = AllNodesNames;
                B = cellstr(MATDSS.TableData.VI(1:v_size,1));
                % k_v = MATDSS_StrComp(AllNodesNames,cellstr(MATDSS.TableData.VI(1:v_size,1)));
                [C, iA, iB] = intersect(A,B, 'stable');
                k_v = iA;
            end
            MATDSS.Meas.k_v = k_v;

            MATDSS.Sim.Meas.AllNodesNames = AllNodesNames;
            MATDSS.Sim.Meas.AllBusNames = AllBusNames;
        else
            MATDSS.Meas.k_v = MATDSS.LoadedSimData.MATDSS.Meas.k_v;
            MATDSS.Meas.VBases = MATDSS.LoadedSimData.MATDSS.Meas.VBases;
            MATDSS.Sim.Meas.AllNodesNames = MATDSS.LoadedSimData.MATDSS.Sim.Meas.AllNodesNames;
            MATDSS.Sim.Meas.AllBusNames = MATDSS.LoadedSimData.MATDSS.Sim.Meas.AllBusNames;
        end
        k_v = MATDSS.Meas.k_v;
        MATDSS.Meas.VProfile = zeros(length(v(k_v)),length(MATDSS.Meas.at));
        MATDSS.Meas.VProfileMag = zeros(length(abs(v(k_v))), length(MATDSS.Meas.at));
        MATDSS.Meas.VProfileAng = zeros(length(angle(v(k_v))),length(MATDSS.Meas.at));

        % Use OpenDSS measurements and calculations
        MATDSS.Meas.VMagProfile = zeros(length(AllBusVmag(k_v)),length(MATDSS.Meas.at));
        MATDSS.Meas.VMagProfilePu = zeros(length(AllBusVmagPu(k_v)),length(MATDSS.Meas.at));

        %         if size(MATDSS.Meas.VBases,2) > size(MATDSS.Meas.VBases,1)
        %             MATDSS.Meas.VBases = MATDSS.Meas.VBases';
        %         end



        MATDSS.Sim.Meas.VMagProfile = zeros(length(AllBusVmag), length(MATDSS.Sim.Meas.at));
        MATDSS.Sim.Meas.VMagProfilePu =  zeros(length(AllBusVmagPu), length(MATDSS.Sim.Meas.at));

        MATDSS.Sim.Meas.VProfile =  zeros(length(v), length(MATDSS.Sim.Meas.at));
        MATDSS.Sim.Meas.VProfileMag =  zeros(length(abs(v)), length(MATDSS.Sim.Meas.at));
        MATDSS.Sim.Meas.VProfileAng =  zeros(length(angle(v).*180/pi), length(MATDSS.Sim.Meas.at));

    end

    at_index = find(MATDSS.Meas.at == at);
    % Save only required voltage-phases (from 'MATDSS_VI.xlsx')
    k_v = MATDSS.Meas.k_v;

    if k_v > 0
        % Manually calculate the voltages mag. and angles
        MATDSS.Meas.VProfile(:,at_index) = [v(k_v)];
        MATDSS.Meas.VProfileMag(:,at_index) = [abs(v(k_v))];
        MATDSS.Meas.VProfileAng(:,at_index) = [angle(v(k_v)).*180/pi];

        % Use OpenDSS measurements and calculations
        MATDSS.Meas.VMagProfile(:,at_index) = [AllBusVmag(k_v)];
        MATDSS.Meas.VMagProfilePu(:,at_index) = [AllBusVmagPu(k_v)];
    else
        % MATDSS.Meas.VProfile = [MATDSS.Meas.VProfile, 0];
        % MATDSS.Meas.VProfileMag = [MATDSS.Meas.VProfileMag, 0];
        % MATDSS.Meas.VProfileAng = [MATDSS.Meas.VProfileAng, 0];
        %
        % % Use OpenDSS measurements and calculations
        % MATDSS.Meas.VMagProfile = [MATDSS.Meas.VMagProfile, 0];
        % MATDSS.Meas.VMagProfilePu = [MATDSS.Meas.VMagProfilePu, 0];
    end
end


% Save all measurements in MATDSS.Sim.Meas if the flag is 1
% MATDSS.Sim.Meas.at = [MATDSS.Sim.Meas.at, at];
MATDSS.Sim.Meas.VMagProfile(:,t) = AllBusVmag;
MATDSS.Sim.Meas.VMagProfilePu(:,t) = AllBusVmagPu;


MATDSS.Sim.Meas.VProfile(:,t) = v;
MATDSS.Sim.Meas.VProfileMag(:,t) = abs(v);
MATDSS.Sim.Meas.VProfileAng(:,t) =angle(v).*180/pi;


end

%}