function [MATDSS] = MATDSS_IPQ(MATDSS)
% [MATDSS] = MATDSS_IPQ(MATDSS)
% This function calculates the I (current), P (active power), and Q (reactive power) 
% profiles of the entire network. It requires the Y_bus matrix (from MATDSS_YBusMatrix) 
% and the V_Profile (from MATDSS_VProfile).
%
% In addition to obtaining network-wide P and Q profiles, it calculates:
%   - P0 and Q0 at the interface bus
%   - Current profile (I Profile)
%   - Apparent power profile (S Profile)
%
% Parameters:
%   - MATDSS: The MATDSS structure, containing the simulation data and configurations
%
% Notes:
%   - The function updates or initializes variables depending on the current time step.
%   - For each branch in the network, the current (I) and its magnitude (Imax) are computed.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com




% Extract the current time step and initial values for calculation
t = MATDSS.Sim.t;  % current time step
Bus0 = MATDSS.Sim.Bus0;  % Interface bus identifier
at = MATDSS.Sim.at;  % simulation 'at' time
ybus = MATDSS.Sim.DSSYbus;  % Y-bus matrix
v = MATDSS.Sim.Meas.VProfile(:, t);  % Voltage profile at current time step

% Calculate current injection profile using Y-bus and voltage
i = ybus * v;  % Current injection profile
s = v .* conj(i);  % Power profile (S = V * conj(I))

% If it's the first time step, initialize certain parameters
if t == 1
    if ~MATDSS.UseSimDataFlag
        % Get all line names from the DSSCircuit object
        MyLines = MATDSS.Sim.DSSCircuit.Lines;
        MyLinesAllNames = MyLines.AllNames;  % All line names
        MATDSS.Sim.MyLinesAllNames = MyLinesAllNames;  % Save line names to the structure

        % Initialize current profiles and handle missing data
        i_size = -1;
        for j = 1:size(MATDSS.TableData.VI, 1)
            if MATDSS.TableData.VI(j, 3) == "" || ismissing(MATDSS.TableData.VI(j, 3))
                break;
            end
            i_size = j;
        end
        
        % If no valid data is found, set the current indices to zero
        if i_size <= 0
            k_i = 0;
            MATDSS.Meas.IProfileNames = [];
        else
            % Get current indices from the table
            k_i = str2double(MATDSS.TableData.VI(1:i_size, 4));
        end
        MATDSS.Meas.k_i = k_i;

        % Initialize phase names for all branches
        for j = 1:size(MATDSS.TableData.VI, 1)
            if MATDSS.TableData.VI(j, 6) == "" || ismissing(MATDSS.TableData.VI(j, 6))
                break;
            end
            i_size = j;
        end
        MATDSS.Sim.Meas.AllBranchesPhasesNames = MATDSS.TableData.VI(1:i_size, 6);

        % Initialize current matrix size
        ImatSize = zeros(length(MATDSS.Sim.Meas.AllBranchesPhasesNames), 1);
    else
        % Load previous simulation data if available
        MATDSS.Sim.MyLinesAllNames = MATDSS.LoadedSimData.MATDSS.Sim.MyLinesAllNames;
        MATDSS.Meas.k_i = MATDSS.LoadedSimData.MATDSS.Meas.k_i;
        MATDSS.Sim.Meas.AllBranchesPhasesNames = MATDSS.LoadedSimData.MATDSS.Sim.Meas.AllBranchesPhasesNames;
    end
end

% Initialize the current and max current profiles
MyLinesAllNames = MATDSS.Sim.MyLinesAllNames;
k_i = MATDSS.Meas.k_i;
I = cell(length(MyLinesAllNames), 1);  % Initialize current cell
Imax = cell(length(MyLinesAllNames), 1);  % Initialize max current cell

% Loop over each line/branch to calculate the current
for j = 1:length(MyLinesAllNames)
    success = MATDSS.Sim.DSSCircuit.SetActiveElement(['Line.' MyLinesAllNames{j}]);  % Set active element
    if ~success
        error('Error in setting active element in branch current measurements in MATDSS_IPQ.m');
    end
    
    iLine = MATDSS.Sim.DSSCircuit.ActiveElement;  % Get the active element (line)
    iLineNPhases = iLine.NumPhases;  % Number of phases in the line
    iLineCurr = iLine.Currents;  % Get the current on both terminals
    iLineCurr = iLineCurr(1:2:end) + 1i * iLineCurr(2:2:end);  % Convert to complex numbers
    iLineCurrMag = abs(iLineCurr);  % Get the magnitude of current
    iLineCurrMag = [iLineCurrMag(1:iLineNPhases); iLineCurrMag(iLineNPhases+1:end)]; % align each terminal values to pick the max for each phase
    iLineCurrMag = max(iLineCurrMag).';  % Pick the max per phase
    
    % Store current magnitudes
    I(j,:) = {iLineCurrMag}; % save the currents in I cell
    if t == 1
        Imax(j,:) = {iLine.NormalAmps.*ones(iLineNPhases,1)};  % Store normal current values for first step
    end
end

% At the first time step, initialize maximum currents and measurement names
if t == 1
    MATDSS.Sim.Meas.Imax = Imax;
    if k_i > 0
        Imaxmat = cell2mat(Imax);
        Imaxmat(Imaxmat == 0) = 1e99;  % Set unlimited lines to a very high value
        MATDSS.Meas.Imax = Imaxmat(k_i);
    end
    MATDSS.Sim.Meas.DSSIProfileNames = MyLinesAllNames;  % Save branch names

    % Initialize profiles for current, power, and apparent power
    MATDSS.Meas.IinjProfile = zeros(length(i), length(MATDSS.Meas.at));  
    MATDSS.Meas.SProfile = zeros(length(s), length(MATDSS.Meas.at));
    MATDSS.Meas.PProfile = zeros(length(real(s)), length(MATDSS.Meas.at));
    MATDSS.Meas.QProfile = zeros(length(imag(s)), length(MATDSS.Meas.at));

    % Initialize interface bus power profiles
    MATDSS.Meas.P0 = zeros(length(MATDSS.Meas.PProfile(Bus0, end)), length(MATDSS.Meas.at));
    MATDSS.Meas.Q0 = zeros(length(MATDSS.Meas.QProfile(Bus0, end)), length(MATDSS.Meas.at));

    % Initialize branch current profiles
    Imat = cell2mat(I);
    MATDSS.Meas.IProfile = zeros(length(Imat(k_i)), length(MATDSS.Meas.at));


    % Do the same for MATDSS.Sim.Meas
    MATDSS.Sim.Meas.IinjProfile = zeros(length(i),length(MATDSS.Sim.Meas.at)); % save current injection profile
    MATDSS.Sim.Meas.SProfile = zeros(length(s),length(MATDSS.Sim.Meas.at));
    MATDSS.Sim.Meas.PProfile = zeros(length(real(s)),length(MATDSS.Sim.Meas.at));
    MATDSS.Sim.Meas.QProfile = zeros(length(imag(s)),length(MATDSS.Sim.Meas.at));

    MATDSS.Sim.Meas.P0 = zeros(length(MATDSS.Meas.PProfile(Bus0,end)),length(MATDSS.Sim.Meas.at));
    MATDSS.Sim.Meas.Q0 = zeros(length(MATDSS.Meas.QProfile(Bus0,end)),length(MATDSS.Sim.Meas.at));
    Imat = cell2mat(I);
    MATDSS.Sim.Meas.IProfile = zeros(length(Imat),length(MATDSS.Sim.Meas.at)); % branch currents calculations.


end

% Update measurements at the current time step
if t == 1 || mod(at, MATDSS.Time.Meas.TimeStep) == 0
    at_index = find(MATDSS.Meas.at == at);
    % Save values to the MATDSS struct
    MATDSS.Meas.IinjProfile(:, at_index) = i;
    MATDSS.Meas.SProfile(:, at_index) = s;
    MATDSS.Meas.PProfile(:, at_index) = real(s);
    MATDSS.Meas.QProfile(:, at_index) = imag(s);
    MATDSS.Meas.P0(:, at_index) = MATDSS.Meas.PProfile(Bus0, at_index);
    MATDSS.Meas.Q0(:, at_index) = MATDSS.Meas.QProfile(Bus0, at_index);
    Imat = cell2mat(I);
    MATDSS.Meas.IProfile(:, at_index) = Imat(k_i);
end

% Save values to the MATDSS.Sim structure for the current time step
MATDSS.Sim.Meas.IinjProfile(:, t) = i;
MATDSS.Sim.Meas.SProfile(:, t) = s;
MATDSS.Sim.Meas.PProfile(:, t) = real(s);
MATDSS.Sim.Meas.QProfile(:, t) = imag(s);
MATDSS.Sim.Meas.P0(:, t) = MATDSS.Sim.Meas.PProfile(Bus0, t);
MATDSS.Sim.Meas.Q0(:, t) = MATDSS.Sim.Meas.QProfile(Bus0, t);
MATDSS.Sim.Meas.IProfile(:, t) = cell2mat(I);

end


%% Old function code

%{

function [MATDSS] = MATDSS_IPQ(MATDSS)
% This function will obtain the I, P and Q profile of the whole network. The
% function requries:
% * MATDSS.Sim.Y_bus (can be obtained by running MATDSS_YBusMatrix function)
% * MATDSS.V_Profile (can be obtained by running MATDSS_VProfile function)
%
% In addition to network P and Q profiles, this function also obtains the
% following:
% * P0 and Q0 at the interface bus
% * Current Profile (I Profile)
% * Apparent Power Profile (S Profile)



% Calculating the Power profile at the last time step (last measured
% voltage)

t = MATDSS.Sim.t; % current time step;
Bus0 = MATDSS.Sim.Bus0;
at = MATDSS.Sim.at;
ybus = MATDSS.Sim.DSSYbus; % ybus matrix
v = MATDSS.Sim.Meas.VProfile(:,t); % last voltage measurements
i = ybus*v; % Current injection Profile
s = v.*conj(i); % Power Profile



if t == 1
    if ~MATDSS.UseSimDataFlag
        % Get all lines/branches names
        MyLines = MATDSS.Sim.DSSCircuit.Lines;
        MyLinesAllNames = MyLines.AllNames;
        MATDSS.Sim.MyLinesAllNames = MyLinesAllNames;
        i_size = -1;
        for j = 1:size(MATDSS.TableData.VI,1)
            if MATDSS.TableData.VI(j,3) == "" || ismissing(MATDSS.TableData.VI(j,3))
                break;
            end
            i_size = j;
        end
        if i_size <= 0
            k_i = 0;
            MATDSS.Meas.IProfileNames = [];
        else
            k_i = str2double(MATDSS.TableData.VI(1:i_size,4));
            %         MATDSS.Meas.IProfileNames = MATDSS.TableData.VI(1:i_size,3);
        end
        MATDSS.Meas.k_i = k_i;


        for j = 1:size(MATDSS.TableData.VI,1)
            if MATDSS.TableData.VI(j,6) == "" || ismissing(MATDSS.TableData.VI(j,6))
                break;
            end
            i_size = j;
        end
        MATDSS.Sim.Meas.AllBranchesPhasesNames = MATDSS.TableData.VI(1:i_size,6);


        ImatSize = zeros(length(MATDSS.Sim.Meas.AllBranchesPhasesNames),1);

    else
        MATDSS.Sim.MyLinesAllNames = MATDSS.LoadedSimData.MATDSS.Sim.MyLinesAllNames;
        MATDSS.Meas.k_i = MATDSS.LoadedSimData.MATDSS.Meas.k_i;
        MATDSS.Sim.Meas.AllBranchesPhasesNames = MATDSS.LoadedSimData.MATDSS.Sim.Meas.AllBranchesPhasesNames;
    end

end

MyLinesAllNames = MATDSS.Sim.MyLinesAllNames;
k_i = MATDSS.Meas.k_i;
% disp('i')
I = cell(length(MyLinesAllNames),1);
Imax =  cell(length(MyLinesAllNames),1);

% DSSLines = MATDSS.Sim.DSSCircuit.Lines;
% i_Line = DSSLines.First;
% i_Line_phases = {};
% i_Line_Names = {};
% while i_Line > 0
for j = 1:length(MyLinesAllNames) % loop over each branch to measure/calculate the current
    success = MATDSS.Sim.DSSCircuit.SetActiveElement(['Line.' MyLinesAllNames{j}]); % set active element to line i
    if ~success
        error('Error in setting active element in Branch current measurements in MATDSS_IPQ.m');
    end

    iLine = MATDSS.Sim.DSSCircuit.ActiveElement; % handle of the current line i
    %     iLine.NodeOrder
    iLineNPhases = iLine.NumPhases; % number of phases in this branch/line
    % i_Line_phases = [i_Line_phases; iLineNPhases];
    % i_Line_Names = [i_Line_Names; {iLine.Name}];
    iLineCurr = iLine.Currents; % get the current measurements (on both sides/terminals)
    iLineCurr = iLineCurr(1:2:end) + 1i.*iLineCurr(2:2:end); % convert it to complex
    iLineCurrMag = abs(iLineCurr); % get the magnitude
    iLineCurrMag = [iLineCurrMag(1:iLineNPhases); iLineCurrMag(iLineNPhases+1:end)]; % align each terminal values to pick the max for each phase
    iLineCurrMag = max(iLineCurrMag).';   % get the max per phase



    I(j,:) = {iLineCurrMag}; % save the currents in I cell
    if t == 1
        Imax(j,:) = {iLine.NormalAmps.*ones(iLineNPhases,1)};
    end
    %     Imax = cellfun(@(x) x.*0.95, Imax,'UniformOutput',false);
    % i_Line = DSSLines.Next;
end

% disp('ii')

%
if t == 1
    MATDSS.Sim.Meas.Imax = Imax;
    if k_i > 0
        Imaxmat = cell2mat(Imax);
        Imaxmat(Imaxmat == 0) = 1e99; % for any unlimited line, set it to 1e99
        MATDSS.Meas.Imax = Imaxmat(k_i);
    end

    MATDSS.Sim.Meas.DSSIProfileNames = MyLinesAllNames; % Branch names



    % Initiallizing the variables
    MATDSS.Meas.IinjProfile = zeros(length(i), length(MATDSS.Meas.at)); % save current injection profile
    MATDSS.Meas.SProfile = zeros(length(s), length(MATDSS.Meas.at));
    MATDSS.Meas.PProfile = zeros(length(real(s)), length(MATDSS.Meas.at));
    MATDSS.Meas.QProfile = zeros(length(imag(s)), length(MATDSS.Meas.at));

    MATDSS.Meas.P0 = zeros(length(MATDSS.Meas.PProfile(Bus0,end)), length(MATDSS.Meas.at));
    MATDSS.Meas.Q0 = zeros(length(MATDSS.Meas.QProfile(Bus0,end)), length(MATDSS.Meas.at));
    Imat = cell2mat(I);
    MATDSS.Meas.IProfile = zeros(length(Imat(k_i)), length(MATDSS.Meas.at)); % branch currents calculations.



    MATDSS.Sim.Meas.IinjProfile = zeros(length(i),length(MATDSS.Sim.Meas.at)); % save current injection profile
    MATDSS.Sim.Meas.SProfile = zeros(length(s),length(MATDSS.Sim.Meas.at));
    MATDSS.Sim.Meas.PProfile = zeros(length(real(s)),length(MATDSS.Sim.Meas.at));
    MATDSS.Sim.Meas.QProfile = zeros(length(imag(s)),length(MATDSS.Sim.Meas.at));

    MATDSS.Sim.Meas.P0 = zeros(length(MATDSS.Meas.PProfile(Bus0,end)),length(MATDSS.Sim.Meas.at));
    MATDSS.Sim.Meas.Q0 = zeros(length(MATDSS.Meas.QProfile(Bus0,end)),length(MATDSS.Sim.Meas.at));
    Imat = cell2mat(I);
    MATDSS.Sim.Meas.IProfile = zeros(length(Imat),length(MATDSS.Sim.Meas.at)); % branch currents calculations.


end

% disp('iii')
if t == 1 || mod(at,MATDSS.Time.Meas.TimeStep) == 0  % time to update the measurements
    at_index = find(MATDSS.Meas.at == at);
    % Save the values to the MATDSS Struct
    MATDSS.Meas.IinjProfile(:,at_index) = [i]; % save current injection profile
    MATDSS.Meas.SProfile(:,at_index)  = [s];
    MATDSS.Meas.PProfile(:,at_index)  = [real(s)];
    MATDSS.Meas.QProfile(:,at_index)  = [imag(s)];

    MATDSS.Meas.P0(:,at_index)  = [MATDSS.Meas.PProfile(Bus0,at_index)];
    MATDSS.Meas.Q0(:,at_index)  = [MATDSS.Meas.QProfile(Bus0,at_index)];
    Imat = cell2mat(I);
    MATDSS.Meas.IProfile(:,at_index)  = [Imat(k_i)]; % branch currents calculations.
end
% disp('iv')
MATDSS.Sim.Meas.IinjProfile(:,t) = [i]; % save current injection profile
MATDSS.Sim.Meas.SProfile(:,t) = [s];
MATDSS.Sim.Meas.PProfile(:,t) = [real(s)];
MATDSS.Sim.Meas.QProfile(:,t) = [imag(s)];

MATDSS.Sim.Meas.P0(:,t) = [MATDSS.Sim.Meas.PProfile(Bus0,t)];
MATDSS.Sim.Meas.Q0(:,t) = [MATDSS.Sim.Meas.QProfile(Bus0,t)];

MATDSS.Sim.Meas.IProfile(:,t) = [cell2mat(I)]; % branch currents calculations.
% disp('v')

end

%}