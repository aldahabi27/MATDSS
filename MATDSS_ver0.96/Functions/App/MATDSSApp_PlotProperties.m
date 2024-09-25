function MATDSSApp_PlotProperties(app)
% MATDSSApp_PlotProperties(app)
% This function is responsible for displaying the correct list of plot
% options that should appear in the Plot Properties panel of the MATDSS
% Application Main Window.
%
% If initialization is missing, this function will display a message to
% prompt the user to either load SimData or run initialization first.
%
% Parameters:
%   - app: The MATDSS application instance.
%
% Last Update for this function was on MATDSS App Ver 0.96 (20 Sept. 2024)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% by Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

% Initialize empty cell arrays for measurement and simulation item names
MeasItemNames = {};
SimItemNames = {};

% Access relevant application data
MATDSS = app.MyRun.MATDSS; % Get the MATDSS object from the application
DER = app.MyRun.DER; % Get the Distributed Energy Resources (DER) object
nVDER = 0; % Initialize count for Virtual Distributed Energy Resources (VDERs)

% Get the number of DERs and control areas
nDER = length(DER); % Count of DERs in the simulation
CA = MATDSS.Cont.CA; % Control Areas in the MATDSS
nCA = length(CA); % Count of Control Areas

% Reset item name lists
MeasItemNames = {};
SimItemNames = {};

% Determine which tab is currently selected to generate corresponding plot options
switch app.PPPTabGroup.SelectedTab.Title
    case [char(916) 'P0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                DeltaP0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generating tree nodes list for DeltaP0 Tab Plot (P0Tiles)

        % Generating the list of curves under measured DeltaP0 Plot
        % Add a group option for P0 that will auto-select all elements within this group
        MeasItemNames = [MeasItemNames; [char(916) 'P0 (Group)']];
        MeasItemNames = [MeasItemNames; [char(916) 'P0']]; % Add P0 curve
        MeasItemNames = [MeasItemNames; [char(916) 'P0,set']]; % Add P0 set point
        MeasItemNames = [MeasItemNames; [char(916) 'P0,set - limits']]; % Add limits for P0 set
        MeasItemNames = [MeasItemNames; [char(916) 'P0,set - Controller limits']]; % Add controller limits for P0 set

        % Indices for measured P0 group
        iMeasP0Group = 1:5; % Index range for the above items

        % Check if any DERs exist and populate measurement item names accordingly
        if nDER >= 1
            MeasItemNames = [MeasItemNames; 'P_DER (Group)']; % Add a group option for DERs
            for i = 1:nDER
                MeasItemNames = [MeasItemNames; ['P_' DER(i).DSSName]]; % Append each DER name to the measurement list
            end
        end
        % Update indices for measured DER group
        iMeasP_DERGroup = [1:nDER+1] + iMeasP0Group(end); % Adjust indices based on existing groups

        % Update the ListBox with measurement items
        app.P0MeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox
        app.P0MeasListBox.Value = app.P0MeasListBox.Items(iMeasP0Group); % Select default values for the measurement ListBox

        % Populate simulation item names for P0
        SimItemNames = [SimItemNames; [char(916) 'P0']]; % Add P0 simulation item
        SimItemNames = [SimItemNames; [char(916) 'P0,set']]; % Add P0 set point for simulation
        app.P0SimListBox.Items = SimItemNames; % Set the items in the simulation ListBox
        app.P0SimListBox.Value = app.P0SimListBox.Items; % Select default values for the simulation ListBox

    case ['CA-' char(916) 'P0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %           CA - DeltaP0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generating tree nodes list for DeltaP0 Tab Plot (CA_P0Tiles)
        % This section creates measurement item names specific to Control Areas (CAs).

        % Loop through each Control Area (CA) to create measurement item names
        for i = 1:nCA
            % Add a group option for each CA, which will auto-select all elements within this group
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0 (group)']];
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0']]; % Add P0 for the CA
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0,set']]; % Add P0 set point for the CA
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0,set - E-limits']]; % Add limits for P0 set
        end

        % Check if any DERs exist and populate measurement item names accordingly
        if nDER >= 1
            MeasItemNames = [MeasItemNames; 'P_DER (Group)']; % Add a group option for DERs
            for i = 1:nDER
                MeasItemNames = [MeasItemNames; ['P_' DER(i).DSSName]]; % Append each DER name to the measurement list
            end
        end

        % Create index groups for measured P0 items
        iMeasP0Group = 1:nCA*4; % Generate index for all items created for CAs
        iMeasP0Group = reshape(iMeasP0Group, 4, [])'; % Reshape indices to group by 4 for each CA

        % Update indices for measured DER group
        iMeasP_DERGroup = [1:nDER+1] + iMeasP0Group(end); % Adjust indices based on existing groups

        % Update the ListBox with measurement items for Control Areas
        app.CAP0MeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox for CAs
        app.CAP0MeasListBox.Value = app.CAP0MeasListBox.Items(reshape(iMeasP0Group(:, [2, 3]), [], 1)); % Select default values for the CA measurement ListBox

        % Populate simulation item names for CA P0
        SimItemNames = [SimItemNames; [char(916) 'P0']]; % Add P0 simulation item for CA
        SimItemNames = [SimItemNames; [char(916) 'P0,set']]; % Add P0 set point for simulation
        app.CAP0SimListBox.Items = SimItemNames; % Set the items in the simulation ListBox for CAs
        % app.CAP0SimListBox.Value = app.CAP0SimListBox.Items; % Optionally set default values (currently commented out)



    case ['P_DER']


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 P_DER Plot Trees                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                DeltaP0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generating tree nodes list for DeltaP0 Tab Plot (P0Tiles)
        % This section creates measurement item names for the DeltaP0 plot.

        % Add a group option for P0, which will auto-select all elements within this group
        MeasItemNames = [MeasItemNames; [char(916) 'P0 (Group)']];
        MeasItemNames = [MeasItemNames; [char(916) 'P0']]; % Add P0 measurement
        MeasItemNames = [MeasItemNames; [char(916) 'P0,set']]; % Add P0 set point
        MeasItemNames = [MeasItemNames; [char(916) 'P0,set - limits']]; % Add limits for P0 set
        MeasItemNames = [MeasItemNames; [char(916) 'P0,set - Controller limits']]; % Add controller limits for P0 set

        % Create an index for measurement items related to P0
        iMeasP0Group = 1:5; % The first five items correspond to P0 measurements

        % Check if any DERs exist and populate measurement item names accordingly
        if nDER >= 1
            MeasItemNames = [MeasItemNames; 'P_DER (Group)']; % Add a group option for DERs
            for i = 1:nDER
                MeasItemNames = [MeasItemNames; ['P_' DER(i).DSSName]]; % Append each DER name to the measurement list
                % Check if the current DER is a VDER and count it
                if strcmp(lower(DER(i).DSSName(1:4)), 'vder')
                    nVDER = nVDER + 1; % Increment VDER count
                end
            end
        end

        % Create index for DER measurement items
        iMeasP_DERGroup = [1:nDER+1] + iMeasP0Group(end); % Adjust indices based on existing groups

        % Update the ListBox with measurement items for DERs
        app.P_DERMeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox for DERs
        app.P_DERMeasListBox.Value = app.P_DERMeasListBox.Items(iMeasP_DERGroup(1:end-nVDER)); % Select default values for the DER measurement ListBox, excluding VDERs

        % Populate simulation item names for P0
        SimItemNames = [SimItemNames; [char(916) 'P0']]; % Add P0 simulation item
        SimItemNames = [SimItemNames; [char(916) 'P0,set']]; % Add P0 set point for simulation
        app.P_DERSimListBox.Items = SimItemNames; % Set the items in the simulation ListBox for DERs
        % app.P_DERSimListBox.Value = app.P_DERSimListBox.Items; % Optionally set default values (currently commented out)

    case [char(916) 'Q0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                DeltaQ0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generating tree nodes list for DeltaQ0 Tab Plot (Q0Tiles)
        % This section creates measurement item names for the DeltaQ0 plot.

        % Add a group option for Q0, which will auto-select all elements within this group
        MeasItemNames = [MeasItemNames; [char(916) 'Q0 (Group)']];
        MeasItemNames = [MeasItemNames; [char(916) 'Q0']]; % Add Q0 measurement
        MeasItemNames = [MeasItemNames; [char(916) 'Q0,set']]; % Add Q0 set point
        MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - limits']]; % Add limits for Q0 set
        MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - Controller limits']]; % Add controller limits for Q0 set

        % Create an index for measurement items related to Q0
        iMeasQ0Group = 1:5; % The first five items correspond to Q0 measurements

        % Check if any DERs exist and populate measurement item names accordingly
        if nDER >= 1
            MeasItemNames = [MeasItemNames; 'Q_DER (Group)']; % Add a group option for DERs
            for i = 1:nDER
                MeasItemNames = [MeasItemNames; ['Q_' DER(i).DSSName]]; % Append each DER name to the measurement list for Q
            end
        end

        % Create index for DER measurement items
        iMeasQ_DERGroup = [1:nDER+1] + iMeasQ0Group(end); % Adjust indices based on existing groups

        % Update the ListBox with measurement items for Q0
        app.Q0MeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox for Q0
        app.Q0MeasListBox.Value = app.Q0MeasListBox.Items(iMeasQ0Group); % Select default values for the Q0 measurement ListBox

        % Populate simulation item names for Q0
        SimItemNames = [SimItemNames; [char(916) 'Q0']]; % Add Q0 simulation item
        SimItemNames = [SimItemNames; [char(916) 'Q0,set']]; % Add Q0 set point for simulation
        app.Q0SimListBox.Items = SimItemNames; % Set the items in the simulation ListBox for Q0
        app.Q0SimListBox.Value = app.Q0SimListBox.Items; % Optionally set default values for the Q0 simulation ListBox


    case ['CA-' char(916) 'Q0']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %           CA - DeltaQ0 Plot Trees                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generating tree nodes list for DeltaQ0 Tab Plot (CA_Q0Tiles)
        % This section creates measurement item names for the CA DeltaQ0 plot.

        % Loop through each CA (Control Area) and generate measurement item names
        for i = 1:nCA
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0 (group)']]; % Group for each CA's Q0
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0']]; % Individual Q0 measurement for each CA
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0,set']]; % Set point for Q0
            MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0,set - E-limits']]; % Error limits for Q0 set
        end

        % Check if any DERs exist and populate measurement item names accordingly
        if nDER >= 1
            MeasItemNames = [MeasItemNames; 'Q_DER (Group)']; % Add a group option for DERs
            for i = 1:nDER
                MeasItemNames = [MeasItemNames; ['Q_' DER(i).DSSName]]; % Append each DER's Q measurement to the list
            end
        end

        % Create index for measurement items related to Q0
        iMeasQ0Group = 1:nCA*4;
        iMeasQ0Group = reshape(iMeasQ0Group, 4, [])'; % Reshape to organize measurement indices

        % Create index for DER measurement items
        iMeasQ_DERGroup = [1:nDER+1] + iMeasQ0Group(end); % Adjust indices based on existing groups

        % Update the ListBox with measurement items for CA Q0
        app.CAQ0MeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox for CA Q0
        app.CAQ0MeasListBox.Value = app.CAQ0MeasListBox.Items(reshape(iMeasQ0Group(:, [2, 3]), [], 1)); % Select default values

        % Populate simulation item names for CA Q0
        SimItemNames = [SimItemNames; [char(916) 'Q0']]; % Add Q0 simulation item
        SimItemNames = [SimItemNames; [char(916) 'Q0,set']]; % Add Q0 set point for simulation
        app.CAQ0SimListBox.Items = SimItemNames; % Set the items in the simulation ListBox for CA Q0
        % app.CAQ0SimListBox.Value = app.Q0SimListBox.Items; % Optionally set default values (currently commented out)

    case ['Q_DER']
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 Q_DER Plot Trees                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This section generates measurement item names for the Q_DER plot.

        % Generating the list of curves under measured DeltaQ0 Plot
        MeasItemNames = [MeasItemNames; [char(916) 'Q0 (Group)']]; % Add group option for Q0 measurements
        MeasItemNames = [MeasItemNames; [char(916) 'Q0']]; % Add Q0 measurement
        MeasItemNames = [MeasItemNames; [char(916) 'Q0,set']]; % Add Q0 set point
        MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - limits']]; % Add limits for Q0 set
        MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - Controller limits']]; % Add controller limits for Q0 set

        % Create an index for measurement items related to Q0
        iMeasQ0Group = 1:5; % The first five items correspond to Q0 measurements

        % Check if any DERs exist and populate measurement item names accordingly
        if nDER >= 1
            MeasItemNames = [MeasItemNames; 'Q_DER (Group)']; % Add a group option for DERs
            for i = 1:nDER
                MeasItemNames = [MeasItemNames; ['Q_' DER(i).DSSName]]; % Append each DER's Q measurement to the list
                % Check if the current DER is a VDER and count it
                if strcmp(lower(DER(i).DSSName(1:4)), 'vder')
                    nVDER = nVDER + 1; % Increment VDER count
                end
            end
        end

        % Create index for DER measurement items
        iMeasQ_DERGroup = [1:nDER+1] + iMeasQ0Group(end); % Adjust indices based on existing groups

        % Update the ListBox with measurement items for Q_DER
        app.Q_DERMeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox for Q_DER
        app.Q_DERMeasListBox.Value = app.Q_DERMeasListBox.Items(iMeasQ_DERGroup(1:end-nVDER)); % Select default values excluding VDERs

        % Populate simulation item names for Q_DER
        SimItemNames = [SimItemNames; [char(916) 'Q0']]; % Add Q0 simulation item for DER
        SimItemNames = [SimItemNames; [char(916) 'Q0,set']]; % Add Q0 set point for DER simulation
        app.Q_DERSimListBox.Items = SimItemNames; % Set the items in the simulation ListBox for Q_DER
        % app.Q_DERSimListBox.Value = app.Q_DERSimListBox.Items; % Optionally set default values (currently commented out)



    case 'VP'
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %               V_Profile Plot Trees               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This section generates measurement item names for the Voltage Profile plot.

        AllNodesNames = MATDSS.Sim.Meas.AllNodesNames; % Retrieve names of all nodes from the simulation data
        k_v = MATDSS.Meas.k_v; % Indexes for measured voltage phases

        % Check if there are any voltage measurements
        if length(k_v) >= 1
            MeasItemNames = [MeasItemNames; 'Measured Phases (Group)']; % Add group option for measured phases
            for i = 1:length(MATDSS.Meas.k_v)
                MeasItemNames = [MeasItemNames; AllNodesNames{k_v(i)}]; % Append each measured phase's name
            end
        end

        % Update the ListBox with measurement items for voltage profile
        app.VMeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox for voltage profile
        app.VMeasListBox.Value = app.VMeasListBox.Items; % Select all items by default

    case 'IP'
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %               I_Profile Plot Trees               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This section generates measurement item names for the Current Profile plot.

        AllLinesNames = MATDSS.Sim.Meas.AllBranchesPhasesNames; % Retrieve names of all branches/lines from simulation data
        k_i = MATDSS.Meas.k_i; % Indexes for measured current phases

        % Check if there are any current measurements
        if length(k_i) >= 1
            MeasItemNames = [MeasItemNames; 'Measured Branches/Lines']; % Add group option for measured branches/lines
            for i = 1:length(k_i)
                MeasItemNames = [MeasItemNames; AllLinesNames{k_i(i)}]; % Append each measured branch/line's name
            end
        end

        % Update the ListBox with measurement items for current profile
        app.IMeasListBox.Items = MeasItemNames; % Set the items in the measurement ListBox for current profile
        app.IMeasListBox.Value = app.IMeasListBox.Items; % Select all items by default

end

end
% end % function end


%% Old function code

%{


    function MATDSSApp_PlotProperties(app)
        % MATDSSApp_PlotProperties(app)
        % This function is responsible for showing up the correct list of plot
        % options that should be shown in the Plot Properties panel in MATDSS
        % Application Main Window.
        %
        %
        % If Initialization is missing, this function will display a msg to ask for
        % either load SimData or run initiallization first.


        MeasItemNames = {};
        SimItemNames = {};

        MATDSS = app.MyRun.MATDSS;
        DER = app.MyRun.DER;
        nVDER = 0;

        nDER = length(DER);
        CA = MATDSS.Cont.CA;
        nCA = length(CA);





        MeasItemNames = {};
        SimItemNames = {};
        switch app.PPPTabGroup.SelectedTab.Title
            case [char(916) 'P0']
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                DeltaP0 Plot Trees                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generating Treenodes list for DeltaP0 Tab Plot (P0Tiles)

                % Generating the list of curves under measured DeltaP0 Plot

                MeasItemNames = [MeasItemNames; [char(916) 'P0 (Group)']]; % selecting this will auto-select all elements within this group
                MeasItemNames = [MeasItemNames; [char(916) 'P0']];
                MeasItemNames = [MeasItemNames; [char(916) 'P0,set']];
                MeasItemNames = [MeasItemNames; [char(916) 'P0,set - limits']];
                MeasItemNames = [MeasItemNames; [char(916) 'P0,set - Controller limits']];

                iMeasP0Group = 1:5;


                % SelMeasItemNames = MeasItemNames;


                if nDER >= 1
                    MeasItemNames = [MeasItemNames; 'P_DER (Group)'];
                    for i = 1:nDER
                        MeasItemNames = [MeasItemNames; ['P_' DER(i).DSSName]];
                    end
                end
                iMeasP_DERGroup = [1:nDER+1] + iMeasP0Group(end);


                app.P0MeasListBox.Items = MeasItemNames;
                app.P0MeasListBox.Value = app.P0MeasListBox.Items(iMeasP0Group);



                SimItemNames = [SimItemNames; [char(916) 'P0']];
                SimItemNames = [SimItemNames; [char(916) 'P0,set']];
                app.P0SimListBox.Items = SimItemNames;
                app.P0SimListBox.Value = app.P0SimListBox.Items;

            case ['CA-' char(916) 'P0']




                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %           CA - DeltaP0 Plot Trees                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generating Treenodes list for DeltaP0 Tab Plot (CA_P0Tiles)
                % MeasItemNames = [MeasItemNames; ['All CAs - ' char(916) 'P0 (Default)']]; % selecting this will auto-select all elements within this group

                for i = 1:nCA
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0 (group)']];
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0']];
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0,set']];
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'P0,set - E-limits']];
                end


                if nDER >= 1
                    MeasItemNames = [MeasItemNames; 'P_DER (Group)'];
                    for i = 1:nDER
                        MeasItemNames = [MeasItemNames; ['P_' DER(i).DSSName]];
                    end
                end


                iMeasP0Group = 1:nCA*4; iMeasP0Group = reshape(iMeasP0Group,4,[])';

                iMeasP_DERGroup = [1:nDER+1] + iMeasP0Group(end);

                app.CAP0MeasListBox.Items = MeasItemNames;
                app.CAP0MeasListBox.Value = app.CAP0MeasListBox.Items(reshape(iMeasP0Group(:,[2,3]),[],1));




                SimItemNames = [SimItemNames; [char(916) 'P0']];
                SimItemNames = [SimItemNames; [char(916) 'P0,set']];
                app.CAP0SimListBox.Items = SimItemNames;
                % app.CAP0SimListBox.Value = app.CAP0SimListBox.Items;



                %%


            case ['P_DER']


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 P_DER Plot Trees                 %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                DeltaP0 Plot Trees                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generating Treenodes list for DeltaP0 Tab Plot (P0Tiles)

                % Generating the list of curves under measured DeltaP0 Plot

                MeasItemNames = [MeasItemNames; [char(916) 'P0 (Group)']]; % selecting this will auto-select all elements within this group
                MeasItemNames = [MeasItemNames; [char(916) 'P0']];
                MeasItemNames = [MeasItemNames; [char(916) 'P0,set']];
                MeasItemNames = [MeasItemNames; [char(916) 'P0,set - limits']];
                MeasItemNames = [MeasItemNames; [char(916) 'P0,set - Controller limits']];

                iMeasP0Group = 1:5;


                % SelMeasItemNames = MeasItemNames;

                if nDER >= 1
                    MeasItemNames = [MeasItemNames; 'P_DER (Group)'];
                    for i = 1:nDER
                        MeasItemNames = [MeasItemNames; ['P_' DER(i).DSSName]];
                        if strcmp(lower(DER(i).DSSName(1:4)),'vder')
                            nVDER = nVDER + 1;
                        end
                    end
                end


                iMeasP_DERGroup = [1:nDER+1] + iMeasP0Group(end);


                app.P_DERMeasListBox.Items = MeasItemNames;
                app.P_DERMeasListBox.Value = app.P_DERMeasListBox.Items(iMeasP_DERGroup(1:end-nVDER));



                SimItemNames = [SimItemNames; [char(916) 'P0']];
                SimItemNames = [SimItemNames; [char(916) 'P0,set']];
                app.P_DERSimListBox.Items = SimItemNames;
                % app.P_DERSimListBox.Value = app.P_DERSimListBox.Items;

            case [char(916) 'Q0']
                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                DeltaQ0 Plot Trees                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generating Treenodes list for DeltaQ0 Tab Plot (Q0Tiles)

                % Generating the list of curves under measured DeltaQ0 Plot

                MeasItemNames = [MeasItemNames; [char(916) 'Q0 (Group)']]; % selecting this will auto-select all elements within this group
                MeasItemNames = [MeasItemNames; [char(916) 'Q0']];
                MeasItemNames = [MeasItemNames; [char(916) 'Q0,set']];
                MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - limits']];
                MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - Controller limits']];

                iMeasQ0Group = 1:5;


                % SelMeasItemNames = MeasItemNames;


                if nDER >= 1
                    MeasItemNames = [MeasItemNames; 'Q_DER (Group)'];
                    for i = 1:nDER
                        MeasItemNames = [MeasItemNames; ['Q_' DER(i).DSSName]];
                    end
                end
                iMeasQ_DERGroup = [1:nDER+1] + iMeasQ0Group(end);


                app.Q0MeasListBox.Items = MeasItemNames;
                app.Q0MeasListBox.Value = app.Q0MeasListBox.Items(iMeasQ0Group);



                SimItemNames = [SimItemNames; [char(916) 'Q0']];
                SimItemNames = [SimItemNames; [char(916) 'Q0,set']];
                app.Q0SimListBox.Items = SimItemNames;
                app.Q0SimListBox.Value = app.Q0SimListBox.Items;

            case ['CA-' char(916) 'Q0']




                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %           CA - DeltaQ0 Plot Trees                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generating Treenodes list for DeltaQ0 Tab Plot (CA_Q0Tiles)
                % MeasItemNames = [MeasItemNames; ['All CAs - ' char(916) 'Q0 (Default)']]; % selecting this will auto-select all elements within this group

                for i = 1:nCA
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0 (group)']];
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0']];
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0,set']];
                    MeasItemNames = [MeasItemNames; ['CA' num2str(CA(i).Area) '_' char(916) 'Q0,set - E-limits']];
                end


                if nDER >= 1
                    MeasItemNames = [MeasItemNames; 'Q_DER (Group)'];
                    for i = 1:nDER
                        MeasItemNames = [MeasItemNames; ['Q_' DER(i).DSSName]];
                    end
                end


                iMeasQ0Group = 1:nCA*4; iMeasQ0Group = reshape(iMeasQ0Group,4,[])';

                iMeasQ_DERGroup = [1:nDER+1] + iMeasQ0Group(end);

                app.CAQ0MeasListBox.Items = MeasItemNames;
                app.CAQ0MeasListBox.Value = app.CAQ0MeasListBox.Items(reshape(iMeasQ0Group(:,[2,3]),[],1));




                SimItemNames = [SimItemNames; [char(916) 'Q0']];
                SimItemNames = [SimItemNames; [char(916) 'Q0,set']];
                app.CAQ0SimListBox.Items = SimItemNames;
                % app.CAQ0SimListBox.Value = app.Q0SimListBox.Items;



                %%


            case ['Q_DER']


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 Q_DER Plot Trees                 %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                DeltaQ0 Plot Trees                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Generating Treenodes list for DeltaQ0 Tab Plot (Q0Tiles)

                % Generating the list of curves under measured DeltaQ0 Plot

                MeasItemNames = [MeasItemNames; [char(916) 'Q0 (Group)']]; % selecting this will auto-select all elements within this group
                MeasItemNames = [MeasItemNames; [char(916) 'Q0']];
                MeasItemNames = [MeasItemNames; [char(916) 'Q0,set']];
                MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - limits']];
                MeasItemNames = [MeasItemNames; [char(916) 'Q0,set - Controller limits']];

                iMeasQ0Group = 1:5;


                % SelMeasItemNames = MeasItemNames;


                if nDER >= 1
                    MeasItemNames = [MeasItemNames; 'Q_DER (Group)'];
                    for i = 1:nDER
                        MeasItemNames = [MeasItemNames; ['Q_' DER(i).DSSName]];
                        if strcmp(lower(DER(i).DSSName(1:4)),'vder')
                            nVDER = nVDER + 1;
                        end
                    end
                end


                iMeasQ_DERGroup = [1:nDER+1] + iMeasQ0Group(end);


                app.Q_DERMeasListBox.Items = MeasItemNames;
                app.Q_DERMeasListBox.Value = app.Q_DERMeasListBox.Items(iMeasQ_DERGroup(1:end-nVDER));



                SimItemNames = [SimItemNames; [char(916) 'Q0']];
                SimItemNames = [SimItemNames; [char(916) 'Q0,set']];
                app.Q_DERSimListBox.Items = SimItemNames;
                % app.Q_DERSimListBox.Value = app.Q_DERSimListBox.Items;


            case 'VP'

                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %               V_Profile Plot Trees               %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                AllNodesNames = MATDSS.Sim.Meas.AllNodesNames;

                k_v = MATDSS.Meas.k_v;

                if length(k_v) >= 1
                    MeasItemNames = [MeasItemNames; 'Measured Phases (Group)']; % selecting this will auto-select all elements within this group
                    for i = 1:length(MATDSS.Meas.k_v)
                        MeasItemNames = [MeasItemNames; AllNodesNames{k_v(i)}];
                    end
                end
                % PPP.Meas.VP.CheckedNodes = PPP.Meas.VP.Children;

                app.VMeasListBox.Items = MeasItemNames;
                app.VMeasListBox.Value = app.VMeasListBox.Items;

            case 'IP'


                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %               I_Profile Plot Trees               %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                AllLinesNames = MATDSS.Sim.Meas.AllBranchesPhasesNames;
                k_i = MATDSS.Meas.k_i;

                if length(k_i) >= 1
                    MeasItemNames = [MeasItemNames; 'Measured Branches/Lines']; % selecting this will auto-select all elements within this group
                    for i = 1:length(k_i)
                        MeasItemNames = [MeasItemNames; AllLinesNames{k_i(i)}];
                    end
                end

                app.IMeasListBox.Items = MeasItemNames;
                app.IMeasListBox.Value = app.IMeasListBox.Items;


        end
    end % function end


%}