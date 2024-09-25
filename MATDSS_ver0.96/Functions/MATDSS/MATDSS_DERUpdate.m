function [MATDSS, DER] = MATDSS_DERUpdate(MATDSS, DER, DERIndex)
% MATDSS_DERUpdate(app)
% This function is responsible for updating the load details and information
% related to Distributed Energy Resources (DER) in the OpenDSS circuit.
% The function takes into account the last setpoints stored in DER(i).kW_setpoint
% and DER(i).kvar_setpoint, which serve as inputs for dynamic simulations.
%
% Parameters:
%   - MATDSS: The main MATDSS structure containing simulation, control, and measurement data.
%   - DER: Structure containing all details of the Distributed Energy Resources (DERs) involved in the system.
%   - DERIndex: The index of the DER being updated.
%
% The function updates the status of the specific DER, calculates gradients for
% active (P) and reactive (Q) power, and considers the control area's parameters
% (like dual variables and multipliers) to adjust the setpoints accordingly.
%
% Error handling ensures that if an error occurs during the update process,
% it will not crash the simulation. This function is critical to maintaining
% real-time dynamic control over the DERs in the MATDSS application.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application




% Assigning the current simulation time and active time step
t = MATDSS.Sim.t;         % Current simulation time step
at = MATDSS.Sim.at;       % Active simulation time step
i = DERIndex;             % Index of the DER being updated
Tau = DER(i).Tau;         % Time constant of the specific DER

% Creating the time window for this update, using the current time step
tw = [MATDSS.Time.Sim.TimeSpan(t); MATDSS.Time.Sim.TimeSpan(t) + MATDSS.Time.Sim.TimeStep];

% Fetch the active and reactive power setpoints (W and Var) from the DER structure
i_setpoint = DER(i).i_setpoint; % Current setpoint for the DER
i_WVar = DER(i).i_WVar;         % Indexes for the W and Var setpoints
x = [DER(i).W(i_WVar,2); DER(i).Var(i_WVar,2)]; % Store the previous setpoints

% If the DER's status is marked for an update, execute the following block
if DER(i).UpdateStatus
    DER(i).UpdateStatus = false;  % Reset the update status to prevent redundant updates

    % Loop through all control areas (CAs) to find the one this DER belongs to
    for j = 1:length(MATDSS.Cont.CA)
        CADERs = MATDSS.Cont.CA(j).DERIndex;   % Retrieve the DER indexes in the current CA
        if ismember(i, CADERs)                 % Check if this DER belongs to the current CA
            CAi = find(CADERs == i);           % Find the DER index within the CA
            break;                             % Exit the loop once the CA is found
        end
    end

    % Extract the control area index and control parameters
    ix = (CAi - 1) * 2 + 1;        % Index offset for control area
    CAIndex = DER(i).CAIndex;      % Control area index for this DER
    CA = MATDSS.Cont.CA(CAIndex);  % Retrieve the control area structure

    Contt = CA.TimeIndex;          % Time index for the control area

    % Initialize gradients for active and reactive power
    f_Grad_Px = 0;  % Gradient for active power
    f_Grad_Qx = 0;  % Gradient for reactive power

    % Compute the gradient for active power (P) using previous setpoints
    for j = 1:length(DER(i).Px) - 1
        f_Grad_Px = f_Grad_Px + DER(i).Px(j) * (length(DER(i).Px) - j);  % Sum the weighted setpoints
    end

    % Compute the gradient for reactive power (Q) using previous setpoints
    for j = 1:length(DER(i).Qx) - 1
        f_Grad_Qx = f_Grad_Qx + DER(i).Qx(j) * (length(DER(i).Qx) - j);  % Sum the weighted setpoints
    end

    % Retrieve control area parameters related to optimization
    alpha = CA.alpha;          % Alpha control parameter
    rp = CA.rp;                % rp multiplier in the control scheme

    % Retrieve the dual variables used for constraint satisfaction in optimization
    rho = CA.Duals.rho(:, Contt);
    lambda = CA.Duals.lambda(:, Contt);
    mu = CA.Duals.mu(:, Contt);
    gamma = CA.Duals.gamma(:, Contt);
    nu = CA.Duals.nu(:, Contt);
    zeta = CA.Duals.zeta(:, Contt);
    sigma = CA.Duals.sigma(:, Contt);
    eta = CA.Duals.eta(:, Contt);
    psi = CA.Duals.psi(:, Contt);
    Errors = CA.Duals.Errors;
    ErrorsY = CA.Duals.ErrorsY;
    DControlError = Errors;    % Error term in dynamic control

    % Extract the control area matrices used for control law computations
    A = CA.ABMH.A;
    B = CA.ABMH.B;
    M = CA.ABMH.M;
    H = CA.ABMH.H;
    s = MATDSS.Cont.s;         % Some state variable or control multiplier

    % Compute the modified control parameters for the current DER
    alpha_xi = DER(i).ax .* alpha;   % Control parameter alpha modified for the DER
    rp_xi = DER(i).cx .* rp;         % Control parameter rp modified for the DER

    % Retrieve setpoints for active (P) and reactive (Q) power at the current time step
    P0Set = CA.P0Set(Contt);    % Active power setpoint for the control area at time t
    Q0Set = CA.Q0Set(Contt);    % Reactive power setpoint for the control area at time t

    % Determine if reactive power tracking should be disabled for local control areas (llc)
    if strcmpi(MATDSS.Cont.CA(CAIndex).Type, 'llc')
        DisableQTrackingFlag = 0;  % Keep reactive power tracking enabled
    else
        DisableQTrackingFlag = 1;  % Disable reactive power tracking for non-llc control areas
    end

    % This code is to fix the case when no current or voltage control.
    if size(A, 2) < ix + 1
        Aix = size(A, 2);
    else
        Aix = ix:ix + 1;
    end
    if size(B, 2) < ix + 1
        Bix = size(B, 2);
    else
        Bix = ix:ix + 1;
    end


    % Check if P Controller is enabled and if the simulation time is greater than 0.2
    if MATDSS.Cont.PID.EnablePFlag && at > 0.2
        % P Controller gain and parameter settings
        Kp = MATDSS.Cont.PID.Kp;
        kappa_p_lambda = MATDSS.Cont.PID.kappa_p.lambda;
        kappa_p_mu = MATDSS.Cont.PID.kappa_p.mu;
        kappa_p_eta = MATDSS.Cont.PID.kappa_p.eta;
        kappa_p_psi = MATDSS.Cont.PID.kappa_p.psi;
        kappa_p_gamma = MATDSS.Cont.PID.kappa_p.gamma;
        kappa_p_nu = MATDSS.Cont.PID.kappa_p.nu;
        kappa_p_zeta = MATDSS.Cont.PID.kappa_p.zeta;
    else
        % If not enabled, set the gain and parameters to zero
        Kp = 0;
        kappa_p_lambda = 0;
        kappa_p_mu = 0;
        kappa_p_eta = 0;
        kappa_p_psi = 0;
        kappa_p_gamma = 0;
        kappa_p_nu = 0;
        kappa_p_zeta = 0;
    end

    % Check if D Controller is enabled and if the simulation time is greater than 0.2
    if MATDSS.Cont.PID.EnableDFlag && at > 0.2
        % D Controller gain and parameter settings
        Kd = MATDSS.Cont.PID.Kd;
        kappa_d_lambda = MATDSS.Cont.PID.kappa_d.lambda;
        kappa_d_mu = MATDSS.Cont.PID.kappa_d.mu;
        kappa_d_eta = MATDSS.Cont.PID.kappa_d.eta;
        kappa_d_psi = MATDSS.Cont.PID.kappa_d.psi;
        kappa_d_gamma = MATDSS.Cont.PID.kappa_d.gamma;
        kappa_d_nu = MATDSS.Cont.PID.kappa_d.nu;
        kappa_d_zeta = MATDSS.Cont.PID.kappa_d.zeta;

        % Determine the source of D control error signal
        switch lower(MATDSS.Cont.PID.Dsignal)
            case 'error'
                DControlError = Errors; % Use general control error
            case 'y'
                DControlError = ErrorsY; % Use specific Y error
            otherwise
                DControlError = Errors; % Default to general error
        end
    else
        % If not enabled, set the gain and parameters to zero
        Kd = 0;
        kappa_d_lambda = 0;
        kappa_d_mu = 0;
        kappa_d_eta = 0;
        kappa_d_psi = 0;
        kappa_d_gamma = 0;
        kappa_d_nu = 0;
        kappa_d_zeta = 0;
    end

    % Time step for the control
    Ts = MATDSS.Time.Cont.TimeStep;

    % Initialize new setpoints for W and Var
    xnew = x;

    % Determine whether P and/or D control is enabled for this DER/VDER
    PEnable = 0;
    DEnable = 0;
    switch true
        case strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'DERs & VDERs') && strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'DERs & VDERs')
            PEnable = 1; % Enable both P and D control for DERs and VDERs
            DEnable = 1;
        case strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'VDERs only') && strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'DERs & VDERs')
            DEnable = 1; % Enable D control for both, P control only for VDERs
            if strcmpi(DER(i).Type, 'vder')
                PEnable = 1;
            end
        case strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'VDERs only') && strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'DERs & VDERs')
            PEnable = 1; % Enable P control for both, D control only for VDERs
            if strcmpi(DER(i).Type, 'vder')
                DEnable = 1;
            end
        case strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'VDERs only') && strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'VDERs only')
            if strcmpi(DER(i).Type, 'vder')
                DEnable = 1; % Enable both P and D control only for VDERs
                PEnable = 1;
            end
        otherwise
            PEnable = 0; % Disable both P and D control
            DEnable = 0;
    end

    % PID Controller adjustments for the control variables
    PIDrho = rho;
    PIDsigma = sigma;

    % Update PID-controlled variables (lambda, mu, eta, psi, gamma, nu, zeta) based on P and D controllers
    PIDlambda = PEnable*Kp*kappa_p_lambda*Errors.lambda(Contt) + lambda + DEnable*Kd*kappa_d_lambda*(DControlError.lambda(Contt)-DControlError.lambda(Contt-1))./Ts;
    PIDmu = PEnable*Kp*kappa_p_mu*Errors.mu(Contt) + mu + DEnable*Kd*kappa_d_mu*(DControlError.mu(Contt)-DControlError.mu(Contt-1))./Ts;
    PIDeta = PEnable*Kp*kappa_p_eta*Errors.eta(Contt) + eta + DEnable*Kd*kappa_d_eta*(DControlError.eta(Contt)-DControlError.eta(Contt-1))./Ts;
    PIDpsi = PEnable*Kp*kappa_p_psi*Errors.psi(Contt) + psi + DEnable*Kd*kappa_d_psi*(DControlError.psi(Contt)-DControlError.psi(Contt-1))./Ts;
    PIDgamma = PEnable*Kp*kappa_p_gamma*Errors.gamma(:,Contt) + gamma + DEnable*Kd*kappa_d_gamma*(DControlError.gamma(:,Contt)-DControlError.gamma(:,Contt-1))./Ts;
    PIDnu = PEnable*Kp*kappa_p_nu*Errors.nu(:,Contt) + nu + DEnable*Kd*kappa_d_nu*(DControlError.nu(:,Contt)-DControlError.nu(:,Contt-1))./Ts;
    PIDzeta = PEnable*Kp*kappa_p_zeta*Errors.zeta(:,Contt) + zeta + DEnable*Kd*kappa_d_zeta*(DControlError.zeta(:,Contt)-DControlError.zeta(:,Contt-1))./Ts;

    % Only perform the control loop if the condition for VDER type and time step is met
    if (strcmpi(DER(i).Type, 'vder') && mod(round(at,4),round(MATDSS.Time.Cont.TimeStep*MATDSS.Cont.RoU,4)) == 0) || ~strcmpi(DER(i).Type, 'vder')
        % Iterate over the control loop size
        LoopSize = MATDSS.Cont.LoopSize;
        for k = 1:LoopSize
            % Adjust the control variables for W and Var
            xnew = xnew - [DER(i).W(1,2); DER(i).Var(1,2)];
            xnew = xnew - alpha_xi .* ([f_Grad_Px; f_Grad_Qx].* xnew + ...
                (M(:,ix:ix+1).' * (s.*(1e-0*(PIDlambda - PIDmu)) * ones(size(rho)) + PIDrho*sign(-P0Set))) + ...
                (1-DisableQTrackingFlag) .* (H(:,ix:ix+1).' * (s.*(1e-0.*(PIDeta-PIDpsi)) * ones(size(sigma)) + PIDsigma*sign(-Q0Set))) + ...
                (A(:,Aix).' * (1e-0.*(PIDgamma - PIDnu))) + ...
                (-B(:,Bix).' * (1e-0.*PIDzeta)) + ...
                rp_xi .* xnew);

            % Ensure the values do not exceed the DER's max/min setpoints
            xnew = xnew + [DER(i).W(1,2); DER(i).Var(1,2)];
            if xnew(1) > DER(i).Pmax
                xnew(1) = sign(xnew(1)) * DER(i).Pmax;
            elseif xnew(1) < DER(i).Pmin
                xnew(1) = sign(xnew(1)) * DER(i).Pmin;
            end

            if xnew(2) > DER(i).Qmax
                xnew(2) = sign(xnew(2)) * DER(i).Qmax;
            elseif xnew(2) < DER(i).Qmin
                xnew(2) = sign(xnew(2)) * DER(i).Qmin;
            end
        end
    end


    i_setpoint = i_setpoint + 1; % Increment the setpoint index for the DER
    DER(i).i_setpoint = i_setpoint; % Update the DER object with the new setpoint index

    % If controllers are enabled and either the DER is a VDER or not (general case)
    if ~MATDSS.Cont.DisableControllersFlag && (strcmpi(DER(i).Type, 'vder') || ~strcmpi(DER(i).Type, 'vder'))

        % If the DER is of type 'vder' (Virtual DER)
        if strcmpi(DER(i).Type, 'vder')
            % Store the unfiltered active power (W) and reactive power (Var) setpoints
            DER(i).WUnfilteredSetpoint(i_setpoint) = xnew(1);
            DER(i).VarUnvfilteredSetpoint(i_setpoint) = xnew(2);

            % Check if the Low Pass Filter (LPF) is enabled and the simulation time is greater than 0.2
            if MATDSS.Cont.LPF.EnableLPFFlag && at > 0.2
                % Retrieve previous power output and setpoints
                y1 = [DER(i).WSetpoint(i_setpoint-1,2); DER(i).VarSetpoint(i_setpoint-1,2)]; % Previous output power
                x1 = [DER(i).WUnfilteredSetpoint(i_setpoint-1); DER(i).VarUnvfilteredSetpoint(i_setpoint-1)]; % Previous unfiltered setpoints
                x0 = [DER(i).WUnfilteredSetpoint(i_setpoint); DER(i).VarUnvfilteredSetpoint(i_setpoint)]; % Current unfiltered setpoints

                % LPF coefficients for filtering
                d1d0 = MATDSS.Cont.LPF.d1d0;
                n0d0 = MATDSS.Cont.LPF.n0d0;
                n1d0 = MATDSS.Cont.LPF.n1d0;

                % Apply the low-pass filter to smoothen the control signals
                y0 = y1 .* -d1d0 + x0.*n0d0 + x1.*n1d0;

                % Update the setpoint with the filtered values
                xnew = y0;
            end
        end

        % Store the final active and reactive power setpoints with the current time
        DER(i).WSetpoint(i_setpoint,:) = [at, xnew(1)];
        DER(i).VarSetpoint(i_setpoint,:) = [at, xnew(2)];

    else % If the RoU cycle for the VDER is not complete
        % Use the previous setpoints for both active and reactive power
        DER(i).WSetpoint(i_setpoint,:) = [at, DER(i).WSetpoint(i_setpoint-1,2)];
        DER(i).VarSetpoint(i_setpoint,:) = [at, DER(i).VarSetpoint(i_setpoint-1,2)];
    end

end

% If the DER is not of type 'vder', handle its dynamics
if ~strcmpi(DER(i).Type,'vder')

    % Get a handle to the OpenDSS circuit and load objects
    my_circut = MATDSS.Sim.DSSCircuit; % Handle for DSS_Circuit
    myloads = my_circut.Loads; % Handle for DSS loads
    myloads.Name = DER(i).DSSName; % Select the current DER based on its name in DSS

    % Define the new setpoints for active (W) and reactive (Var) power
    x_set = [DER(i).WSetpoint(i_setpoint,2); DER(i).VarSetpoint(i_setpoint,2)]; % New setpoints
    x0 = [DER(i).W(i_WVar,2); DER(i).Var(i_WVar,2)]; % Initial power values

    % Use the ODE solver to model the dynamic response of the DER to the new setpoints
    [ODETime, ODEx] = ode45(@(t,x) FOODE(t,x,x0, Tau, x_set), tw, x0);

    % Store the ODE results (time and state) in the DER object
    DER(i).x.ODETime = [DER(i).x.ODETime; ODETime];
    DER(i).x.ODEx = [DER(i).x.ODEx; ODEx];

    % If controllers are not disabled, update the output powers from the ODE results
    if ~MATDSS.Cont.DisableControllersFlag
        DER(i).W(i_WVar+1,2) = ODEx(end,1); % Update active power (W)
        DER(i).Var(i_WVar+1,2) = ODEx(end,2); % Update reactive power (Var)
    else
        % If controllers are disabled, keep the previous power values
        DER(i).W(i_WVar+1,2) = DER(i).W(i_WVar,2);
        DER(i).Var(i_WVar+1,2) = DER(i).Var(i_WVar,2);
    end

    % Increment the power index for tracking
    DER(i).i_WVar = i_WVar + 1;

    % Update the OpenDSS load with the new power values
    myloads.Name = DER(i).DSSName; % Re-select the DER to be updated
    myloads.kW = -DER(i).W(i_WVar+1,2) / 1e3; % Update active power in OpenDSS (convert to kW)
    myloads.kvar = -DER(i).Var(i_WVar+1,2) / 1e3; % Update reactive power in OpenDSS (convert to kvar)

else

    % If the DER is of type 'vder', handle the power output based on the CA (Control Area)
    i_ = strfind(DER(i).DSSName,'_'); % Find the index of underscores in the DSS name
    CAIndex = str2num(DER(i).DSSName(i_(end)+3:end)); % Extract the Control Area index from the name
    CA = MATDSS.Cont.CA(CAIndex); % Get the Control Area information

    Contt = CA.TimeIndex; % Get the time index for the control area
    DER(i).W(i_WVar+1,2) = [-sum(CA.P0(:,Contt))]; % Update active power based on CA data
    DER(i).Var(i_WVar+1,2) = [-sum(CA.Q0(:,Contt))]; % Update reactive power based on CA data

    % Increment the power index for tracking
    DER(i).i_WVar = i_WVar + 1;

end

end

% ODE function to model the change in power setpoints
function dxdt = FOODE(t,x,x0, Tau, x_set)
% t and x are ODE variables (time and state)
% x0 is the initial state (previous power values)
% Tau is the time constant for the ODE (controls the rate of change)
% x_set is the 2x1 setpoints vector [P;Q] for active and reactive power

% Ensure Tau is a diagonal matrix (if it's not already)
if size(Tau,1) ~= size(Tau,2)
    Tau = diag(Tau);
end

% Define the ODE's rate of change based on the difference between current and setpoint
dxdt = -inv(Tau) * (x - x_set);

end




%% Old Function Code

%{

function [MATDSS, DER] = MATDSS_DERUpdate(MATDSS,DER,DERIndex)
%
%This function is responsible for updating the load details/information in
%OpenDSS circuit.

%The function will consider the last setpoints values in DER(i).kW_setpoint
%and DER(i).kvar_setpoint as the input/set point for the dyanmics ODE.
%


% Setting some easy callback variables
t = MATDSS.Sim.t;
at = MATDSS.Sim.at;
i = DERIndex; % index of DER to be updated
Tau = DER(i).Tau;
tw = [MATDSS.Time.Sim.TimeSpan(t); MATDSS.Time.Sim.TimeSpan(t) + MATDSS.Time.Sim.TimeStep]; % time window for the simulation


i_setpoint = DER(i).i_setpoint;
i_WVar = DER(i).i_WVar;
x = [DER(i).W(i_WVar,2); DER(i).Var(i_WVar,2)]; % W and Var previous setpoints



if DER(i).UpdateStatus
    DER(i).UpdateStatus = false;
    for j = 1:length(MATDSS.Cont.CA)
        CADERs = MATDSS.Cont.CA(j).DERIndex;
        if ismember(i,CADERs)
            CAi = find(CADERs == i);
            break;
        end
    end
    ix = (CAi-1)*2 + 1;
    CAIndex = DER(i).CAIndex;
    CA = MATDSS.Cont.CA(CAIndex);

    Contt = CA.TimeIndex;

    % calculating gradient of f
    f_Grad_Px = 0;
    f_Grad_Qx = 0;


    for j = 1:length(DER(i).Px)-1
        f_Grad_Px = f_Grad_Px + DER(i).Px(j)*(length(DER(i).Px)-j);
    end
    for j = 1:length(DER(i).Qx)-1
        f_Grad_Qx = f_Grad_Qx + DER(i).Qx(j)*(length(DER(i).Qx)-j);
    end


    alpha = CA.alpha;
    rp = CA.rp;

    rho = CA.Duals.rho(:,Contt);
    lambda = CA.Duals.lambda(:,Contt);
    mu = CA.Duals.mu(:,Contt);
    gamma = CA.Duals.gamma(:,Contt);
    nu = CA.Duals.nu(:,Contt);
    zeta = CA.Duals.zeta(:,Contt);

    sigma = CA.Duals.sigma(:,Contt);
    eta = CA.Duals.eta(:,Contt);
    psi = CA.Duals.psi(:,Contt);
    Errors = CA.Duals.Errors;
    ErrorsY = CA.Duals.ErrorsY;
    DControlError = Errors;

    A = CA.ABMH.A;
    B = CA.ABMH.B;
    M = CA.ABMH.M;
    H = CA.ABMH.H;
    s = MATDSS.Cont.s;



    alpha_xi = DER(i).ax.*alpha;
    rp_xi = DER(i).cx.*rp;

    P0Set = CA.P0Set(Contt); % P setpoint at time t
    Q0Set = CA.Q0Set(Contt);
    if strcmpi(MATDSS.Cont.CA(CAIndex).Type,'llc')
        DisableQTrackingFlag = 0;
    else
        DisableQTrackingFlag = 1;
    end


    if size(A,2) < ix+1
        Aix = size(A,2);
    else
        Aix = ix:ix+1;
    end
    if size(B,2) < ix+1
        Bix = size(B,2);
    else
        Bix = ix:ix+1;
    end


    % Check P Controller settings!
    if MATDSS.Cont.PID.EnablePFlag && at > 0.2
        Kp = MATDSS.Cont.PID.Kp;
        kappa_p_lambda = MATDSS.Cont.PID.kappa_p.lambda;
        kappa_p_mu = MATDSS.Cont.PID.kappa_p.mu;
        kappa_p_eta = MATDSS.Cont.PID.kappa_p.eta;
        kappa_p_psi = MATDSS.Cont.PID.kappa_p.psi;
        kappa_p_gamma = MATDSS.Cont.PID.kappa_p.gamma;
        kappa_p_nu = MATDSS.Cont.PID.kappa_p.nu;
        kappa_p_zeta = MATDSS.Cont.PID.kappa_p.zeta;
    else
        Kp = 0;
        kappa_p_lambda = 0;
        kappa_p_mu = 0;
        kappa_p_eta = 0;
        kappa_p_psi = 0;
        kappa_p_gamma = 0;
        kappa_p_nu = 0;
        kappa_p_zeta = 0;
    end
    % Check D Controller settings!
    if MATDSS.Cont.PID.EnableDFlag && at > 0.2
        Kd = MATDSS.Cont.PID.Kd;
        kappa_d_lambda = MATDSS.Cont.PID.kappa_d.lambda;
        kappa_d_mu = MATDSS.Cont.PID.kappa_d.mu;
        kappa_d_eta = MATDSS.Cont.PID.kappa_d.eta;
        kappa_d_psi = MATDSS.Cont.PID.kappa_d.psi;
        kappa_d_gamma = MATDSS.Cont.PID.kappa_d.gamma;
        kappa_d_nu = MATDSS.Cont.PID.kappa_d.nu;
        kappa_d_zeta = MATDSS.Cont.PID.kappa_d.zeta;

        switch lower(MATDSS.Cont.PID.Dsignal)
            case 'error'
                DControlError = Errors;
            case 'y'
                DControlError = ErrorsY;
            otherwise
                DControlError = Errors;
        end
    else
        Kd = 0;
        kappa_d_lambda = 0;
        kappa_d_mu = 0;
        kappa_d_eta = 0;
        kappa_d_psi = 0;
        kappa_d_gamma = 0;
        kappa_d_nu = 0;
        kappa_d_zeta = 0;
    end





    Ts = MATDSS.Time.Cont.TimeStep;
    xnew = x;



    % Checking if PControl or Dcontrol is enabled for this DER/VDER

    PEnable = 0;
    DEnable = 0;
    switch true
        case strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'DERs & VDERs') && strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'DERs & VDERs')
            PEnable = 1;
            DEnable = 1;
        case strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'VDERs only') && strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'DERs & VDERs')
            DEnable = 1;
            if strcmpi(DER(i).Type, 'vder')
                PEnable = 1;
            end
        case strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'VDERs only') && strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'DERs & VDERs')
            PEnable = 1;
            if strcmpi(DER(i).Type, 'vder')
                DEnable = 1;
            end
        case strcmpi(MATDSS.Cont.PID.DControlAppliedto, 'VDERs only') && strcmpi(MATDSS.Cont.PID.PControlAppliedto, 'VDERs only')
            if strcmpi(DER(i).Type, 'vder')
                DEnable = 1;
                PEnable = 1;
            end
        otherwise
            PEnable = 0;
            DEnable = 0;
    end


    PIDrho = rho;
    PIDsigma = sigma;

    PIDlambda = PEnable*Kp*kappa_p_lambda*Errors.lambda(Contt)    + lambda    + DEnable*Kd*kappa_d_lambda*(DControlError.lambda(Contt)-DControlError.lambda(Contt-1))./Ts;
    PIDmu = PEnable*Kp*kappa_p_mu*Errors.mu(Contt)                + mu        + DEnable*Kd*kappa_d_mu*(DControlError.mu(Contt)-DControlError.mu(Contt-1))./Ts;
    PIDeta = PEnable*Kp*kappa_p_eta*Errors.eta(Contt)             + eta       + DEnable*Kd*kappa_d_eta*(DControlError.eta(Contt)-DControlError.eta(Contt-1))./Ts;
    PIDpsi = PEnable*Kp*kappa_p_psi*Errors.psi(Contt)             + psi       + DEnable*Kd*kappa_d_psi*(DControlError.psi(Contt)-DControlError.psi(Contt-1))./Ts;
    PIDgamma = PEnable*Kp*kappa_p_gamma*Errors.gamma(:,Contt)     + gamma     + DEnable*Kd*kappa_d_gamma*(DControlError.gamma(:,Contt)-DControlError.gamma(:,Contt-1))./Ts;
    PIDnu = PEnable*Kp*kappa_p_nu*Errors.nu(:,Contt)              + nu        + DEnable*Kd*kappa_d_nu*(DControlError.nu(:,Contt)-DControlError.nu(:,Contt-1))./Ts;
    PIDzeta = PEnable*Kp*kappa_p_zeta*Errors.zeta(:,Contt)        + zeta      + DEnable*Kd*kappa_d_zeta*(DControlError.zeta(:,Contt)-DControlError.zeta(:,Contt-1))./Ts;



    temp2 = [0;0];
    if (strcmpi(DER(i).Type, 'vder') && mod(round(at,4),round(MATDSS.Time.Cont.TimeStep*MATDSS.Cont.RoU,4)) == 0) || ~strcmpi(DER(i).Type, 'vder')
        LoopSize = MATDSS.Cont.LoopSize;
        for k = 1:LoopSize
            xnew = xnew - [DER(i).W(1,2); DER(i).Var(1,2)];
            xnew = xnew - alpha_xi.*([f_Grad_Px;f_Grad_Qx].* xnew+ ...
                (M(:,ix:ix+1).'*(s.*(1e-0*(PIDlambda - PIDmu))*ones(size(rho)) + PIDrho*sign(-P0Set))) + ...
                (1-DisableQTrackingFlag).*(H(:,ix:ix+1).'*(s.*(1e-0.*(PIDeta-PIDpsi))*ones(size(sigma)) + PIDsigma*sign(-Q0Set))) + ...
                (A(:,Aix).'*(1e-0.*(PIDgamma - PIDnu))) +...
                (-B(:,Bix).'*(1e-0.*PIDzeta)) +...
                rp_xi.*xnew);
            %     x = xnew;
            temp2 = [temp2,xnew];
            xnew = xnew + [DER(i).W(1,2); DER(i).Var(1,2)];
            if xnew(1) > DER(i).Pmax
                xnew(1) = sign(xnew(1))*DER(i).Pmax;
            elseif xnew(1) < DER(i).Pmin
                xnew(1) = sign(xnew(1))*DER(i).Pmin;
            end

            if xnew(2) > DER(i).Qmax
                xnew(2) = sign(xnew(2))*DER(i).Qmax;
            elseif xnew(2) < DER(i).Qmin
                xnew(2) = sign(xnew(2))*DER(i).Qmin;
            end

        end
        temp2;
    end

    %{
Original code for the loop
if (strcmpi(DER(i).Type, 'vder') && mod(round(at,4),round(MATDSS.Time.Cont.TimeStep*MATDSS.Cont.RoU,4)) == 0) || ~strcmpi(DER(i).Type, 'vder')
        LoopSize = MATDSS.Cont.LoopSize;
        for k = 1:LoopSize
            xnew = xnew - [DER(i).W(1,2); DER(i).Var(1,2)];
            xnew = xnew - alpha_xi.*([f_Grad_Px;f_Grad_Qx].* xnew+ ...
                (M(:,ix:ix+1).'*(s.*(1e-0*(lambda - mu))*ones(size(rho)) + rho*sign(-P0Set))) + ...
                (1-DisableQTrackingFlag).*(H(:,ix:ix+1).'*(s.*(1e-0.*(eta-psi))*ones(size(sigma)) + sigma*sign(-Q0Set))) + ...
                (A(:,Aix).'*(1e-0.*(gamma - nu))) +...
                (-B(:,Bix).'*(1e-0.*zeta)) +...
                rp_xi.*xnew);
            %     x = xnew;
            temp2 = [temp2,xnew];
            xnew = xnew + [DER(i).W(1,2); DER(i).Var(1,2)];
            if xnew(1) > DER(i).Pmax
                xnew(1) = sign(xnew(1))*DER(i).Pmax;
            elseif xnew(1) < DER(i).Pmin
                xnew(1) = sign(xnew(1))*DER(i).Pmin;
            end

            if xnew(2) > DER(i).Qmax
                xnew(2) = sign(xnew(2))*DER(i).Qmax;
            elseif xnew(2) < DER(i).Qmin
                xnew(2) = sign(xnew(2))*DER(i).Qmin;
            end

        end
        temp2;
    end
    %}
    i_setpoint = i_setpoint + 1;
    DER(i).i_setpoint = i_setpoint;
    if ~MATDSS.Cont.DisableControllersFlag && (strcmpi(DER(i).Type, 'vder') || ~strcmpi(DER(i).Type, 'vder'))
        if strcmpi(DER(i).Type, 'vder')
            DER(i).WUnfilteredSetpoint(i_setpoint) = xnew(1);
            DER(i).VarUnvfilteredSetpoint(i_setpoint) = xnew(2);


            % Check if LPF is enabled!
            if MATDSS.Cont.LPF.EnableLPFFlag && at > 0.2
                y1 = [DER(i).WSetpoint(i_setpoint-1,2);DER(i).VarSetpoint(i_setpoint-1,2)]; % previous output power
                x1 = [DER(i).WUnfilteredSetpoint(i_setpoint-1); DER(i).VarUnvfilteredSetpoint(i_setpoint-1)]; % previous setpoint
                x0 = [DER(i).WUnfilteredSetpoint(i_setpoint); DER(i).VarUnvfilteredSetpoint(i_setpoint)]; % current setpoint

                d1d0 = MATDSS.Cont.LPF.d1d0;
                n0d0 = MATDSS.Cont.LPF.n0d0;
                n1d0 = MATDSS.Cont.LPF.n1d0;


                y0 = y1 .* -d1d0 + x0.*n0d0 + x1.*n1d0;

                xnew = y0;
            end
        end
        DER(i).WSetpoint(i_setpoint,:) = [at, xnew(1)];
        DER(i).VarSetpoint(i_setpoint,:) = [at, xnew(2)];
    else % VDER with RoU cycle not complete
        DER(i).WSetpoint(i_setpoint,:) = [at,DER(i).WSetpoint(i_setpoint-1,2)];
        DER(i).VarSetpoint(i_setpoint,:) = [at,DER(i).VarSetpoint(i_setpoint-1,2)];
    end
end

if ~strcmpi(DER(i).Type,'vder') % DER Dynamics
    my_circut = MATDSS.Sim.DSSCircuit; % handle of DSS_Circuit
    myloads = my_circut.Loads; % handle of loads in DSS_Circuit
    myloads.Name = DER(i).DSSName; %Select the DER to be updated

    % apply DER Dynamics
    x_set = [DER(i).WSetpoint(i_setpoint,2); DER(i).VarSetpoint(i_setpoint,2)]; % new setpoints
    x0 = [DER(i).W(i_WVar,2); DER(i).Var(i_WVar,2)]; % initial point


    % Call ODE to model the change
    [ODETime, ODEx] = ode45(@(t,x) FOODE(t,x,x0, Tau, x_set), tw, x0);



    % Save the dynamics
    DER(i).x.ODETime = [DER(i).x.ODETime; ODETime];
    DER(i).x.ODEx = [DER(i).x.ODEx; ODEx];

    if ~MATDSS.Cont.DisableControllersFlag
        % Update the output powers
        DER(i).W(i_WVar+1,2) = ODEx(end,1);
        DER(i).Var(i_WVar+1,2) = ODEx(end,2);
    else
        DER(i).W(i_WVar+1,2) = DER(i).W(i_WVar,2);
        DER(i).Var(i_WVar+1,2) = DER(i).Var(i_WVar,2);
    end

    DER(i).i_WVar = i_WVar + 1;

    % Update the output power in OpenDSS
    myloads.Name = DER(i).DSSName; %Select the DER to be updated

    % Update the DER P and Q setpoints
    % [-myloads.kW, -myloads.kvar;...
    %     DER(i).kW(t+1), DER(i).kvar(t+1)];
    myloads.kW = -DER(i).W(i_WVar+1,2)./1e3;
    myloads.kvar = -DER(i).Var(i_WVar+1,2)./1e3;

else

    i_ = strfind(DER(i).DSSName,'_');
    CAIndex = str2num(DER(i).DSSName(i_(end)+3:end));
    CA = MATDSS.Cont.CA(CAIndex);

    Contt = CA.TimeIndex;
    DER(i).W(i_WVar+1,2) = [-sum(CA.P0(:,Contt))];
    DER(i).Var(i_WVar+1,2) = [-sum(CA.Q0(:,Contt))];

    DER(i).i_WVar = i_WVar + 1;

end
end

% ODE function
function dxdt = FOODE(t,x,x0, Tau, x_set)
% t and x are ode variables
% x0 is initial state
% tau is diagonal time constant for P and Q (ordered)
% x_set is a 2x1 input(set-points) [P;Q] ordered

% diagolize tau if it is not diag. matrix
if size(Tau,1) ~= size(Tau,2)
    Tau = diag(Tau);
end

dxdt = -inv(Tau)*(x-x_set);

end





%}