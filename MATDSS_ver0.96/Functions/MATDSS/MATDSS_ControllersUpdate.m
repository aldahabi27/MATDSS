function [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS,DER)
% MATDSS_ControllersUpdate(MATDSS,DER)
% This function updates the controllers in the MATDSS system, ensuring that
% their measurements and control actions are properly synchronized with the
% simulation time. It accounts for measurement delays and applies corrective
% measures based on controller settings.
%
% The function is responsible for calculating control variables such as
% active and reactive power (P0, Q0), voltage limits, and other duals like
% rho, lambda, mu, gamma, nu, and zeta. These are used to ensure the system
% operates within the desired parameters.
%
% This is a core function for real-time feedback and optimization, ensuring
% that the Distributed Energy Resources (DER) and controllers work in
% harmony, adjusting their operations based on the latest measurements and
% setpoints.
%
% Parameters:
%   - MATDSS: The main structure holding the simulation state, including
%             time, measurements, control settings, etc.
%   - DER: A structure representing the Distributed Energy Resources, which
%          interact with the controllers.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com



% Extract current time (t) and adjusted time (at) from simulation state
t = MATDSS.Sim.t;
at = MATDSS.Sim.at;

% Define a time range for the period of interest (used later)
MyPer = 0.75; MyPerEnd = 0.90;
MyPerSec = (MyPerEnd - MyPer)/1; % Calculate period section

% Check if the current time step matches the controller's time step
if mod(at, MATDSS.Time.Cont.TimeStep) == 0
    % Calculate the delay index based on measurement delay and time step
    DelayIndex = ceil((MATDSS.Time.Meas.Delay/MATDSS.Time.Meas.TimeStep));

    % Determine positions of available measurements relative to the adjusted time
    pos = MATDSS.Meas.at - at;
    pos(pos > 0) = []; % Remove future measurements
    Meast = length(pos); % Count available past measurements
    Meast = Meast - DelayIndex; % Account for delay

    % Find the controller's time index (Contt) corresponding to the current time
    Contt = find(MATDSS.Cont.at == at);

    % If there are no valid measurements, display a warning and disable the controller
    if Meast < 1
        disp('Warning, delay is high and measurements are missing. Controller is disabled!')
    else
        % Iterate over each controller in the system
        for i = 1:length(MATDSS.Cont.CA)
            CA = MATDSS.Cont.CA(i); % Get the current controller

            pause(0.001); % Small pause for performance reasons

            % Check if the controller's A, B, M, and H matrices exist; if not, calculate them
            if ~(isfield(CA.ABMH,'A') && isfield(CA.ABMH,'B') && isfield(CA.ABMH,'M') && isfield(CA.ABMH,'H'))
                % Call MATDSS_ABMH to compute the matrices
                [A, B, M, H, Mv, MvIndex, Mi, MiBuses] = MATDSS_ABMH(MATDSS,DER,CA);
                CA.ABMH.A = A; % Store A matrix
                CA.ABMH.B = B; % Store B matrix
                CA.ABMH.M = M; % Store M matrix
                CA.ABMH.H = H; % Store H matrix
            end

            % Retrieve voltage (V) and base voltage (Vbase) for the controller if k_v exists
            if isempty(CA.k_v)
                V = 0; % Default voltage to 0 if k_v is empty
                Vbase = 0; % Default base voltage to 0 if k_v is empty
            else
                V = MATDSS.Meas.VMagProfile(CA.k_vIndex, Meast); % Measured voltage magnitude
                Vbase = MATDSS.Meas.VBases(CA.k_v); % Base voltage
            end

            % Retrieve current (I) and upper limit (Iul) for the controller if k_i exists
            if isempty(CA.k_i)
                I = 0; % Default current to 0 if k_i is empty
                Iul = 1e99; % Set a large default current upper limit
            else
                I = MATDSS.Meas.IProfile(CA.k_iIndex, Meast); % Measured current
                Iul = MATDSS.Meas.Imax(CA.k_iIndex); % Maximum allowed current
            end

            % Calculate P0 and Q0 (active and reactive power) using the helper function
            [P0, Q0, P0Set, Q0Set] = MATDSS_CAP0Q0Calculation(MATDSS, DER, CA, Meast);

            % Update controller time index and store calculated power values
            CA.TimeIndex = Contt + 1;
            CA.P0Set(Contt + 1) = [P0Set]; % Store set active power
            CA.Q0Set(Contt + 1) = [Q0Set]; % Store set reactive power
            CA.P0(:,Contt + 1) = [P0]; % Store actual active power
            CA.Q0(:,Contt + 1) = [Q0]; % Store actual reactive power

            % Retrieve error values for the controller
            E = CA.E;

            % Adjust base voltage if k_v exists
            if length(CA.k_v) < 1
                Vbase = 0; % Set to 0 if no voltage base defined
            else
                Vbase = MATDSS.Meas.VBases(CA.k_v) * 1e3; % Convert base voltage to volts
            end

            % Calculate upper (Vul) and lower (Vll) voltage limits
            Vul = CA.vul * Vbase; % Upper voltage limit
            Vll = CA.vll * Vbase; % Lower voltage limit

            % Extract dual variables for optimization calculations
            rho = CA.Duals.rho(:,Contt); % Dual variable for active power control
            lambda = CA.Duals.lambda(:,Contt); % Dual variable for active power setpoint
            mu = CA.Duals.mu(:,Contt); % Dual variable for reactive power control
            gamma = CA.Duals.gamma(:,Contt); % Dual variable for upper voltage limit
            nu = CA.Duals.nu(:,Contt); % Dual variable for lower voltage limit
            zeta = CA.Duals.zeta(:,Contt); % Dual variable for current limit control

            % Depending on controller type, adjust dual variables for L2C-type controllers
            if strcmpi(CA.Type, 'l2c')
                sigma = zeros(size(CA.Duals.sigma(:,Contt))); % Initialize sigma
                eta = zeros(size(CA.Duals.eta(:,Contt))); % Initialize eta
                psi = zeros(size(CA.Duals.psi(:,Contt))); % Initialize psi
            else
                sigma = CA.Duals.sigma(:,Contt); % Use existing sigma
                eta = CA.Duals.eta(:,Contt); % Use existing eta
                psi = CA.Duals.psi(:,Contt); % Use existing psi
            end

            % Update rho (dual variable for power control) with its new value
            alpharho = CA.arho * CA.alpha; % Scaling factor for rho update (depends on alpha and specific CA parameters)
            rdrho = CA.crho * CA.rbard; % Correction term for rho based on system parameters
            rhoNew = rho + alpharho * (sign(-P0Set) * P0 - rdrho * rho);
            % Calculate new rho using a correction term based on the difference between
            % the setpoint P0Set and actual power P0, adjusted by rdrho to control the update rate
            rhoError = 0; % Placeholder for error in rho, set to 0 for now
            rhoY = 0; % Placeholder for additional tracking value (Y) for rho, initialized to 0

            % Update lambda (dual variable for balancing power mismatch)
            alphalambda = CA.alambda * CA.alpha; % Scaling factor for lambda update
            rdlambda = CA.clambda * CA.rbard; % Correction term for lambda update
            lambdaNew = lambda + alphalambda * ((sum(P0) - P0Set) - E - rdlambda * lambda);
            % Update lambda based on the difference between total power sum(P0) and the setpoint P0Set,
            % considering energy mismatch E and applying correction via rdlambda
            lambdaError = (sum(P0) - P0Set) - E - rdlambda * lambda;
            % Calculate the error for lambda, reflecting power imbalance and correction term
            % lambdaError(lambdaError < 0) = 0; % Uncomment this to prevent negative lambda errors (optional)
            lambdaY = sum(P0); % Track total power (sum of P0) for lambda

            % Update mu (dual variable for reactive power balancing)
            alphamu = CA.amu * CA.alpha; % Scaling factor for mu update
            rdmu = CA.cmu * CA.rbard; % Correction term for mu update
            muNew = mu + alphamu * ((P0Set - sum(P0)) - E - rdmu * mu);
            % Update mu based on the difference between the setpoint P0Set and the total power sum(P0),
            % accounting for energy mismatch E and adjusting by rdmu
            muError = (P0Set - sum(P0)) - E - rdmu * mu;
            % Calculate error for mu based on setpoint deviation and correction
            % muError(muError < 0) = 0; % Uncomment this to prevent negative mu errors (optional)
            muY = -sum(P0); % Track negative total power for mu (could indicate power deficit)

            % Update gamma (dual variable for voltage lower limit violation)
            alphagamma = CA.agamma * CA.alpha; % Scaling factor for gamma update
            rdgamma = CA.cgamma * CA.rbard; % Correction term for gamma update
            gammaNew = gamma + alphagamma * ((V - Vul) - rdgamma * gamma);
            % Update gamma based on the difference between actual voltage V and the lower voltage limit Vul,
            % applying correction via rdgamma
            gammaError = (V - Vul) - rdgamma * gamma;
            % Calculate error for gamma, capturing voltage violation and correction
            % gammaError(gammaError < 0) = 0; % Uncomment this to prevent negative gamma errors (optional)
            gammaY = V; % Track actual voltage V for gamma

            % Update nu (dual variable for voltage upper limit violation)
            alphanu = CA.anu * CA.alpha; % Scaling factor for nu update
            rdnu = CA.cnu * CA.rbard; % Correction term for nu update
            nuNew = nu + alphanu * ((Vll - V) - rdnu * nu);
            % Update nu based on the difference between upper voltage limit Vll and the actual voltage V,
            % applying correction via rdnu
            nuError = (Vll - V) - rdnu * nu;
            % Calculate error for nu, capturing voltage violation and applying correction
            % nuError(nuError < 0) = 0; % Uncomment this to prevent negative nu errors (optional)
            nuY = -V; % Track negative voltage for nu (to capture deviation from Vll)

            % Update zeta (dual variable for current limit violation)
            alphazeta = CA.azeta * CA.alpha; % Scaling factor for zeta update
            rdzeta = CA.czeta * CA.rbard; % Correction term for zeta update
            zetaNew = zeta + alphazeta * (I - Iul - rdzeta * zeta);
            % Update zeta based on the difference between actual current I and the upper current limit Iul,
            % applying correction via rdzeta
            zetaError = I - Iul - rdzeta * zeta;
            % Calculate error for zeta, capturing current limit violation and applying correction
            % zetaError(zetaError < 0) = 0; % Uncomment this to prevent negative zeta errors (optional)
            zetaY = I; % Track actual current I for zeta

            % Check if the control algorithm type is 'l2c' (specific control type)
            if strcmpi(CA.Type,'l2c')
                % For 'l2c', directly use the old values without updating
                sigmaNew = sigma; sigmaError = 0; sigmaY = 0; % No update, use existing sigma
                etaNew = eta; etaError = 0; etaY = 0; % No update, use existing eta
                psiNew = psi; psiError = 0; psiY = 0; % No update, use existing psi
            else
                % For other types, update sigma, eta, and psi with new values
                alphasigma = CA.asigma * CA.alpha; % Scaling factor for sigma update
                rdsigma = CA.csigma * CA.rbard; % Correction term for sigma update
                sigmaNew = sigma + alphasigma * (sign(-Q0Set) * Q0 - rdsigma * sigma);
                % Update sigma based on the difference between Q0 and its setpoint, applying correction
                sigmaError = 0; % Placeholder for sigma error, set to 0 for now
                sigmaY = 0; % Placeholder for tracking value of sigma, initialized to 0

                % Update eta (dual variable for reactive power balancing)
                alphaeta = CA.aeta * CA.alpha; % Scaling factor for eta update
                rdeta = CA.ceta * CA.rbard; % Correction term for eta update
                etaNew = eta + alphaeta * ((sum(Q0) - Q0Set) - E - rdeta * eta);
                % Update eta based on the difference between sum of reactive power Q0 and its setpoint,
                % adjusted by rdeta and error term E
                etaError = (sum(Q0) - Q0Set) - E - rdeta * eta;
                % Calculate eta error based on reactive power mismatch and correction
                etaY = sum(Q0); % Track total reactive power Q0 for eta

                % Update psi (dual variable for reactive power balancing)
                alphapsi = CA.apsi * CA.alpha; % Scaling factor for psi update
                rdpsi = CA.cpsi * CA.rbard; % Correction term for psi update
                psiNew = psi + alphapsi * ((Q0Set - sum(Q0)) - E - rdpsi * psi);
                % Update psi based on the difference between reactive power setpoint Q0Set and sum(Q0),
                % adjusted by rdpsi and error term E
                psiError = (Q0Set - sum(Q0)) - E - rdpsi * psi;
                % Calculate psi error based on reactive power setpoint deviation
                psiY = -sum(Q0); % Track negative sum of Q0 for psi
            end

            % Fix projection to R+ (Ensure all dual variables are non-negative)
            rhoNew(rhoNew < 0) = 0; % Clamp rhoNew to non-negative values
            lambdaNew(lambdaNew < 0) = 0; % Clamp lambdaNew to non-negative values
            muNew(muNew < 0) = 0; % Clamp muNew to non-negative values
            gammaNew(gammaNew < 0) = 0; % Clamp gammaNew to non-negative values
            nuNew(nuNew < 0) = 0; % Clamp nuNew to non-negative values
            zetaNew(zetaNew < 0) = 0; % Clamp zetaNew to non-negative values

            sigmaNew(sigmaNew < 0) = 0; % Clamp sigmaNew to non-negative values
            etaNew(etaNew < 0) = 0; % Clamp etaNew to non-negative values
            psiNew(psiNew < 0) = 0; % Clamp psiNew to non-negative values

            % Save the updated dual variables for the current time step (Contt+1)
            CA.Duals.rho(:,Contt+1) = [rhoNew]; % Save new rho values
            CA.Duals.lambda(:,Contt+1) = [lambdaNew]; % Save new lambda values
            CA.Duals.mu(:,Contt+1) = [muNew]; % Save new mu values
            CA.Duals.gamma(:,Contt+1) = [gammaNew]; % Save new gamma values
            CA.Duals.nu(:,Contt+1) = [nuNew]; % Save new nu values
            CA.Duals.zeta(:,Contt+1) = [zetaNew]; % Save new zeta values

            CA.Duals.sigma(:,Contt+1) = [sigmaNew]; % Save new sigma values
            CA.Duals.eta(:,Contt+1) = [etaNew]; % Save new eta values
            CA.Duals.psi(:,Contt+1) = [psiNew]; % Save new psi values

            % Save the computed errors for dual variables
            CA.Duals.Errors.rho(:,Contt+1) = [rhoError]; % Save error in rho
            CA.Duals.Errors.lambda(:,Contt+1) = [lambdaError]; % Save error in lambda
            CA.Duals.Errors.mu(:,Contt+1) = [muError]; % Save error in mu
            CA.Duals.Errors.gamma(:,Contt+1) = [gammaError]; % Save error in gamma
            CA.Duals.Errors.nu(:,Contt+1) = [nuError]; % Save error in nu
            CA.Duals.Errors.zeta(:,Contt+1) = [zetaError]; % Save error in zeta

            CA.Duals.Errors.sigma(:,Contt+1) = [sigmaError]; % Save error in sigma
            CA.Duals.Errors.eta(:,Contt+1) = [etaError]; % Save error in eta
            CA.Duals.Errors.psi(:,Contt+1) = [psiError]; % Save error in psi

            % Save the auxiliary error tracking variables (Y)
            CA.Duals.ErrorsY.rho(:,Contt+1) = [rhoY]; % Save Y value for rho
            CA.Duals.ErrorsY.lambda(:,Contt+1) = [lambdaY]; % Save Y value for lambda
            CA.Duals.ErrorsY.mu(:,Contt+1) = [muY]; % Save Y value for mu
            CA.Duals.ErrorsY.gamma(:,Contt+1) = [gammaY]; % Save Y value for gamma
            CA.Duals.ErrorsY.nu(:,Contt+1) = [nuY]; % Save Y value for nu
            CA.Duals.ErrorsY.zeta(:,Contt+1) = [zetaY]; % Save Y value for zeta

            CA.Duals.ErrorsY.sigma(:,Contt+1) = [sigmaY]; % Save Y value for sigma
            CA.Duals.ErrorsY.eta(:,Contt+1) = [etaY]; % Save Y value for eta
            CA.Duals.ErrorsY.psi(:,Contt+1) = [psiY]; % Save Y value for psi

            % Loop over all DERs to update their status and setpoints
            for j = 1:CA.nDER
                DER(CA.DERIndex(j)).UpdateStatus = true; % Mark DERs for setpoint update
            end

            % Store the updated CA (control algorithm) structure back to MATDSS
            MATDSS.Cont.CA(i) = CA;

        end % End of for loop (controlling CA systems)

    end % End of if Meast < 1 - else

end % End of if mod (at) --> controller update timestep

end % End of the function

%% Old function Code
%{

function [MATDSS, DER] = MATDSS_ControllersUpdate(MATDSS,DER)

t = MATDSS.Sim.t;
at = MATDSS.Sim.at;
MyPer = 0.75; MyPerEnd = 0.90;
MyPerSec = (MyPerEnd - MyPer)/1;
if mod(at,MATDSS.Time.Cont.TimeStep) == 0
    DelayIndex = ceil((MATDSS.Time.Meas.Delay/MATDSS.Time.Meas.TimeStep));
    pos = MATDSS.Meas.at - at;
    pos(pos > 0) = [];
    Meast = length(pos);
    Meast = Meast - DelayIndex;
    Contt = find(MATDSS.Cont.at == at);
    if Meast < 1
        disp('Warning, delay is high and measurements are missing. Controller is disabled!')
    else
        for i = 1:length(MATDSS.Cont.CA)
            CA = MATDSS.Cont.CA(i);
            pause(0.001);
            if ~(isfield(CA.ABMH,'A') && isfield(CA.ABMH,'B') && isfield(CA.ABMH,'M') && isfield(CA.ABMH,'H'))
                [A, B, M, H, Mv, MvIndex, Mi, MiBuses] = MATDSS_ABMH(MATDSS,DER,CA);
                CA.ABMH.A = A;
                CA.ABMH.B = B;
                CA.ABMH.M = M;
                CA.ABMH.H = H;
            end

            if isempty(CA.k_v)
                V = 0;
                Vbase = 0;
            else
                V = MATDSS.Meas.VMagProfile(CA.k_vIndex,Meast);
                Vbase = MATDSS.Meas.VBases(CA.k_v);
            end

            if isempty(CA.k_i)
                I = 0;
                Iul = 1e99;
            else
                I = MATDSS.Meas.IProfile(CA.k_iIndex,Meast);
                Iul = MATDSS.Meas.Imax(CA.k_iIndex);
            end


            [P0, Q0, P0Set, Q0Set] = MATDSS_CAP0Q0Calculation(MATDSS,DER,CA,Meast);
            CA.TimeIndex = Contt+1;
            CA.P0Set(Contt+1) = [P0Set];
            CA.Q0Set(Contt+1) = [Q0Set];
            CA.P0(:,Contt+1) = [P0];
            CA.Q0(:,Contt+1) = [Q0];

            E = CA.E;
            if length(CA.k_v) < 1
                Vbase = 0;
            else
                Vbase = MATDSS.Meas.VBases(CA.k_v).*1e3; %Vbase in V
            end
            Vul = CA.vul.*Vbase;
            Vll = CA.vll.*Vbase;

            rho = CA.Duals.rho(:,Contt);
            lambda = CA.Duals.lambda(:,Contt);
            mu = CA.Duals.mu(:,Contt);
            gamma = CA.Duals.gamma(:,Contt);
            nu = CA.Duals.nu(:,Contt);
            zeta = CA.Duals.zeta(:,Contt);



            if strcmpi(CA.Type,'l2c')
                sigma = zeros(size(CA.Duals.sigma(:,Contt)));
                eta = zeros(size(CA.Duals.eta(:,Contt)));
                psi = zeros(size(CA.Duals.psi(:,Contt)));
            else
                sigma = CA.Duals.sigma(:,Contt);
                eta = CA.Duals.eta(:,Contt);
                psi = CA.Duals.psi(:,Contt);
            end


            alpharho = CA.arho.*CA.alpha;
            rdrho = CA.crho.*CA.rbard;
            rhoNew = rho + alpharho.*(sign(-P0Set)*P0 - rdrho.*rho);
            rhoError = 0;%sign(-P0Set)*P0 - rdrho.*rho;
            rhoY = 0;


            alphalambda = CA.alambda.*CA.alpha;
            rdlambda = CA.clambda.*CA.rbard;
            lambdaNew = lambda + alphalambda.*((sum(P0) - P0Set) - E - rdlambda.*lambda);
            lambdaError = (sum(P0) - P0Set) - E - rdlambda.*lambda;
            % lambdaError(lambdaError < 0) = 0;
            lambdaY = sum(P0);

            alphamu = CA.amu.*CA.alpha;
            rdmu = CA.cmu.*CA.rbard;
            muNew = mu + alphamu.*((P0Set - sum(P0)) - E - rdmu.*mu);
            muError = (P0Set - sum(P0)) - E - rdmu.*mu;
            % muError(muError < 0) = 0;
            muY = - sum(P0);

            alphagamma = CA.agamma.*CA.alpha;
            rdgamma = CA.cgamma.*CA.rbard;
            gammaNew = gamma + alphagamma.*((V - Vul) - rdgamma.*gamma);
            gammaError = (V - Vul) - rdgamma.*gamma;
            % gammaError(gammaError < 0) = 0;
            gammaY = V;

            alphanu = CA.anu.*CA.alpha;
            rdnu = CA.cnu.*CA.rbard;
            nuNew = nu + alphanu.*((Vll - V) - rdnu.*nu);
            nuError = (Vll - V) - rdnu.*nu;
            % nuError(nuError < 0) = 0;
            nuY = -V;

            alphazeta = CA.azeta.*CA.alpha;
            rdzeta = CA.czeta.*CA.rbard;
            zetaNew = zeta + alphazeta.*(I - Iul - rdzeta.*zeta);
            zetaError = I - Iul - rdzeta.*zeta;
            % zetaError(zetaError < 0) = 0;
            zetaY = I;

            if strcmpi(CA.Type,'l2c')
                sigmaNew = sigma; sigmaError = 0; sigmaY = 0;
                etaNew = eta; etaError = 0; etaY = 0;
                psiNew = psi; psiError = 0; psiY = 0;
            else

                alphasigma = CA.asigma.*CA.alpha;
                rdsigma = CA.csigma.*CA.rbard;
                sigmaNew = sigma + alphasigma.*(sign(-Q0Set)*Q0 - rdsigma.*sigma);
                sigmaError = 0; %sign(-Q0Set)*Q0 - rdsigma.*sigma
                sigmaY = 0;

                alphaeta = CA.aeta.*CA.alpha;
                rdeta = CA.ceta.*CA.rbard;
                etaNew = eta + alphaeta.*((sum(Q0) - Q0Set) - E - rdeta.*eta);
                etaError = (sum(Q0) - Q0Set) - E - rdeta.*eta;
                % etaError(etaError < 0) = 0;
                etaY = sum(Q0);


                alphapsi = CA.apsi.*CA.alpha;
                rdpsi = CA.cpsi.*CA.rbard;
                psiNew = psi + alphapsi.*((Q0Set - sum(Q0)) - E - rdpsi.*psi);
                psiError = (Q0Set - sum(Q0)) - E - rdpsi.*psi;
                % psiError(psiError < 0) = 0;
                psiY = -sum(Q0);
            end

            % Fix projection to R+
            rhoNew(rhoNew < 0) = 0;
            lambdaNew(lambdaNew < 0) = 0;
            muNew(muNew < 0) = 0;
            gammaNew(gammaNew < 0) = 0;
            nuNew(nuNew < 0) = 0;
            zetaNew(zetaNew < 0) = 0;

            sigmaNew(sigmaNew < 0) = 0;
            etaNew(etaNew < 0) = 0;
            psiNew(psiNew < 0) = 0;


            CA.Duals.rho(:,Contt+1) = [rhoNew];
            CA.Duals.lambda(:,Contt+1) = [lambdaNew];
            CA.Duals.mu(:,Contt+1) = [muNew];
            CA.Duals.gamma(:,Contt+1) = [gammaNew];
            CA.Duals.nu(:,Contt+1) = [nuNew];
            CA.Duals.zeta(:,Contt+1) = [zetaNew];

            CA.Duals.sigma(:,Contt+1) = [sigmaNew];
            CA.Duals.eta(:,Contt+1) = [etaNew];
            CA.Duals.psi(:,Contt+1) = [psiNew];




            % Saving the errors
            CA.Duals.Errors.rho(:,Contt+1) = [rhoError];
            CA.Duals.Errors.lambda(:,Contt+1) = [lambdaError];
            CA.Duals.Errors.mu(:,Contt+1) = [muError];
            CA.Duals.Errors.gamma(:,Contt+1) = [gammaError];
            CA.Duals.Errors.nu(:,Contt+1) = [nuError];
            CA.Duals.Errors.zeta(:,Contt+1) = [zetaError];

            CA.Duals.Errors.sigma(:,Contt+1) = [sigmaError];
            CA.Duals.Errors.eta(:,Contt+1) = [etaError];
            CA.Duals.Errors.psi(:,Contt+1) = [psiError];

            % Saving the errors
            CA.Duals.ErrorsY.rho(:,Contt+1) = [rhoY];
            CA.Duals.ErrorsY.lambda(:,Contt+1) = [lambdaY];
            CA.Duals.ErrorsY.mu(:,Contt+1) = [muY];
            CA.Duals.ErrorsY.gamma(:,Contt+1) = [gammaY];
            CA.Duals.ErrorsY.nu(:,Contt+1) = [nuY];
            CA.Duals.ErrorsY.zeta(:,Contt+1) = [zetaY];

            CA.Duals.ErrorsY.sigma(:,Contt+1) = [sigmaY];
            CA.Duals.ErrorsY.eta(:,Contt+1) = [etaY];
            CA.Duals.ErrorsY.psi(:,Contt+1) = [psiY];
            for j = 1:CA.nDER
                DER(CA.DERIndex(j)).UpdateStatus = true; % Tell DER to update the setpoint
            end

            MATDSS.Cont.CA(i) = CA;

        end

    end

end

end


%}