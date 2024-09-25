function [MATDSS, DER]= MATDSS_L2C(app,MATDSS, DER)
% MATDSS = Level2Controller_Init(MATDSS,DER) is function that will
% initiallize the Level2Controller with all necessary parameters and
% initial values of "zero". This function should be called only once before
% the "time" analysis has started (t = 0).

%{
%% Obtaining A, B, M and H now

A = QLv*[KWye KDelta]*QR;  % Unit is [unitless] * [1/A] * [unitless]
% B = B;  % Unit is [1/V] * [unitless] = 1/V
M = real(GBar);  % [unitless] * [unitless] = unitless
H = imag(GBar);  % [unitless] * [unitless] = unitless
%}
t = MATDSS.Sim.t;
at = MATDSS.Sim.at;

if mod(at,MATDSS.Time.L2C.TimeStep) == 0
    DelayIndex = ceil((MATDSS.Time.Meas.Delay/MATDSS.Time.Meas.TimeStep));
    [~, Meast] = min(abs(MATDSS.Meas.at - at));
    Meast = Meast - DelayIndex;
    MATDSS.L2C.at = at;


    if Meast < 1
        disp('Warning, delay is high and measurements are missing. Controller is disabled!')
    else
        nDER = MATDSS.Sim.nDER;
        

        if ~isfield(MATDSS.L2C,'alpha')
           
            [A, B, M, H, Mv, MvIndex, Mi, MiIndex, MiIndexRef,MiBuses] = MATDSS_ABMH(app,MATDSS,DER); % Obtain A, B, M and H matrices
            
            
            L2CTableData = MATDSS.TableData.L2C;
            if MATDSS_StrComp(L2CTableData(1), 'auto')
                MATDSS.L2C.alpha = MATDSS.L2C.Gain*MATDSS.Time.L2C.TimeStep;% MATDSS.Time.Sim.TimeStep./MATDSS.Time.L2C.TimeStep; %alpha = 0.2
            else
                MATDSS.L2C.alpha = str2double(L2CTableData(1));
            end

            MATDSS.L2C.rp = str2double(L2CTableData(2));
            MATDSS.L2C.rbard = str2double(L2CTableData(3));
            MATDSS.L2C.E = str2double(L2CTableData(4));
            MATDSS.L2C.Vul = str2double(L2CTableData(5)); %upper limit in p.u.
            MATDSS.L2C.Vll = str2double(L2CTableData(6)); %lower limit in p.u.
            MATDSS.L2C.Iul = MATDSS.Meas.Imax;
            MATDSS.L2C.Iul(MATDSS.L2C.Iul <=0) = 1e99; % any unspecidfied limit is set to inf.

            MATDSS.L2C.alphaDual = str2double(convertCharsToStrings(L2CTableData(8:13)'))'.*MATDSS.L2C.alpha; % All dual variables' alphas
            MATDSS.L2C.rdDual = str2double(convertCharsToStrings(L2CTableData(14:19)'))'.*MATDSS.L2C.rbard; % All dual variables' alphas
            



            % Define the dual variables and initiallize them
            MATDSS.L2C.rho(:,1) = zeros(3,1);
            MATDSS.L2C.lambda(:,1) = zeros(1,1);
            MATDSS.L2C.mu(:,1) = zeros(1,1);
            MATDSS.L2C.gamma(:,1) = zeros(max(size(Mv)),1);
            MATDSS.L2C.nu(:,1) = zeros(max(size(Mv)),1);
            MATDSS.L2C.zeta(:,1) = zeros(sum(MiIndex(:,2)),1);

            MATDSS.L2C.ABMH.A = A;
            MATDSS.L2C.ABMH.B = B;
            MATDSS.L2C.ABMH.M = M;
            MATDSS.L2C.ABMH.H = H;
            MATDSS.L2C.ABMH.Mv = Mv;
            MATDSS.L2C.ABMH.MvIndex = MvIndex;
            MATDSS.L2C.ABMH.Mi = Mi;
            MATDSS.L2C.ABMH.MiIndex = MiIndex;
            MATDSS.L2C.ABMH.MiIndexRef = MiIndexRef;
            MATDSS.L2C.ABMH.MiBuses = MiBuses;
            
        else
            % Level 2 Controller algorithm
            % step 1 - Updating the dual variables
            %             t = at/MATDSS.Time.L2C.TimxeStep + 1;
            %             Meast = MATDSS.Time.Meas.QueueIndex; % index for measurements to be used
            A = MATDSS.L2C.ABMH.A;
            B = MATDSS.L2C.ABMH.B;
            M = MATDSS.L2C.ABMH.M;
            H = MATDSS.L2C.ABMH.H;
            Mv = MATDSS.L2C.ABMH.Mv;
            MvIndex = MATDSS.L2C.ABMH.MvIndex;
            Mi = MATDSS.L2C.ABMH.Mi;
            MiIndex = MATDSS.L2C.ABMH.MiIndex;
            MiIndexRef = MATDSS.L2C.ABMH.MiIndexRef;
            
            [A, B, M, H, Mv, MvIndex, Mi, MiIndex, MiIndexRef,MiBuses] = MATDSS_ABMH(app,MATDSS,DER); % Obtain A, B, M and H matrices


            
            if MvIndex < 1
                Vbase = 0;
            else
                Vbase = MATDSS.Meas.VBases(MvIndex).*1e3; %Vbase in V
            end

            P0 = MATDSS.Meas.P0(:,Meast);
            V = MATDSS.Meas.VMagProfile(:,Meast); % we need to make sure we put here only the voltages we have considered in our matrix A
            I = MATDSS.Meas.IProfile(:,Meast);
            P0Set = MATDSS.ControlSignals.P0Set(t); % P setpoint at time t
            alphaDual = MATDSS.L2C.alphaDual;
            rdDual = MATDSS.L2C.rdDual;
            E = MATDSS.L2C.E;
            Vul = MATDSS.L2C.Vul.*Vbase;
            Vll = MATDSS.L2C.Vll.*Vbase;
            Iul = MATDSS.L2C.Iul;


            rho = MATDSS.L2C.rho(:,end);
            lambda = MATDSS.L2C.lambda(:,end);
            mu = MATDSS.L2C.mu(:,end);
            gamma = MATDSS.L2C.gamma(:,end);
            nu = MATDSS.L2C.nu(:,end);
            zeta = MATDSS.L2C.zeta(:,end);



            % update the dual variable here
            rhoNew      = rho      + alphaDual(1).*(sign(-P0Set)*sum(P0) - rdDual(1).*rho);
            lambdaNew   = lambda   + alphaDual(2).*((sum(P0) - P0Set) - E - rdDual(2).*lambda);
            muNew       = mu       + alphaDual(3).*((P0Set - sum(P0)) - E - rdDual(3).*mu);
            gammaNew    = gamma    + alphaDual(4).*((V - Vul) - rdDual(4).*gamma);
            nuNew       = nu       + alphaDual(5).*((Vll - V) - rdDual(5).*nu);
            zetaNew     = zeta     + alphaDual(6).*(I - Iul - rdDual(6).*zeta);

            % Fix projection to R+
            rhoNew(rhoNew < 0) = 0;
            lambdaNew(lambdaNew < 0) = 0;
            muNew(muNew < 0) = 0;
            gammaNew(gammaNew < 0) = 0;
            nuNew(nuNew < 0) = 0;
            zetaNew(zetaNew < 0) = 0;

            %% Debuging L2C
            
%             rhoNew (:) = 0;

            %%
            %     debugvar = []
            MATDSS.L2C.rho = [MATDSS.L2C.rho, rhoNew];
            MATDSS.L2C.lambda = [MATDSS.L2C.lambda, lambdaNew];
            MATDSS.L2C.mu = [MATDSS.L2C.mu, muNew];
            MATDSS.L2C.gamma = [MATDSS.L2C.gamma, gammaNew];
            MATDSS.L2C.nu = [MATDSS.L2C.nu, nuNew];
            MATDSS.L2C.zeta = [MATDSS.L2C.zeta, zetaNew];
            
            for i = 1:MATDSS.Sim.nDER
                DER(i).UpdateStatus = true; % Tell DER to update the setpoint
            end
            
        end
    end
end