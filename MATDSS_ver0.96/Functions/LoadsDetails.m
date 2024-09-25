% Loop over loads and get their status

myloads = MATDSS.Sim.DSSCircuit.Loads;
myloadsnames = myloads.AllNames;

MyLoadsSummary = cell(length(myloadsnames),5);
i_load = myloads.First;

while i_load
    MyLoadsSummary(i_load, :) = {myloads.Name, myloads.kw, myloads.kvar, myloads.PF, myloads.kva};
    i_load = myloads.Next;
end


%%
clear MyLoadsSummary_t5t0
MyLoadsSummary_t5t0 = MyLoadsSummary_t0(:,1);
MyLoadsSummary_t5t0 = [MyLoadsSummary_t5t0, mat2cell([[MyLoadsSummary_t5{:,2}]-[MyLoadsSummary_t0{:,2}]]',ones(length(myloadsnames),1))]
MyLoadsSummary_t5t0 = [MyLoadsSummary_t5t0, mat2cell([[MyLoadsSummary_t5{:,3}]-[MyLoadsSummary_t0{:,3}]]',ones(length(myloadsnames),1))]
MyLoadsSummary_t5t0 = [MyLoadsSummary_t5t0, mat2cell([[MyLoadsSummary_t5{:,4}]-[MyLoadsSummary_t0{:,4}]]',ones(length(myloadsnames),1))]
MyLoadsSummary_t5t0 = [MyLoadsSummary_t5t0, mat2cell([[MyLoadsSummary_t5{:,5}]-[MyLoadsSummary_t0{:,5}]]',ones(length(myloadsnames),1))]

%%
DERs_kw = [MyLoadsSummary_t5t0{1178:end,2}]';
DERs_kvar = [MyLoadsSummary_t5t0{1178:end,3}]';
DERs_kva = [MyLoadsSummary_t5t0{1178:end,5}]';


%%
sum(DERs_kw)
sum(DERs_kvar)
sum(DERs_kva)

%% Check the change in voltage profile

Vt0 = MATDSS.Sim.Meas.VMagProfilePu(:,1);
Vt5 = MATDSS.Sim.Meas.VMagProfilePu(:,t-1);

Vt5t0 = Vt5-Vt0;



%% check the change in power profile!

Pt0 = MATDSS.Sim.Meas.PProfile(:,1);
Pt5 = MATDSS.Sim.Meas.PProfile(:,t-1);

Pt5t0 = Pt5-Pt0;

sum(Pt5t0(1:3))
