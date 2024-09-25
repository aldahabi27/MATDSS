v = MATDSS.Sim.DSSCircuit.AllBusVolts;
v = v(1:2:end) + 1i*v(2:2:end);
v = reshape(v,MATDSS.Sim.DSSCircuit.NumNodes,[]);

ANM = AllNodesNames;
v_base = MATDSS.Meas.VBases;
%comparing v with allbusvmag and allbusvmagpu
vpu = v./v_base;
v_mag = abs(v);

v_magpu = v_mag./(v_base*1e3);

[v_mag, AllBusVmag, v_mag-AllBusVmag];
disp('v diff')
max(v_mag-AllBusVmag)



[v_magpu, AllBusVmagPu, v_magpu-AllBusVmagPu];

disp('v pu diff')
max(v_magpu-AllBusVmagPu)



% Now we have checked that the voltages are correct. let's see how can we
% get the powers calculated correctly as in opendss. We will focus on the
% sourcebus now at first!

dssallbus = MATDSS.Sim.DSSCircuit.AllNodeNames;
% we see that sourcebus is the first three nodes in the list. hence the
% first three rows/cols are for the source bus in the ybus matrix!
ybus = MATDSS.Sim.Ybus;
s = v.*(conj(ybus*v));
v0 = v(1:3);
abs_v0 = abs(v0);
s0 = s(1:3);
p0 = real(s0);

%
mybus = MATDSS.Sim.DSSCircuit.Buses('sourcebus');
b_v = mybus.Voltages;
b_v = b_v(1:2:end) + 1i*b_v(2:2:end);
b_v = b_v';
abs_b_v = abs(b_v);
abs_v0;

mybuspce = mybus.AllPCEatBus
mybuspde = mybus.AllPDEatBus

mybus_cktelm = MATDSS.Sim.DSSCircuit.CktElements('Vsource.source')
mybus_cktelm.Name
reactor_s = mybus_cktelm.Powers;
reactor_s = reactor_s(1:2:end) + 1i*reactor_s(2:2:end);
reactor_s = reactor_s';
reactor_p = real(reactor_s);
reactor_q = imag(reactor_s);
reactor_p(1:3)'
[p0./1e3]'
disp('ratio')
[(p0./1e3)./reactor_p(1:3)]'


reactor_s = mybus_cktelm.TotalPowers;
reactor_s = reactor_s(1:2:end) + 1i*reactor_s(2:2:end);
reactor_s = reactor_s';
reactor_p = real(reactor_s)'
reactor_q = imag(reactor_s);
disp('ratio')
(sum(p0)./1e3)./reactor_p

%% using vpu instead of v (did not give equal results)





% Now we have checked that the voltages are correct. let's see how can we
% get the powers calculated correctly as in opendss. We will focus on the
% sourcebus now at first!

dssallbus = MATDSS.Sim.DSSCircuit.AllNodeNames;
% we see that sourcebus is the first three nodes in the list. hence the
% first three rows/cols are for the source bus in the ybus matrix!
ybus = MATDSS.Sim.Ybus;
v = vpu;
s = v.*(conj(ybus*v));
v0 = v(1:3);
abs_v0 = abs(v0);
s0 = s(1:3);
p0 = real(s0);

%
mybus = MATDSS.Sim.DSSCircuit.Buses('sourcebus');
b_v = mybus.Voltages;
b_v = b_v(1:2:end) + 1i*b_v(2:2:end);
b_v = b_v';
abs_b_v = abs(b_v);
abs_v0;

mybuspce = mybus.AllPCEatBus
mybuspde = mybus.AllPDEatBus

mybus_cktelm = MATDSS.Sim.DSSCircuit.CktElements('Vsource.source')
mybus_cktelm.Name
reactor_s = mybus_cktelm.Powers;
reactor_s = reactor_s(1:2:end) + 1i*reactor_s(2:2:end);
reactor_s = reactor_s';
reactor_p = real(reactor_s);
reactor_q = imag(reactor_s);
reactor_p(1:3)'
[p0./1e3]'
disp('ratio')
[(p0./1e3)./reactor_p(1:3)]'


reactor_s = mybus_cktelm.TotalPowers;
reactor_s = reactor_s(1:2:end) + 1i*reactor_s(2:2:end);
reactor_s = reactor_s';
reactor_p = real(reactor_s)'
reactor_q = imag(reactor_s);
disp('ratio')
(sum(p0)./1e3)./reactor_p



%% using default opendss ybus matrix!



Nnodes = MATDSS.Sim.DSSCircuit.NumNodes; % this includes all phases considered.
Ybus = MATDSS.Sim.DSSCircuit.SystemY; % Ymatrix contain only line and shunt admittances
Ybus = Ybus(1:2:end) + 1i*Ybus(2:2:end); % convert it complex Y matrix and then save it
Ybus = reshape(Ybus,Nnodes,[]); % Convert it to square matrix


% Now we have checked that the voltages are correct. let's see how can we
% get the powers calculated correctly as in opendss. We will focus on the
% sourcebus now at first!

dssallbus = MATDSS.Sim.DSSCircuit.AllNodeNames;
% we see that sourcebus is the first three nodes in the list. hence the
% first three rows/cols are for the source bus in the ybus matrix!
ybus = Ybus;
% v = vpu;
s = v.*(conj(ybus*v));
v0 = v(1:3);
abs_v0 = abs(v0);
s0 = s(1:3);
p0 = real(s0);

%
mybus = MATDSS.Sim.DSSCircuit.Buses('sourcebus');
b_v = mybus.Voltages;
b_v = b_v(1:2:end) + 1i*b_v(2:2:end);
b_v = b_v';
abs_b_v = abs(b_v);
abs_v0;

mybuspce = mybus.AllPCEatBus
mybuspde = mybus.AllPDEatBus

mybus_cktelm = MATDSS.Sim.DSSCircuit.CktElements('Vsource.source')
mybus_cktelm.Name
reactor_s = mybus_cktelm.Powers;
reactor_s = reactor_s(1:2:end) + 1i*reactor_s(2:2:end);
reactor_s = reactor_s';
reactor_p = real(reactor_s);
reactor_q = imag(reactor_s);
reactor_p(1:3)'
[p0./1e3]'
disp('ratio')
[(p0./1e3)./reactor_p(1:3)]'


reactor_s = mybus_cktelm.TotalPowers;
reactor_s = reactor_s(1:2:end) + 1i*reactor_s(2:2:end);
reactor_s = reactor_s';
reactor_p = real(reactor_s)'
reactor_q = imag(reactor_s);
disp('ratio')
(sum(p0)./1e3)./reactor_p
