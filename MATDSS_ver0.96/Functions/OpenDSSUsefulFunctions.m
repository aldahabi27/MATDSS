%% Check Line properties
Linename = 'line2';
success = MATDSS.Sim.DSSCircuit.SetActiveElement(['Line.' Linename]); % set active element to line k
if ~success
    disp('Error in setting active element');
end
MyLine = MATDSS.Sim.DSSCircuit.ActiveElement; % handle of the current line i
MyLineBuses = MyLine.BusNames;
MyLineS = MyLine.Powers;
MyLineS = MyLineS(1:2:end)+1i.*MyLineS(2:2:end);
MyLineS = reshape(MyLineS,3,[])';
MyLineS1 = MyLineS(1,1:3);
MyLineP1 = real(MyLineS1);
MyLineQ1 = imag(MyLineS1);



[MyLineP1, sum(MyLineP1);MyLineQ1, sum(MyLineQ1)]