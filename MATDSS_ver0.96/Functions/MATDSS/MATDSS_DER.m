function [MATDSS, DER] = MATDSS_DER(MATDSS)
% MATDSS_DER(app,MATDSS,DER)
% This function is responsible for creating Distributed Energy Resources (DERs)
% in the MATDSS application. Originally designed to interface with OpenDSS,
% it has been refactored to handle DER creation based on a table format
% configured through the application's settings window.
%
% Parameters:
%   - MATDSS: Structure containing simulation data and configurations.
%
% Outputs:
%   - MATDSS: Updated MATDSS structure after DER creation.
%   - DER: Structure containing details of created DERs.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
% Contact the developer at ilyas.farhat@outlook.com
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

DER = []; % Initialize DER structure
DER = struct; % Initialize as structure
MATDSS.Sim.nDER = 0; % Initialize number of DERs

% Loop through each row in the DER table to create DERs
for i = 1:size(MATDSS.TableData.DER, 1)
    [MATDSS, DER] = MATDSS_DERNew(MATDSS, DER, MATDSS.TableData.DER{i, :});
end

end % End of MATDSS_DER function


%% Old function code
%{

function [MATDSS, DER] = MATDSS_DER(MATDSS)
% [MATDSS, DER] = MATDSS_DER(app,MATDSS,DER) function is developed
% initially to creat DERs in OpenDSS. Since we moved to Table format of
% DERs information (through configuration window), this function has
% reformed into a simple loop that loops over the table entries.
%
%
% This function calls MATDSS_DERNew which creats the actual DERs loads in
% OpenDSS.
DER = [];
DER = struct;
MATDSS.Sim.nDER = 0;
for i = 1:size(MATDSS.TableData.DER,1)
    [MATDSS, DER] = MATDSS_DERNew(MATDSS,DER,MATDSS.TableData.DER{i,:});
end



end
%}