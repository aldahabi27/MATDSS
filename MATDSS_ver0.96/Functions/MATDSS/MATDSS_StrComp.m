% MATDSS_StrComp(A,B)
% This function compares two strings or string arrays based on specific cases.
% The function can handle different scenarios including:
% 
% Case 1: Both A and B are individual strings.
%   - Output: eqflag = 1 if A = B, otherwise eqflag = 0.
%
% Case 2: A is a cell array of strings and B is a single string.
%   - Output: eqflag = index of the matching string in A; -1 if no match is found.
%   - If there are multiple matches, only the first index is returned.
%
% Case 3: Both A and B are cell arrays of strings.
%   - Output: eqflag = an array of indices indicating the positions of elements in B that match elements in A.
%   - Entries in eqflag will be -1 for any element in B that does not exist in A.
%
% The function also includes checks to handle empty inputs gracefully.
%
% Last Update for this function was on MATDSS App Ver 0.93 (29 Sept. 2023)
%
% MATDSS Application
% Copyright (c) 2024, Ilyas Farhat
%
% This file is part of MATDSS Application
% Contact the developer at ilyas.farhat@outlook.com

function eqflag = MATDSS_StrComp(A,B)
    % Initialize default string case assumption and output flag
    StrCase = 1; % By default, A and B are considered to be individual strings
    eqflag = -1; % Default value indicating no match found

    % Check if A is a cell array with more than one element
    if iscell(A) && max(size(A)) > 1
        % Determine if B is also a cell array
        if iscell(B) && max(size(B)) > 1
            StrCase = 3; % Both A and B are cell arrays
        else
            StrCase = 2; % A is a cell array of strings and B is a single string
        end
    end

    % Check for empty inputs
    if min(min(size(A),min(size(B)))) < 1
        eqflag = -1; % If either A or B is empty, set eqflag to -1
        return; % Exit the function
    end

    % Main switch case to handle different string comparison scenarios
    switch StrCase
        case 1 % A and B are individual strings
            if iscell(A) % If A is a cell, convert it to string
                A = A{1}; % Extract the string from the cell
            end
            if iscell(B) % If B is a cell, convert it to string
                B = B{1}; % Extract the string from the cell
            end
            eqflag = strcmpi(A,B); % Compare strings case-insensitively
            
        case 2 % A is a cell array and B is an individual string
            eqflag = find(strcmpi(A, B),1); % Find the index of the matching string in A
            if isempty(eqflag)
                eqflag = -1; % If no match is found, set eqflag to -1
            end
            
        case 3 % A and B are cell arrays
            [~, eqflag] = ismember(B,A); % Get indices of B elements that are found in A
    end

    % clearvars -except eqflag % Clear unnecessary variables while keeping eqflag
end

%% Test code
%{
A = {'hi', 'pi', 'two', 'tahseen', 'ilyas'}; % Test array A
B = {'ilyas','hi','ilyas','saeed'}; % Test array B
MATDSS_StrComp(A,B) % Test function call
B = {'hi','ilyas','saeed'}; % Another test array B
MATDSS_StrComp(A,B) % Test function call
B = {'ilyas'}; % Single test case for B
MATDSS_StrComp(A,B) % Test function call
B = {'saeed'}; % Non-matching case for B
MATDSS_StrComp(A,B) % Test function call
A = 'hi'; B = 'hi'; MATDSS_StrComp(A,B) % Individual string comparison
A = {'hi'}; B = 'hi'; MATDSS_StrComp(A,B) % Cell to string comparison
A = 'hi'; B = {'hi'}; MATDSS_StrComp(A,B) % String to cell comparison
%}


%% Old function code

%{

function eqflag = MATDSS_StrComp(A,B)
% eqflag = MATDSS_StrComp(A,B) compares two strings as follows:
% --------------------------------------------------------------------------
% Case 1: A and B are both strings (each is one string)
% 
% Output: eqflag = 1 if A = B and 0 else
%
%
%
% Case 2: A is a cell array of strings and B is an individual string
%
% Output: eqflag = index, where index is the string index in A that matches
% B. If B is not equal to any A elements, eqflag = -1. if more than one
% element in A is equal to B, then the function will return the first one
% only!
%
%
%
% Case 3: A and B are both cell array of strings
%
% Output: eqflag = index_array, where index array will be a vector of
% integers pointing to the index of elements in A that matches the elements
% in B; in the order of there appearance in B. index_array will have -1
% entries for any element in B that is not in A.

StrCase = 1; %by default, A and B are considered to be individual strings
eqflag = -1; % default value

if iscell(A) && max(size(A)) > 1
    if iscell(B) && max(size(B)) > 1
        StrCase = 3; % A and B are both cell arrays
    else
        StrCase = 2; % A is cell array of strings and B is single string
    end
end

if min(min(size(A),min(size(B)))) < 1
    eqflag = -1;
    return;
end


switch StrCase
    case 1 % A and B are individual string
        if iscell(A) % if A is cell, convert it to string
            A = A{1};
        end
        if iscell(B) % if B is cell, convert it to string
            B = B{1};
        end
        eqflag = strcmpi(A,B);
        
    
    
    case 2 % A is cell array and B is an individual string
        eqflag = find(strcmpi(A, B),1);
        if isempty(eqflag)
            eqflag = -1;
        end
        % if iscell(B) % if B is cell, convert it to string
        %     B = B{1};
        % end
        % for i = 1:max(size(A))
        %     if strcmp(A{i}, B)
        %         eqflag = i;
        %         break;
        %     end
        % end

    case 3 % A and B are cell arrays
       
        [~, eqflag] = ismember(B,A);
        % eqflag = eqflag.*ones(size(B));
        % for j = 1:max(size(B))
        %     for i = 1:max(size(A))
        %         if strcmp(A{i},B{j})
        %             eqflag(j) = i;
        %             break;
        %         end
        %     end
        % end
end

% clearvars -except eqflag
end



%% Test code
%{
A = {'hi', 'pi', 'two', 'tahseen', 'ilyas'}
B = {'ilyas','hi','ilyas','saeed'};MATDSS_StrComp(A,B)
B = {'hi','ilyas','saeed'};MATDSS_StrComp(A,B)
B = {'ilyas'};MATDSS_StrComp(A,B)
B = {'saeed'};MATDSS_StrComp(A,B)
A = 'hi'; B = 'hi'; MATDSS_StrComp(A,B)
A = {'hi'}; B = 'hi'; MATDSS_StrComp(A,B)
A = 'hi'; B = {'hi'}; MATDSS_StrComp(A,B)
%}


%}