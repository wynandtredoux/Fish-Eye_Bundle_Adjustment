%% find setting in CFG file
% Inputs:
    % CFG - n by 2 string matrix containing setting names/values from .cfg file
    % str - string containing name of setting who's value to find
    % error_tally - integer used for keeping track of how many errors occurred while reading settings
    % Check_truefalse - Optional boolean. If true function will check if setting value is 1/0, and return an error if not
function [value,error_tally] = findSetting(CFG,str,error_tally,Check_truefalse)
if nargin < 4
    Check_truefalse = 0;
end
value = -1;
settingfound = false;

% find setting
for i=1:size(CFG,1)
    if strcmp(CFG(i,1),str)
        settingfound = true;
        value = str2double(CFG(i,2));
    end
end

% if setting cannot be found
if ~settingfound
    error_tally = error_tally + 1;
    disp(['Error:findSetting() could not find setting ' str]);
    return
end

% check for NaN
if isnan(value)
    error_tally = error_tally + 1;
    disp(['Error:findSetting() ' str ' cannot be NaN. This may be caused by an invalid setting value (setting values must contain numbers only)']);
    return
end

% check if setting value is 1/0 (if applicable)
if Check_truefalse
    if value ~= 1 && value ~= 0
        error_tally = error_tally + 1;
        disp(['Error:findSetting() ' str ' must be 1 or 0']);
        return
    end
end

end

