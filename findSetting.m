function [value,error] = findSetting(CFG,str)
value = -1;
error = 0;
settingfound = false;
for i=1:size(CFG,1)
    if strcmp(CFG(i,1),str)
        settingfound = true;
        value = str2double(CFG(i,2));
    end
end

if ~settingfound
    error = 1;
    disp(['Error:findSetting could not find setting ' str]);
    return
end

