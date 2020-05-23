%% Read in text files for bundle adjuistment
% Wynand Tredoux, May 2020
% TODO: config files
function [terminate, filecontents] = ReadFiles(files)
%files = {'.pho','.ext','.cnt','.int'};
filecontents = cell(length(files),1);
terminate  = 0;
% files to read in:
    % .PHO: Image point observations Coordinates in pixels
        % format: Point image x y
    % .EXT: Approximate exterior orientation parameters Coordinates in mm Angles in decimal degrees
        % format: Image camera Xc Yc Zc w f k
    % .CNT: Approximate values for object point coordinates Coordinates in mm
        % format: Point X Y Z
    % .INT: Interior orientation parameters
        % format: Camera y_axis_dir x_min y_min x_max y_max
        %         xp yp c
    % .CFG: Defines verious optional parameters for the program
    % .Defaultcfg: Default cfg options that can be overwritten by options in .CFG
    
% read in .pho .ext .cnt .int
for i=1:length(files)
   file = files{i};
   % get filepath for current file
   dirfiles = dir(strcat('*',file));
   if length(dirfiles)>1 % if more than 1 file is found
       error = errordlg(['Error: multiple files with ' file ' extension, please select 1 to use'],['Too many ' file ' files']);
       uiwait(error)
       filepath = uigetfile(strcat('*',file),['Select ' file]);
   elseif length(dirfiles)<1 % if no file is found
       answer = questdlg(['No ' file ' in current directory. Would you like to browse for it in another location?'],['No ' file], 'Yes','Cancel','Yes');
       if strcmp(answer,'Cancel')
           %disp("TERMINATE");
           terminate = 1;
           break;
       else
           [file, path] = uigetfile(strcat('*',file),['Select ' file]);
           if file == 0
               terminate = 1;
               break;
           end
           filepath = fullfile(path,file); 
       end
   else
       filepath = dirfiles.name;
   end
   % read in file
    input = readmatrix(filepath,'FileType','text','NumHeaderLines',0,'Delimiter',{' ','\t'},'ConsecutiveDelimitersRule','join','LeadingDelimitersRule','ignore','OutputType','string','CommentStyle','#*');
    filecontents(i,1) = {input};
end
end