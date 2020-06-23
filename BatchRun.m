%% Batch run
% Run main.m for n number of datasets
% BatchRun.m will search selected folers and all subfolders for sets of cnt, int, ext, and pho files
% output files from main.m will be created in the subfolders

% Usage: 
%   - Run BatchRun.m
%   - select folder(s) containing datasets

% Folder structure example:
% ->Dataset1
%   -> Dataset1.cnt
%   -> Dataset1.int
%   -> Dataset1.ext
%   -> Dataset1.pho
%   -> Cam0
%       -> Cam0.cnt
%       -> Cam0.int
%       -> Cam0.ext
%       -> Cam0.pho
% ->Dataset2
%   -> Cam0
%       -> Cam0.cnt
%       -> Cam0.int
%       -> Cam0.ext
%       -> Cam0.pho

% This folder structure will produce the following output in the data directories
% ->Dataset1
%   -> Dataset1_Dataset1.out
%   -> Cam0
%       -> Dataset1_Cam0.out
% ->Dataset2
%   -> Cam0
%       -> Dataset2_Cam0.out

% makes use of uigetfile_n_dir
% Reference: Tiago (2020). uigetfile_n_dir : select multiple files and directories
% (https://www.mathworks.com/matlabcentral/fileexchange/32555-uigetfile_n_dir-select-multiple-files-and-directories), MATLAB Central File Exchange.
% Retrieved June 23, 2020.


%selpath = uigetdir2(pwd,'Select Data Folders');
% for debugging:
clear all
clc
selpath = {'E:\OneDrive - University of Calgary\Summer2020\More Data\#Raw_Data_Only\GoPro dataset','E:\OneDrive - University of Calgary\Summer2020\More Data\#Raw_Data_Only\Ladybug dataset 1','E:\OneDrive - University of Calgary\Summer2020\More Data\#Raw_Data_Only\Ladybug dataset 2'};
exts = {'.pho','.ext','.cnt','.int'};

% search through folders and subfolders for all valid sets of data files
allfolders  = [];
for i = 1:length(selpath)
    folders = findfiles(selpath{i},exts);
    allfolders = [allfolders folders];
end

% loop through folder list
for i = 1:length(allfolders)
    folder = allfolders{i};
    % run main script with batch argument
    main(folder)
end














% find a specific set of files based on their extensions
function [folders] = findfiles(folder,file_exts)
% folder - starting folder (string)
% file_exts - set of file extensions to look for (cell of strings)

% filepaths - list of folders containing file sets

folders = [];

dirinfo = dir(folder); % get files/folders in folder
dirinfo=dirinfo(~ismember({dirinfo.name},{'.','..'})); % remove '.' and '..'

% check folder for file_exts
files_only = dirinfo(~[dirinfo.isdir]);
if ~isempty(files_only) % if the folder contains files
    %files_found = zeros(length(file_exts),1); % array to keep track of which files in file_exts have been found
    found = []; % keep track of which files were found
    for index = 1:length(files_only) % for each file in folder
        [~,~,ext] = fileparts(files_only(index).name);
        % check if ext is one of the extensions in file_exts
        for index_j = 1:length(file_exts)
            if strcmp(ext,file_exts{index_j}) % if extension matches
                for index_k = 1:length(found) % check that extension hasn't already been found
                    if strcmp(found{index_k},ext)
                        % if more than 1 file of a certain type was found, output warning to console
                        disp(['Warning: More than 1 ' ext ' file was found in ' folder]);
                        found = []; % reset found to supress 2nd output to console
                        break;
                    end
                end
                found = [found file_exts(index_j)];% update found array
                break;
            end
        end
    end
    
    % if no files were found, do nothing
    
    % if some files were found but not all, output warning to console
    if ~isempty(found) && length(found) < length(file_exts)
        % build error string
        
        str = 'Error: ';
        if length(found) > 1
            for index = 1:(length(found) - 1)
                str = [str found{index} ', '];
            end
            if length(found) == 2
                str(end-1) = '';
            end
            str = [str 'and ' found{end} ' were found in "' folder '" but not '];
        else
            str = [str found{1} ' was found in "' folder '" but not'];
        end
        notfound = setdiff(file_exts,found);
        if length(notfound) > 1
            for index = 1:(length(notfound) - 1)
                str = [str notfound{index} ', '];
            end
            if length(found) == 2
                str(end-1) = '';
            end
            str = [str 'and ' notfound{end}];
        else
            str = [str notfound{1}];
        end
        % display error
        disp([str '. This folder will be skipped']);   
      
    % if all files were found, add folder to folders[] list
    elseif length(found) == length(file_exts)
        folders = [folders {folder}];
    end
end

% run this function for each subfolder
subfolders = dirinfo([dirinfo.isdir]);
for index = 1:length(subfolders)
    folders_tmp = findfiles([subfolders(index).folder '\' subfolders(index).name],file_exts);
    folders = [folders folders_tmp];
end

%disp('end');
end