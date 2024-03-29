%% Bundle adjustment script using equidistant fish-eye model
% Wynand Tredoux, May 2020

% Usage:  main('folder',plot) - runs the adjustment in 'folder'. 'plot' 
%             defines if the program should draw plots and can be true (1) or false (0)
%         main() - uses project folder with plots enabled
%         main('filepath/myfolder') - runs the adjustment on the given folder
%         main('', 0) - runs the adjustment on project folder without plots

function main_error = main(varargin)
%% define parser object
p = inputParser;
addOptional(p, 'folder', '', @ischar);
addOptional(p, 'plot', true, @(x) islogical(x));
parse(p,varargin{:})

%% Setup
%clear all
close all
format longg
%clc
celldisp(varargin)
main_error = 0;

enable_plots = p.Results.plot
batch = 0; 
folder = p.Results.folder

% set batch mode if folder is not ''
if ~strcmp(folder,'')
    batch = 1;
end
        
projectDir = pwd;
%addpath([projectDir '\functions']); % add functions folder to path

date = char(datetime); %date
% get version of code using the "git describe" command
[gite, version] = system('git describe --dirty'); % "dirty" indicates that the local version has been modified and doesn't fully match the version on github
if gite>0
    disp('warning: could not get git version with "git describe"')
    disp('Version will be set to "unknown"');
    version = 'Unknown\n';
end
mfiles = '';
% if version has been modified
if contains(version,'dirty')
    % get list of modified files
    [~, mfiles] = system('git ls-files -m');
end
%% read in files
%PHO image measurements [pointID image xmm ymm]
%EXT EOPs [imageID CamreaID Xc Yc Zc omega phi kappa]
%CNT object coordinates [TargetID X Y Z]
%INT IOPs [CameraID yaxis_dir xmin ymin xmax ymax]
%          [xp yp c k1 k2 k3 ... p1 p2] <- if k1,k2,k3... and p1,p2 are omitted, they will be assumed to be 0
%TIE = []; % list of tie point target IDs (read in later if needed)
%CFG config file

% if main is being run in batch mode
if batch
    % need to add the project directory to MATLAB's path so functions still work
    addpath(projectDir)
    
    % check to see if there is a cfg file in the directory if given
    files = dir(folder);
    files = files(~[files.isdir]);
    cfgfound = 0;
    for i = 1:length(files)
        [~, ~, ext] = fileparts(files(i).name);
        if strcmp(ext,'.cfg')
            cfgfound = 1;
            break;
        end
    end
    % if there no cfg file, read in the one in project directory
    if ~cfgfound
        [filereaderror, files] = ReadFiles({'.cfg'});
        CFG = files{1};
    else
        % change dir and read the CFG in the data folder
        cd(folder);
        [filereaderror, files] = ReadFiles({'.cfg'});
        CFG = files{1};
    end
    
    cd(folder)
    % read in the rest of the files
    [tmp, files] = ReadFiles({'.pho','.ext','.cnt','.int'});
    filereaderror = filereaderror + tmp;
    
else % if main is being run normally  
    [filereaderror, files] = ReadFiles({'.pho','.ext','.cnt','.int','.cfg'});
    CFG = files{5};
end

if filereaderror >= 1
    disp ('Error reading files')
    main_error = 1;
    return
else
    disp('Files read successfully!')
end

PHO = files{1}; 
EXT = files{2}; 
CNT = files{3}; 
INT = files{4};
TIE = [];
CZE = [];

%% Get settings from cfg file
% Settings that are allowed to be missed (have some default values defines here)
% get output filename from cfg (is allowed to error without exiting the program)
cfg_errors = 0;
[data.settings.Output_Filename,cfg_errors] = findSetting(CFG,'Output_Filename',cfg_errors);
if cfg_errors>0 % default filename if none is provided
    [~, Output_Filename, ~] = fileparts(pwd); % set to folder name
    data.settings.Output_Filename = strcat(Output_Filename,'.out');
end
clear Output_Filename
% get measurement standard deviation (is allowed to error without exiting the program)
cfg_errors = 0;
[data.settings.Meas_std,cfg_errors] = findSetting(CFG,'Meas_std',cfg_errors);
if cfg_errors>0 % if no std is provided, set to 1
    data.settings.Meas_std = 1;
    data.settings.no_std_y = 1;
else
    [data.settings.Meas_std_y,data.settings.no_std_y] = findSetting(CFG,'Meas_std_y',cfg_errors);
end

% get model type (is allowed to error without exiting the program)
cfg_errors = 0;
[data.settings.type,cfg_errors] = findSetting(CFG,'Type',cfg_errors);
if cfg_errors>0 % if no type is provided, set to fisheye
    data.settings.type = 'fisheye';
end
disp(['Type set to ' data.settings.type]);

% Calculate check points
cfg_errors = 0;
[data.settings.Check_Points,cfg_errors] = findSetting(CFG,'Check_Points',cfg_errors,1);
if cfg_errors>0 % default to no check points
    data.settings.Check_Points = 0;
end
 
cfg_errors = 0;
% settings that *must* be provided (program will exit if they are missing):
% General Settings
[data.settings.Iteration_Cap,cfg_errors] = findSetting(CFG,'Iteration_Cap',cfg_errors);
[data.settings.threshold,cfg_errors] = findSetting(CFG,'Threshold_Value',cfg_errors);
% inner constraints
[data.settings.Inner_Constraints,cfg_errors] = findSetting(CFG,'Inner_Constraints',cfg_errors,1);
% Estimate EOPs
[data.settings.Estimate_Xc,cfg_errors] = findSetting(CFG,'Estimate_Xc',cfg_errors,1);
[data.settings.Estimate_Yc,cfg_errors] = findSetting(CFG,'Estimate_Yc',cfg_errors,1);
[data.settings.Estimate_Zc,cfg_errors] = findSetting(CFG,'Estimate_Zc',cfg_errors,1);
[data.settings.Estimate_w,cfg_errors] = findSetting(CFG,'Estimate_Omega',cfg_errors,1);
[data.settings.Estimate_p,cfg_errors] = findSetting(CFG,'Estimate_Phi',cfg_errors,1);
[data.settings.Estimate_k,cfg_errors] = findSetting(CFG,'Estimate_Kappa',cfg_errors,1);
% Estimate IOPs
[data.settings.Estimate_c,cfg_errors] = findSetting(CFG,'Estimate_c',cfg_errors,1);
[data.settings.Estimate_xp,cfg_errors] = findSetting(CFG,'Estimate_xp',cfg_errors,1);
[data.settings.Estimate_yp,cfg_errors] = findSetting(CFG,'Estimate_yp',cfg_errors,1);
% Estimate Distortions
[data.settings.Estimate_radial,cfg_errors] = findSetting(CFG,'Estimate_Radial_Distortions',cfg_errors,1);
[data.settings.Num_Radial_Distortions,cfg_errors] = findSetting(CFG,'Num_Radial_Distortions',cfg_errors);
[data.settings.Estimate_decent,cfg_errors] = findSetting(CFG,'Estimate_Decentering_Distortions',cfg_errors,1);
% Estimate Ground Coordinates of tie points
[data.settings.Estimate_tie,cfg_errors] = findSetting(CFG,'Estimate_tie',cfg_errors,1);
[data.settings.Estimate_AllGCP,cfg_errors] = findSetting(CFG,'Estimate_AllGCP',cfg_errors,1);

if cfg_errors>0
    disp ('Error getting settings')
    main_error = 1;
    return
end

%% Read in tie point list if needed
if data.settings.Estimate_tie == 1 && data.settings.Estimate_AllGCP == 0
    [filereaderror, files] = ReadFiles({'.tie'});
    if filereaderror == 1
        disp ('Error reading files')
        main_error = 1;
        return
    end
    TIE = files{1};
end

%% Estimate_AllGCP
% if data.settings.Estimate_AllGCP == 1
%     TIE = CNT(:,1); % add all targets found in PHO as TIE points
%     data.settings.Estimate_tie = 1; % change Estimate_tie to 1 so all GCPs are estimated
% end

%% Convert files from String to Cell
% PHO
tmp = cell(size(PHO,1),size(PHO,2));
for i = 1:size(PHO,1)
    tmp(i,1) = {char(PHO(i,1))};
    tmp(i,2) = {char(PHO(i,2))};
    tmp(i,3) = {str2double(PHO(i,3))}; 
    tmp(i,4) = {str2double(PHO(i,4))};
end
PHO = tmp;
% EXT
tmp = cell(size(EXT,1),size(EXT,2));
for i = 1:size(EXT,1)
    tmp(i,1) = {char(EXT(i,1))};
    tmp(i,2) = {char(EXT(i,2))};
    tmp(i,3) = {str2double(EXT(i,3))}; 
    tmp(i,4) = {str2double(EXT(i,4))};
    tmp(i,5) = {str2double(EXT(i,5))};
    % convert w p k from degrees to radians
    tmp(i,6) = {str2double(EXT(i,6))*pi()/180}; 
    tmp(i,7) = {str2double(EXT(i,7))*pi()/180};
    tmp(i,8) = {str2double(EXT(i,8))*pi()/180};
end
EXT = tmp;
% CNT
tmp = cell(size(CNT,1),size(CNT,2));
for i = 1:size(CNT,1)
    tmp(i,1) = {char(CNT(i,1))};
    tmp(i,2) = {str2double(CNT(i,2))};
    tmp(i,3) = {str2double(CNT(i,3))}; 
    tmp(i,4) = {str2double(CNT(i,4))};
end
CNT = tmp;
% INT
tmp = cell(size(INT,1),size(INT,2));
for i = 1:2:size(INT,1)
    % row 1 
    tmp(i,1) = {char(INT(i,1))};
    tmp(i,2) = {str2double(INT(i,2))};
    tmp(i,3) = {str2double(INT(i,3))}; 
    tmp(i,4) = {str2double(INT(i,4))};
    tmp(i,5) = {str2double(INT(i,5))};
    tmp(i,6) = {str2double(INT(i,6))};
    % row 2
    for j = 1:3
        tmp(i+1,j) = {str2double(INT(i+1,j))};
    end
    % if no distortion parameters are provided in INT, set them to 0
    for j = 4:(5+data.settings.Num_Radial_Distortions)
        if size(INT,2)>=j % check that elements exist in INT 
            if ~ismissing(INT(i+1,j)) % check that they are not missing
                tmp(i+1,j) = {str2double(INT(i+1,j))};
            else
                tmp(i+1,j) = {0};
            end
        else
            tmp(i+1,j) = {0};
        end
    end
    
end
INT = tmp;
clear tmp

%% Estimate_AllGCP
if data.settings.Estimate_AllGCP == 1
    TIE = unique(PHO(:,1)); % add all targets found in PHO as TIE points
    data.settings.Estimate_tie = 1; % change Estimate_tie to 1 so all GCPs are estimated
end

%% Read in check points if needed
if data.settings.Check_Points
    [filereaderror, files] = ReadFiles({'.cze'});
    if filereaderror == 1
        disp ('Error reading files')
        main_error = 1;
        return
    end
    CZE = files{1};
end

%% Build points structure
%points = struct('x',zeros(size(PHO,1),1),'y',zeros(size(PHO,1),1),'targetID');
% Loop through PHO measurements and assign parameters to each measurement
for i = 1:size(PHO,1) % for each measurement
    %% Get photo measurements
    data.points(i).x = PHO{i,3};% x coordinate on the image
    data.points(i).y = PHO{i,4};% y coordinate on the image
    data.points(i).targetID = PHO{i,1}; % target name/id
    data.points(i).imageID = PHO{i,2}; % image ID
    %% find matching EOPs from EXT
    ext_index = -1;
    for j = 1:length(EXT)
        if strcmp(EXT{j,1},data.points(i).imageID)
            ext_index = j;
            break;
        end        
    end
    if ext_index < 0
        err_dlg = errordlg(['Could not find image ' data.points(i).imageID ' from .pho in .ext. Check that the image ID exists in both files'],['Error on image ' data.points(i).imageID]);
        waitfor(err_dlg);  % code execution is stopped till err_dlg is dismissed.
        return
    end
    data.points(i).ext_index = ext_index;
    % get cameraID in EXT
    data.points(i).cameraID = EXT{ext_index,2};
    % get EOPs
    data.points(i).Xc = EXT{ext_index,3};
    data.points(i).Yc = EXT{ext_index,4};
    data.points(i).Zc = EXT{ext_index,5};
    data.points(i).w = EXT{ext_index,6};
    data.points(i).p = EXT{ext_index,7};
    data.points(i).k = EXT{ext_index,8};
    %% find matching IOPs from INT
    int_index = -1;
    for j = 1:2:length(INT)
        if strcmp(INT(j,1),data.points(i).cameraID)
            int_index = j;
            break;
        end        
    end
    if int_index < 0
        err_dlg = errordlg(['Could not find camera ' data.points(i).cameraID ' from .ext in .int. Check that the camera ID exists in both files'],['Error on camera ' data.points(i).cameraID]);
        waitfor(err_dlg);  % code execution is stopped till err_dlg is dismissed.
        return
    end
    data.points(i).int_index = int_index;
    data.points(i).cam_num = (int_index + 1)/2;
    % xp, yp, c
    data.points(i).xp = INT{int_index + 1,1};
    data.points(i).yp = INT{int_index + 1,2};
    data.points(i).c = INT{int_index + 1,3};
    % distortions
    data.points(i).K = [INT{int_index + 1,4:4+data.settings.Num_Radial_Distortions-1}]';
    data.points(i).P = [INT{int_index + 1,4+data.settings.Num_Radial_Distortions:5+data.settings.Num_Radial_Distortions}]';
    % get y_axis_dir
    y_dir = INT{int_index,2};
    % Check that y_dir is +-1
    if y_dir ~= 1 && y_dir ~= -1
        disp('y_dir should be +-1 only')
        return
    end
    % get x/y min/max
    data.points(i).xmin = INT{int_index,3};
    data.points(i).ymin = INT{int_index,4};
    data.points(i).xmax = INT{int_index,5};
    data.points(i).ymax = INT{int_index,6};
    data.points(i).y_dir = y_dir;
    %% find ground coordinates from CNT
    cnt_index = -1;
    for j = 1:length(CNT)
        if strcmp(CNT(j,1),data.points(i).targetID)
            cnt_index = j;
            break;
        end        
    end
    if cnt_index < 0
        err_dlg = errordlg(['Could not find target ' data.points(i).targetID ' from .pho in .cnt. Check that the target ID exists in both files'],['Error on target ' data.points(i).targetID]);
        waitfor(err_dlg);  % code execution is stopped till err_dlg is dismissed.
        return
    end
    % get object coordinates
    data.points(i).cnt_index = cnt_index;
    data.points(i).X = CNT{cnt_index,2};
    data.points(i).Y = CNT{cnt_index,3};
    data.points(i).Z = CNT{cnt_index,4};
    %% check if this is a tie point
    data.points(i).tieIndex = -1;
    % find current targetID in TIE
    for j = 1:size(TIE,1)
        if strcmp(TIE(j),data.points(i).targetID)
            data.points(i).tieIndex = j;
            break;
        end
    end
    
    if data.points(i).tieIndex ~= -1 % if targetID was found in TIE, then it is a TIE point and must be estimated
        data.points(i).isTie = 1;
    else
        data.points(i).isTie = 0;
    end
    
end
data.numImg = length(unique({data.points(:).imageID}));
data.numCam = length(unique({data.points(:).cameraID}));
data.n = size(data.points,2)*2;
data.numGCP = length(unique([data.points(:).cnt_index]));
data.numtie = length(TIE);
clear PHO y_dir ext_index int_index cnt_index tieIndex
%% Main Loop
tic %start time
% build xhat (initial values for unknowns)
[xhaterror, xhat, xhatnames] = Buildxhat(data, EXT, INT, TIE, CNT);
if xhaterror == 1
    disp('Error building xhat');
    main_error = 1;
    return
end


% build P weight matrix from Meas_std and Meas_std_y if provided
if data.settings.no_std_y
    Cl = diag(repmat(data.settings.Meas_std^2,size(data.points,2)*2,1));
    data.settings = rmfield(data.settings,'Meas_std_y');
else
    Cl = diag(repmat([data.settings.Meas_std^2; data.settings.Meas_std_y^2;],size(data.points,2),1));
end
% set a-priori variance factor
priori = 1;
P = priori*Cl^-1;

deltasum = 100;
count = 0;
xhat_arr = xhat';
deltasumarr = [];
% main loop
while deltasum > data.settings.threshold
    count = count + 1;
    disp(['Iteration ' num2str(count) ':']);
    %% build A and w matrices (and G if applicable)
    [Awerror, A, w, G, dist_scaling] = BuildAwG(data, xhat);
    if Awerror == 1
        disp('Error building A and w');
        main_error = 1;
        return
    end
    
    %% Calculate solution
    u = A'*P*w;
    N = A'*P*A;
    %Cx = [];
    
    if data.settings.Inner_Constraints
        NG = [N G;
            G' zeros(size(G,2));];
        uG = [u; zeros(size(G,2),1);];
        Cx = NG^-1;
        
        % separate k and delta
        delta_k = -Cx*uG;
        delta = delta_k(1:size(u,1),:);
        %k = delta_k(size(u,1)+1:end,:);
        
        % separate Cx into only the unknowns part
        Cx = Cx(1:size(A,2),1:size(A,2));
    else
        Cx = N^-1;
        delta = -Cx*u;
    end
    
    % Correlation matrix (this needs to be done before distortion re-scaling)
    Correlation = zeros(size(Cx,1),size(Cx,2));
    for i = 1:size(Cx,1)
        for j = 1:size(Cx,2)
            sigmaX = sqrt(Cx(i,i));
            sigmaY = sqrt(Cx(j,j));
            sigmaXY = Cx(i,j);

            Correlation(i,j) = sigmaXY / (sigmaX*sigmaY);
        end
    end
    
    % scale distortion parameters and variances in delta/Cx with values from dist_scaling
    % for each camera
    for i = 1:size(dist_scaling,1)
        % radial distortion
        if data.settings.Estimate_radial
            % get radial distortion index
            radial_index = dist_scaling(i,1);
            % scale parameters
            for j = 1:data.settings.Num_Radial_Distortions
                delta(radial_index+j-1) = delta(radial_index+j-1)/dist_scaling(i,j+2);
                Cx(radial_index+j-1,radial_index+j-1) = Cx(radial_index+j-1,radial_index+j-1)/(dist_scaling(i,j+2).^2);
            end
        end
        % decentering distortion
        if data.settings.Estimate_decent
            % get decentering distortion index
            decent_index = dist_scaling(i,2);
            % scale P1
            delta(decent_index) = delta(decent_index)/dist_scaling(i,3);
            Cx(decent_index,decent_index) = Cx(decent_index,decent_index)/(dist_scaling(i,3).^2);
            % scale P2
            delta(decent_index+1) = delta(decent_index+1)/dist_scaling(i,3);
            Cx(decent_index+1,decent_index+1) = Cx(decent_index+1,decent_index+1)/(dist_scaling(i,3).^2);
        end
    end
    
    xhat = xhat + delta;
    xhat_arr = [xhat_arr; xhat';];
        
    deltasum = sumabs(delta)
    deltasumarr = [deltasumarr deltasum];
    % constrain loop iterations
    if count >= data.settings.Iteration_Cap
        disp('Iteration Cap reached. This can be changed in the .cfg file')
        break;
    end
end
% end time
time = num2str(toc);
disp(['Elapsed time is ' char(time) ' seconds.'])

[~,name,~] = fileparts(data.settings.Output_Filename);

%plot normal of delta over iterations
if enable_plots
    fig = figure;
    hold on
    plot(deltasumarr)
    title('normal of \delta over time')
    xlabel('Iteration')
    ylabel('normal value')
    
    saveas(fig,strcat('delta_',name,'.png'));
    close(fig);
end

%plot Xc, Yc, Zc over iterations
if enable_plots
    fig = figure;
    hold on
    legend_str = [];
    % Xc
    if data.settings.Estimate_Xc
        plot(xhat_arr(:,1))
        legend_str = [legend_str; 'Xc';];
    end
    % Yc
    if data.settings.Estimate_Yc
        plot(xhat_arr(:,data.settings.Estimate_Xc + 1))
        legend_str = [legend_str; 'Yc';];
    end
    % Zc
    if data.settings.Estimate_Zc
        plot(xhat_arr(:,data.settings.Estimate_Xc + data.settings.Estimate_Yc +1))
        legend_str = [legend_str; 'Zc';];
    end
    legend(legend_str)
    
    saveas(fig,strcat('XcYcZc_',name,'.png'));
    close(fig);
end

%plot w, p, k over iterations
if enable_plots
    fig = figure;
    hold on
    legend_str = [];
    offset = data.settings.Estimate_Xc + data.settings.Estimate_Yc + data.settings.Estimate_Zc;
    % Xc
    if data.settings.Estimate_w
        plot(xhat_arr(:,offset + 1))
        legend_str = [legend_str; 'w';];
    end
    % Yc
    if data.settings.Estimate_p
        plot(xhat_arr(:,offset + data.settings.Estimate_w + 1))
        legend_str = [legend_str; 'p';];
    end
    % Zc
    if data.settings.Estimate_k
        plot(xhat_arr(:,offset + data.settings.Estimate_w + data.settings.Estimate_k + 1))
        legend_str = [legend_str; 'k';];
    end
    legend(legend_str)
    
    saveas(fig,strcat('wpk_',name,'.png'));
    close(fig);
end

%% Residuals
% calculate x,y residuals for all image measurements
v = A*delta + w;
% build RSD and corrected measurements with the format point_id image _id x_obs y_obs radial_dist v_x v_y v_r v_t
RSD = BuildRSD(v, data, xhat);

% create plot of radial component of the residuals - RSD(:,8) - as a function of radial distance - RSD(:,5)
if enable_plots
    fig = figure;
    hold on
    scatter([RSD{:,5}],[RSD{:,8}]);
    title('v_r vs r')
    xlabel('radial distance r')
    ylabel('radial component of the residuals v_r')
    % save figure
    saveas(fig,strcat('RSDvR_',name,'.png'));
    close(fig);
end

%% corrected image coordinates
for i = 1:length(data.points)
    data.points(i).x_corr = data.points(i).x + RSD{i,6}; % correct x coords
    data.points(i).y_corr = data.points(i).y + RSD{i,7}; % correct y coords
end

%% RMS
% RMSx: x residuals are all the odd indices of v
RMSx = rms(v(1:2:end));
% RMSy: y residuals are all the even indices of v
RMSy = rms(v(2:2:end));

RMS = sqrt(RMSx^2 + RMSy^2);

%%  variance factor
sigma02 = v'*P*v/(size(A,1)-size(A,2))
Cx = sigma02.*Cx;

%% Check point differences - if used
if data.settings.Check_Points
    % this is just the difference between the given check point CZE coordinates and
    % the corresponding calculated object points
    cp_diff = cell(size(CZE,1), size(CZE,2));
    for i = 1:size(CZE,1)
        % get name of check point
        cp_diff(i,1) = {CZE{i,1}};
        % find index of X coord in xhat using xhatnames
        index = find(strcmp([xhatnames{:}],['X_' CZE{i,1}]));
        if isempty(index) % if no points were found
            disp(['Warning: Check point not found in xhat -> ' CZE{i,1}])
        elseif length(index) > 1 % if more than 1 point was found
            disp(['Warning: More than 1 match found for -> ' ['X_' CZE{i,1}]])
        else % exactly 1 point was found
            % calculate check point difference
            XYZ = [xhat(index) xhat(index+1) xhat(index+2)]; % estimated XYZ
            XYZ_cp = str2double([CZE(i,2:4)]); % measured XYZ
            cp_diff(i,2:4) = num2cell([XYZ - XYZ_cp]); % (estimated - measured)
        end
    end
    % calculate mean and RMSE of check point differences
    cp_mean = [mean([cp_diff{:,2}]) mean([cp_diff{:,3}]) mean([cp_diff{:,4}])];
    cp_rms = [rms([cp_diff{:,2}]) rms([cp_diff{:,3}]) rms([cp_diff{:,4}])];
end


%% Create output file
padding = 4;
disp('Writing output file...');
line = '*************************************************************************************************************';

% create output file
fileID = fopen(data.settings.Output_Filename,'w');

% heading
fprintf(fileID,['Version: ' version]);
fprintf(fileID,'Fish-eye model Bundle Adjustment\nWynand Tredoux -- University of Calgary -- 2020\n\n');
if ~isempty(mfiles) % if modified files exist
    fprintf(fileID,['Modified files:\n' mfiles]); % print list of modified files
end
fprintf(fileID,line);
fprintf(fileID,['\n\nExecution date:\t' date '\nTime Taken:\t\t' char(time) ' seconds\nIterations:\t\t' num2str(count) '\nModel Used:\t\t' data.settings.type]);

% settings used
fprintf(fileID,'\n\nSettings used:\n');
printCell(fileID, [fieldnames(data.settings) struct2cell(data.settings)], '\t\t', padding);
fprintf(fileID, ['\n' line '\n']);

% Summary of unknowns/observations
fprintf(fileID, '\nObservations/Unknowns Summary\n\n');

tmp = [{'Number of Photos'}	{num2str(data.numImg)}
{'Total EOP unknowns'}	{num2str((data.settings.Estimate_Xc + data.settings.Estimate_Yc + data.settings.Estimate_Zc + data.settings.Estimate_w + data.settings.Estimate_p + data.settings.Estimate_k)*data.numImg)}
{'Number of Cameras'}	{num2str(data.numCam)}
{'Total IOP unknowns'}	{num2str((data.settings.Estimate_c + data.settings.Estimate_xp + data.settings.Estimate_yp)*data.numCam)}
{'Total distortion unknowns'}	{num2str((data.settings.Estimate_radial*data.settings.Num_Radial_Distortions + data.settings.Estimate_decent*2)*data.numCam)}
{'Number of tie/control points'}	{num2str(data.numGCP)}
{'Number of tie/control points to be estimated'}	{num2str(data.numtie)}
{'Number of control/tie point unknowns'}	{num2str(data.numtie*3)}
{'\line'}	{''}
{'Total Unknowns'}	{num2str(length(xhat))}
{'\n'}	{''}
{'Number of image points'}	{num2str(data.n/2)}
{'Total number of observations'}	{num2str(data.n)}
{'Number of Inner Constraints'}	{num2str(7*data.settings.Inner_Constraints)}
{'\line'}	{''}
{'Total Number of Observations'}	{num2str(data.n + 7*data.settings.Inner_Constraints)}
{'\n'}	{''}
{'Total Degrees of Freedom'}	{num2str((data.n + 7*data.settings.Inner_Constraints) - length(xhat))}
{'\n'}	{''}
{'A-Posteriori'}	{num2str(sigma02,10)}
{'RMSx'}	{num2str(RMSx,10)}
{'RMSy'}	{num2str(RMSy,10)}
{'RMS'}	{num2str(RMS,10)}
{'\n'}	{''}
];

printCell(fileID, tmp, '', padding);

fprintf(fileID, [line '\n\n']);


decimals = '5'; % number of decimal places in .out file
% determine longest targetID
unique_IDs = unique({data.points.targetID});
target_width = 0; % width for targetIDs
for i = 1:length(unique_IDs)
    if length(unique_IDs{i}) > target_width
        target_width = length(unique_IDs{i});
    end
end
% determine longest imgID
unique_IDs = unique({data.points.imageID});
img_width = 0; % width for image IDs
for i = 1:length(unique_IDs)
    if length(unique_IDs{i}) > target_width
        img_width = length(unique_IDs{i});
    end
end
clear unique_IDs
num_width = 12; % width for numbers
width = num2str(max([target_width img_width num_width]) + 2);
clear target_width img_width num_width

% Estimated EOPs for each image
fprintf(fileID, 'Estimated EOPs\nEOP Name\tValue\tStandard Deviation\n');

xhat_count = 1;
EOP_IOP_Corr = cell(data.numImg,4); % cell matrix to store EOP/IOP Correlations
for i = 1:size(EXT,1) % for each image
    imageID = EXT{i,1};
    camID = EXT{i,2};
    EOP_IOP_Corr(i,2) = {camID};
    count = countImagePoints(imageID,data);
    fprintf(fileID,'\n');
    tmp = [{'Image'} {imageID}
        {'Camera'} {camID}
        {'Number of image points'} {num2str(count)}
        {'\line'} {''}];
    printCell(fileID, tmp, '', padding);
    
    % Xc
    if data.settings.Estimate_Xc
        printEOP(fileID,'Xc',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals); % print EOP value and std
        EOP_IOP_Corr{i,1} = [EOP_IOP_Corr{i,1} xhat_count]; % save EOP xhatindex for later
        EOP_IOP_Corr{i,3}{end + 1} = 'Xc'; % save EOP name for later
        xhat_count = xhat_count + 1; % increment xhat_count counter
    end
    % Yc
    if data.settings.Estimate_Yc
        printEOP(fileID,'Yc',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        EOP_IOP_Corr{i,1} = [EOP_IOP_Corr{i,1} xhat_count];
        EOP_IOP_Corr{i,3}{end + 1} = 'Yc';
        xhat_count = xhat_count + 1;
    end
    % Zc
    if data.settings.Estimate_Zc
        printEOP(fileID,'Zc',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        EOP_IOP_Corr{i,1} = [EOP_IOP_Corr{i,1} xhat_count];
        EOP_IOP_Corr{i,3}{end + 1} = 'Zc';
        xhat_count = xhat_count + 1;
    end
    % orientation angles need to be converted to degrees form radians
    % w
    if data.settings.Estimate_w
        printEOP(fileID,'Omega',xhat(xhat_count)*180/pi(),sqrt(Cx(xhat_count,xhat_count))*180/pi(),width,decimals);
        EOP_IOP_Corr{i,1} = [EOP_IOP_Corr{i,1} xhat_count];
        EOP_IOP_Corr{i,3}{end + 1} = 'Omega';
        xhat_count = xhat_count + 1;
    end
    % p
    if data.settings.Estimate_p
        printEOP(fileID,'Phi',xhat(xhat_count)*180/pi(),sqrt(Cx(xhat_count,xhat_count))*180/pi(),width,decimals);
        EOP_IOP_Corr{i,1} = [EOP_IOP_Corr{i,1} xhat_count];
        EOP_IOP_Corr{i,3}{end + 1} = 'Phi';
        xhat_count = xhat_count + 1;
    end
    % k
    if data.settings.Estimate_k
        printEOP(fileID,'Kappa',xhat(xhat_count)*180/pi(),sqrt(Cx(xhat_count,xhat_count))*180/pi(),width,decimals);
        EOP_IOP_Corr{i,1} = [EOP_IOP_Corr{i,1} xhat_count];
        EOP_IOP_Corr{i,3}{end + 1} = 'Kappa';
        xhat_count = xhat_count + 1;
    end
end

% Estimates for IOPs and distortions for each camera
fprintf(fileID,['\n' line '\n\nEstimated IOPs and Distortions for each Camera\nIOP Name\tValue\tStandard Deviation\n\n']);
PAR = [{'Created with Fish-eye model Bundle Adjustment version:'} {version} {[]}; % cell matrix to store camera IOPs and their standard deviations
    {'Execution date'} {date} {[]};
    {[]} {[]} {[]};];

for i = 1:size(INT,1)/2 % for each camera
    cameraID = INT{i*2-1,1};
    ydir = INT{i*2-1,2};
    xmin = INT{i*2-1,3};
    ymin = INT{i*2-1,4};
    xmax = INT{i*2-1,5};
    ymax = INT{i*2-1,6};
    tmp = [{'Camera'} {cameraID}
        {'y axis dir'} {num2str(ydir)}
        {'x min'} {num2str(xmin)}
        {'y min'} {num2str(ymin)}
        {'x max'} {num2str(xmax)}
        {'y max'} {num2str(ymax)}
        {'\line'} {''}];
%     Partmp = cell(1,3);
%     Partmp(1,1:2) = [{'Camera'} {cameraID}];
    PAR = [PAR; [{'Camera'} {cameraID} {[]}];];
    printCell(fileID, tmp, '', padding);
    xhat_count_start = xhat_count;
    if data.settings.Estimate_xp
        printEOP(fileID,'xp',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        PAR = [PAR; [{'xp'} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
        xhat_count = xhat_count + 1;
    end
    if data.settings.Estimate_yp
        printEOP(fileID,'yp',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        PAR = [PAR; [{'yp'} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
        xhat_count = xhat_count + 1;
    end
    if data.settings.Estimate_c
        printEOP(fileID,'c',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        PAR = [PAR; [{'c'} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
        xhat_count = xhat_count + 1;
    end
    if data.settings.Estimate_radial
        for j = 1:data.settings.Num_Radial_Distortions
            printDist(fileID,strcat('k',num2str(j)),xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
            PAR = [PAR; [{strcat('k',num2str(j))} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
            xhat_count = xhat_count + 1;
        end
    end
    if data.settings.Estimate_decent
        for j = 1:2 
            printDist(fileID,strcat('p',num2str(j)),xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
            PAR = [PAR; [{strcat('p',num2str(j))} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
            xhat_count = xhat_count + 1;
        end
    end
    xhat_count_end = xhat_count - 1;
    
    fprintf(fileID,'\nIOP Correlation sub-matrix\n-------------------------------\n');
    Corr_sub = Correlation(xhat_count_start:xhat_count_end,xhat_count_start:xhat_count_end); % get correlation sub-matrix
    num_IOPs = (xhat_count_end-xhat_count_start);
    names = PAR((i-1)*(num_IOPs+2)+5:(i-1)*(num_IOPs+2)+5+num_IOPs,1); % get IOP list from PAR
    names = [{''}; names;]; % add whitespace
    fprintf(fileID,strcat('%-6.2s'),names{:}); % print cell names at the top
    fprintf(fileID,'\n');
    
    % print sub-matrix
    for j = 1:size(Corr_sub,1)
        fprintf(fileID,strcat('%-6.2s'),names{j+1}); % print name
        for k=1:j
            fprintf(fileID,strcat('%-+6.2f'),Corr_sub(j,k)); % print only lower triangular part of matrix
        end
        fprintf(fileID,'\n'); % new line
    end
    fprintf(fileID,'\n'); % new line
    
    % loop through EOP_IOP_Corr array and assign correct IOP indices for each image
    for j = 1:size(EOP_IOP_Corr,1) % for each image
        if cameraID == EOP_IOP_Corr{j,2} % if cameraID matches
            % append IOP indices
            EOP_IOP_Corr{j,1} = [EOP_IOP_Corr{j,1} xhat_count_start:xhat_count_end];
            % append IOP names
            EOP_IOP_Corr{j,3} = [EOP_IOP_Corr{j,3} names{:}];
        end
        % get EOP/IOP Correlation sub-matrix
        num_EOPIOP = size(EOP_IOP_Corr{j,1},2);
        indices = EOP_IOP_Corr{j,1};
        sub_mat = zeros(num_EOPIOP);
        for k = 1:num_EOPIOP % for each row (EOP/IOP)
            for l = 1:k % for each col in the lower triangle
                sub_mat(k,l) = Correlation(indices(k), indices(l));
            end
        end
        EOP_IOP_Corr{j,4} = sub_mat;
    end
end

% Estimated Ground Coordinates
if data.settings.Estimate_tie
    fprintf(fileID,['\n' line '\n\nEstimated Ground Coordinates of targets\nTargetID\tnumImages\tX\tY\tZ\tstdX\tstdY\tstdZ\n\n']);
    all_var = [];
    for i = 1:size(TIE,1) % for each tie point/target
        targetID = TIE{i}; % get target name
        numImages = countTargetImages(targetID,data); % get number of images for target
        XYZ = xhat(xhat_count:xhat_count+2); % get estimated XYZ form xhat
        stdxyz = zeros(3,1);
        varxyz = zeros(1,3);
        for j = 1:3 % get estimated standard deviations of XYZ from Cxhat
            stdxyz(j) = sqrt(Cx(xhat_count+j-1,xhat_count+j-1));
            varxyz(j) = Cx(xhat_count+j-1,xhat_count+j-1);
        end
        all_var = [all_var; varxyz;];
        printTIE(fileID,targetID,numImages,XYZ,stdxyz,width,decimals); % print to output file
        xhat_count = xhat_count + 3;    
    end
    % compute mean std in X Y and Z
    avg_std = sqrt(mean(all_var));
    % print to file
    fprintf(fileID,'\n\t\tMeanStd X\tMeanStd Y\tMeanStd Z\n');
    fprintf(fileID,strcat('\t\t%1$-',width,'.',decimals,'f%2$-',width,'.',decimals,'f%3$-',width,'.',decimals,'f\n'),avg_std(1),avg_std(2),avg_std(3));
end

% corrected image measurements
fprintf(fileID,['\n' line '\n\nCorrected Image Measurements\nPointID\tImageID\tCorrected x\tCorrected y\n\n']);
for i = 1:data.n/2 % for each point
    fprintf(fileID, strcat('%1$-',width,'s%2$-',width,'s%3$-',width,'.',decimals,'f%4$-',width,'.',decimals,'f\n'),data.points(i).targetID,data.points(i).imageID,data.points(i).x_corr,data.points(i).y_corr);
end

if xhat_count ~= size(A,2) + 1
    disp("warning: xhat_count didn't end on it's expected value (unknowns + 1)");
end

% Calculate Mean correlation coefficients between EOPs and IOPs using EOP_IOP_Corr array
fprintf(fileID,['\n' line '\n\nAbsolute (positive) mean correlation coefficients between EOPs and IOPs\n\n']);
EOP_IOP_Corr = sortrows(EOP_IOP_Corr,2); % sort by camera name
count = 1;
while true % loop through each camera
    sum_count = 1;
    if count > size(EOP_IOP_Corr,1)
        break; % break when count exceeds EOP_IOP_Corr's rows
    end
    camID = EOP_IOP_Corr{count,2};
    fprintf(fileID, ['Camera ' EOP_IOP_Corr{count,2} '\n']);
    names = [{''}; EOP_IOP_Corr{count,3}';]; % add whitespace
    fprintf(fileID,'%-6.2s',names{:}); % print cell names at the top
    fprintf(fileID,'\n');
    
    mean_corr = abs(EOP_IOP_Corr{count,4});
    count = count + 1;
    while strcmp(EOP_IOP_Corr{count,2}, camID) % while the cameraID stays the same
        mean_corr = mean_corr + abs(EOP_IOP_Corr{count,4}); % sum the absolute value of the correlation matrices
        count = count + 1;
        sum_count = sum_count + 1;
        if count > size(EOP_IOP_Corr,1)
            break; % break when count exceeds EOP_IOP_Corr's rows
        end
    end
    % calculate mean correlation
    mean_corr = mean_corr./sum_count;
    % print mean correlation sub-matrix
    for j = 1:size(mean_corr,1)
        fprintf(fileID,strcat('%-6.2s'),names{j+1}); % print name
        for k=1:j
            fprintf(fileID,strcat('%-+6.2f'),mean_corr(j,k)); % print only lower triangular part of matrix
        end
        fprintf(fileID,'\n'); % new line
    end
    fprintf(fileID,'\n'); % new line
end
clear sum_count count    

% check point differences in cp_diff
if data.settings.Check_Points
    fprintf(fileID,['\n' line '\n\nCheck point differences\n']);
    fprintf(fileID, strcat('%1$-',width,'s%2$-',width,'s%3$-',width,'s%4$-',width,'s\n\n'),'TargetID', 'diff X', 'diff Y', 'diff Z');
    for i = 1:size(cp_diff,1) % for each point
        fprintf(fileID, strcat('%1$-',width,'s%2$-',width,'.',decimals,'f%3$-',width,'.',decimals,'f%4$-',width,'.',decimals,'f\n'),cp_diff{i,1},cp_diff{i,2},cp_diff{i,3},cp_diff{i,4});
    end
    % print mean and rms check point differences
    fprintf(fileID, strcat('\n%1$-',width,'s%2$-',width,'.',decimals,'f%3$-',width,'.',decimals,'f%4$-',width,'.',decimals,'f\n'),'Mean',cp_mean(1),cp_mean(2),cp_mean(3));
    fprintf(fileID, strcat('%1$-',width,'s%2$-',width,'.',decimals,'f%3$-',width,'.',decimals,'f%4$-',width,'.',decimals,'f\n'),'RMS',cp_rms(1),cp_rms(2),cp_rms(3));
end

% close file
fclose(fileID);

% output PAR and RSD files
%[~,name,~] = fileparts(data.settings.Output_Filename);
writecell(RSD,strcat(name,'.rsd'),'Delimiter','tab','FileType','text');
writecell(PAR,strcat(name,'.par'),'Delimiter','tab','FileType','text');

% change back to project directory if in batch mode
if batch
    cd(projectDir)
end

fclose('all'); % just in case files aren't closed properly
disp('Done!');

end

%% Functions
% small functions not worth putting in their own files
function printEOP(fileID,name,value,std,width,decimals)
fprintf(fileID, strcat('%1$-',width,'.',decimals,'s%2$-',width,'.',decimals,'f%3$-',width,'.',decimals,'f\n'),name,value,std);
end
function printDist(fileID,name,value,std,width,decimals)
fprintf(fileID, strcat('%1$-',width,'.',decimals,'s%2$-',width,'.',decimals,'e%3$-',width,'.',decimals,'e\n'),name,value,std);
end
function printTIE(fileID,targetID,numImages,XYZ,stdxyz,width,decimals)
fprintf(fileID, strcat('%1$-',width,'s%2$-',width,'.','0','i%3$-',width,'.',decimals,'f%4$-',width,'.',decimals,'f%5$-',width,'.',decimals,'f%6$-',width,'.',decimals,'f%7$-',width,'.',decimals,'f%8$-',width,'.',decimals,'f\n'),targetID,numImages,XYZ(1),XYZ(2),XYZ(3),stdxyz(1),stdxyz(2),stdxyz(3));
end
function count = countImagePoints(imageID,data)
count = 0;
for function_i=1:size(data.points,2)
    if strcmp(data.points(function_i).imageID,imageID)
        count = count + 1;
    end
end
end
function count = countTargetImages(targetID,data)
count = 0;
for function_i=1:size(data.points,2)
    if strcmp(data.points(function_i).targetID,targetID)
        count = count + 1;
    end
end
end
% function to calculate rms given a vector y
function out = rms(y)
    N = length(y);
    sumsq = sum(y.^2);
    out = sqrt(sumsq/N);
end