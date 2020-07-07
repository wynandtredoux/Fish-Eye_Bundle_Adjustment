%% Bundle adjustment script using equidistant fish-eye model
% Wynand Tredoux, May 2020

function main_error = main(varargin)
%% Setup
%clear all
%close all
format longg
%clc
celldisp(varargin)

main_error = 0;
enable_plots = 0;
batch = 0;
folder = '';

% get data directory if provided
if ~isempty(varargin)
    if length(varargin) > 1
        disp('Error: too many arguments in main');
        main_error = 1;
        return
    else
        batch = 1;
        folder = varargin{1};
    end
end
        
projectDir = pwd;

tic %start time
date = char(datetime); %date
% get version of code using the "git describe" command
[gite, version] = system('git describe --dirty'); % "dirty" indicates that the local version has been modified and doesn't fully match the version on github
if gite>0
    disp('Error: could not get git verion with "git describe"')
    disp('Version will be set to "unknown"');
    version = 'Unknown';
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

%% Get settings from cfg file
% get output filename from cfg (is allowed to error without exiting the program)
cfg_errors = 0;
[Output_Filename,cfg_errors] = findSetting(CFG,'Output_Filename',cfg_errors);
if cfg_errors>0 % default filename if none is provided
    Output_Filename = 'output.out';
end

% get measurement standard deviation (is allowed to error without exiting the program)
cfg_errors = 0;
no_std_y = 0;
[Meas_std,cfg_errors] = findSetting(CFG,'Meas_std',cfg_errors);
if cfg_errors>0 % if no std is provided, set to 1
    Meas_std = 1;
    no_std_y = 1;
else
    [Meas_std_y,no_std_y] = findSetting(CFG,'Meas_std_y',cfg_errors);
end

cfg_errors = 0;
% General Settings
[Iteration_Cap,cfg_errors] = findSetting(CFG,'Iteration_Cap',cfg_errors);
[threshold,cfg_errors] = findSetting(CFG,'Threshold_Value',cfg_errors);
% inner constraints
[Inner_Constraints,cfg_errors] = findSetting(CFG,'Inner_Constraints',cfg_errors,1);
% Estimate EOPs
[Estimate_Xc,cfg_errors] = findSetting(CFG,'Estimate_Xc',cfg_errors,1);
[Estimate_Yc,cfg_errors] = findSetting(CFG,'Estimate_Yc',cfg_errors,1);
[Estimate_Zc,cfg_errors] = findSetting(CFG,'Estimate_Zc',cfg_errors,1);
[Estimate_w,cfg_errors] = findSetting(CFG,'Estimate_Omega',cfg_errors,1);
[Estimate_p,cfg_errors] = findSetting(CFG,'Estimate_Phi',cfg_errors,1);
[Estimate_k,cfg_errors] = findSetting(CFG,'Estimate_Kappa',cfg_errors,1);
% Estimate IOPs
[Estimate_c,cfg_errors] = findSetting(CFG,'Estimate_c',cfg_errors,1);
[Estimate_xp,cfg_errors] = findSetting(CFG,'Estimate_xp',cfg_errors,1);
[Estimate_yp,cfg_errors] = findSetting(CFG,'Estimate_yp',cfg_errors,1);
% Estimate Distortions
[Estimate_radial,cfg_errors] = findSetting(CFG,'Estimate_Radial_Distortions',cfg_errors,1);
[Num_Radial_Distortions,cfg_errors] = findSetting(CFG,'Num_Radial_Distortions',cfg_errors);
[Estimate_decent,cfg_errors] = findSetting(CFG,'Estimate_Decentering_Distortions',cfg_errors,1);
% Estimate Ground Coordinates of tie points
[Estimate_tie,cfg_errors] = findSetting(CFG,'Estimate_tie',cfg_errors,1);
[Estimate_AllGCP,cfg_errors] = findSetting(CFG,'Estimate_AllGCP',cfg_errors,1);

if cfg_errors>0
    disp ('Error getting settings')
    main_error = 1;
    return
end

%% Read in tie point list if needed
if Estimate_tie == 1 && Estimate_AllGCP == 0
    [filereaderror, files] = ReadFiles({'.tie'});
    if filereaderror == 1
        disp ('Error reading files')
        main_error = 1;
        return
    end
    TIE = files{1};
end

%% Estimate_AllGCP
if Estimate_AllGCP == 1
    TIE = CNT(:,1); % add all targets as TIE points
    Estimate_tie = 1; % change Estimate_tie to 1 so all GCPs are estiamted
end

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
    for j = 4:(5+Num_Radial_Distortions)
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
%% Main Loop
% initialize xhat (initial values for unknowns)
[xhaterror, xhat] = Buildxhat(EXT, INT, TIE, CNT, ... % data
    Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c, Estimate_radial, Num_Radial_Distortions, Estimate_decent); % settings
if xhaterror == 1
    disp('Error building xhat');
    main_error = 1;
    return
end


% build P weight matrix from Meas_std and Meas_std_y if provided
if no_std_y
    Cl = diag(repmat(Meas_std^2,size(PHO,1)*2,1));
else
    Cl = diag(repmat([Meas_std^2; Meas_std_y^2;],size(PHO,1),1));
end
P = Cl^-1;

deltasum = 100;
count = 0;
xhat_arr = xhat';
deltasumarr = [];
% main loop
while deltasum > threshold
    count = count + 1;
    disp(['Iteration ' num2str(count) ':']);
    %% build A and w matrices
    [Awerror, A, w, G, dist_scaling] = BuildAwG(PHO, EXT, CNT, INT, TIE, xhat,...
        Inner_Constraints, Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c, Estimate_radial, Num_Radial_Distortions, Estimate_decent);
    if Awerror == 1
        disp('Error building A and w');
        main_error = 1;
        return
    end
    
    %% Calculate solution
    u = A'*P*w;
    N = A'*P*A;
    %Cx = [];
    
    if Inner_Constraints
        NG = [N G;
            G' zeros(size(G,2));];
        uG = [u; zeros(size(G,2),1);];
        Cx = NG^-1;
        
        % seperate k and delta
        delta_k = -Cx*uG;
        delta = delta_k(1:size(u,1),:);
        %k = delta_k(size(u,1)+1:end,:);
        
        % seperate Cx into only the unknowns part
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
        if Estimate_radial
            % get radial distortion index
            radial_index = dist_scaling(i,1);
            % scale paramaters
            for j = 1:Num_Radial_Distortions
                delta(radial_index+j-1) = delta(radial_index+j-1)/dist_scaling(i,j+2);
                Cx(radial_index+j-1,radial_index+j-1) = Cx(radial_index+j-1,radial_index+j-1)/dist_scaling(i,j+2);
            end
        end
        % decentering distortion
        if Estimate_decent
            % get decentering distortion index
            decent_index = dist_scaling(i,2);
            % scale P1
            delta(decent_index) = delta(decent_index)/dist_scaling(i,3);
            Cx(decent_index,decent_index) = Cx(decent_index,decent_index)/dist_scaling(i,3);
            % scale P2
            delta(decent_index+1) = delta(decent_index+1)/dist_scaling(i,3);
            Cx(decent_index+1,decent_index+1) = Cx(decent_index+1,decent_index+1)/dist_scaling(i,3);
        end
    end
    
    xhat = xhat + delta;
    xhat_arr = [xhat_arr; xhat';];
        
    deltasum = sumabs(delta)
    deltasumarr = [deltasumarr deltasum];
    % residuals
    %v = A*delta + w;
    % constrain loop iterations
    if count >= Iteration_Cap
        disp('Iteration Cap reached. This can be changed in the .cfg file')
        break;
    end
end
% end time
time = num2str(toc);
disp(['Elapsed time is ' char(time) ' seconds.'])



% plot normal of delta over time
if enable_plots
    figure;
    hold on
    plot(1:length(deltasumarr),deltasumarr)
    title('normal of \delta over time')
    xlabel('Iteration')
    ylabel('normal value')
end

%% Residuals
% calculate x,y residuals for all image measurements
v = A*delta + w;
% build RSD and corrected measurements with the format point_id image _id x_obs y_obs radial_dist v_x v_y v_r v_t
RSD = BuildRSD(v,PHO,EXT,INT,xhat,... %data
    Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c, Estimate_radial, Num_Radial_Distortions, Estimate_decent); % settings

% create plot of radial component of the residuals - RSD(:,8) - as a function of radial distance - RSD(:,5)
if enable_plots
    figure;
    hold on
    scatter([RSD{:,5}],[RSD{:,8}]);
    title('radial component of the residuals v_r as a function of radial distance r')
    xlabel('radial distance r')
    ylabel('radial component of the residuals v_r')
end

%% corrected image coordinates
PHO_corr = PHO; % get origional coordinates
for i = 1:length(PHO_corr)
    PHO_corr(i,3) = {PHO_corr{i,3} + RSD{i,6}}; % correct x coords
    PHO_corr(i,4) = {PHO_corr{i,4} + RSD{i,7}}; % correct y coords
end

%% RMS
% RMSx/RMSy
sum2x = 0;
sum2y = 0;
n = size(v,1)/2;
for i = 1:n
    sum2x = sum2x + (v(2*i-1))^2;
    sum2y = sum2y + (v(2*i))^2;
end
RMSx = sqrt(1/n*sum2x)
RMSy = sqrt(1/n*sum2y)
RMS = sqrt(RMSx^2 + RMSy^2)

%%  variance factor
sigma0 = sqrt(v'*P*v/(size(A,1)-size(A,2)))
%Cx = sigma0.*Cx;

%% Create output file
padding = 4;
disp('Writing output file...');
line = '*************************************************************************************************************';

% create output file
fileID = fopen(Output_Filename,'w');

% heading
fprintf(fileID,['Version: ' version]);
fprintf(fileID,'Fish-eye model Bundle Adjustment\nWynand Tredoux -- University of Calgary -- 2020\n\n');
if ~isempty(mfiles) % if modified files exist
    fprintf(fileID,['Modified files:\n' mfiles]); % print list of modified files
end
fprintf(fileID,line);
fprintf(fileID,['\n\nExecution date:\t' date '\nTime Taken:\t\t' char(time) ' seconds\nIterations:\t\t' num2str(count)]);

% settings used
fprintf(fileID,'\n\nSettings used:\n');
% tmp = [{'Iteration_Cap'} 	{num2str(Iteration_Cap)}
% {'Threshold_Value'}	{num2str(threshold)}
% {'Inner_Constraints'}	{num2str(Inner_Constraints)}
% {'Estimate_Xc'}	{num2str(Estimate_Xc)}
% {'Estimate_Yc'}	{num2str(Estimate_Yc)}
% {'Estimate_Zc'}	{num2str(Estimate_Zc)}
% {'Estimate_Omega'}	{num2str(Estimate_w)}
% {'Estimate_Phi'}	{num2str(Estimate_p)}
% {'Estimate_Kappa'}	{num2str(Estimate_k)}
% {'Estimate_c'}	{num2str(Estimate_c)}
% {'Estimate_xp'}	{num2str(Estimate_xp)}
% {'Estimate_yp'}	{num2str(Estimate_yp)}
% {'Estimate_Radial_Distortions'}	{num2str(Estimate_radial)}
% {'Num_Radial_Distortions'}	{num2str(Num_Radial_Distortions)}
% {'Estimate_tie'}	{num2str(Estimate_tie)}
% {'Estimate_AllGCP'} {num2str(Estimate_AllGCP)}];

printCell(fileID, CFG, '\t\t', padding);
fprintf(fileID, ['\n' line '\n']);

% Summery of unknowns/observations
fprintf(fileID, '\nObservations/Unknowns Summery\n\n');

tmp = [{'Number of Photos'}	{num2str(size(EXT,1))}
{'Total EOP unknowns'}	{num2str((Estimate_Xc + Estimate_Yc + Estimate_Zc + Estimate_w + Estimate_p + Estimate_k)*size(EXT,1))}
{'Number of Cameras'}	{num2str(size(INT,1)/2)}
{'Total IOP unknowns'}	{num2str((Estimate_c + Estimate_xp + Estimate_yp)*size(INT,1)/2)}
{'Total distortion unknowns'}	{num2str((Estimate_radial*Num_Radial_Distortions + Estimate_decent*2)*size(INT,1)/2)}
{'Number of control points'}	{num2str(size(CNT,1))}
{'Number of tie/control points to be estimated'}	{num2str(size(TIE,1))}
{'Number of control/tie point unknowns'}	{num2str(size(TIE,1)*3)}
{'\line'}	{''}
{'Total Unknowns'}	{num2str(length(xhat))}
{'\n'}	{''}
{'Number of image points'}	{num2str(size(PHO,1))}
{'Total number of observations'}	{num2str(size(PHO,1)*2)}
{'Number of Inner Constraints'}	{num2str(7*Inner_Constraints)}
{'\line'}	{''}
{'Total Number of Observations'}	{num2str(size(PHO,1)*2 + 7*Inner_Constraints)}
{'\n'}	{''}
{'Total Degrees of Freedom'}	{num2str((size(PHO,1)*2 + 7*Inner_Constraints) - length(xhat))}
{'\n'}	{''}
{'Sigma0'}	{num2str(sigma0,10)}
{'RMSx'}	{num2str(RMSx,10)}
{'RMSy'}	{num2str(RMSy,10)}
{'RMS'}	{num2str(RMS,10)}
{'\n'}	{''}
];

printCell(fileID, tmp, '', padding);

fprintf(fileID, [line '\n\n']);

% Estimated EOPs for each image
decimals = '5';
width = '14';
fprintf(fileID, 'Estimated EOPs\nEOP Name\tValue\tStandard Deviation\n');
%u_perimage = Estimate_Xc + Estimate_Yc + Estimate_Zc + Estimate_w + Estimate_p + Estimate_k; %unknowns per image

xhat_count = 1;
for i = 1:size(EXT,1) % for each image
    imageID = EXT{i,1};
    count = countImagePoints(imageID,PHO);
    fprintf(fileID,'\n');
    tmp = [{'Image'} {imageID}
        {'Camera'} {EXT{i,2}}
        {'Number of image points'} {num2str(count)}
        {'\line'} {''}];
    printCell(fileID, tmp, '', padding);
    
    % Xc
    if Estimate_Xc
        %fprintf(fileID, strcat('%1$-',width,'.',decimals,'s%2$-',width,'.',decimals,'f%3$-',width,'.',decimals,'f\n'),'Xc',xhat(xhat_count),Cx(xhat_count,xhat_count));
        printEOP(fileID,'Xc',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        xhat_count = xhat_count + 1;
    end
    % Yc
    if Estimate_Yc
        printEOP(fileID,'Yc',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        xhat_count = xhat_count + 1;
    end
    % Xc
    if Estimate_Zc
        printEOP(fileID,'Zc',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        xhat_count = xhat_count + 1;
    end
    % orentation angles need to be converted to degrees form radians
    % w
    if Estimate_w
        printEOP(fileID,'Omega',xhat(xhat_count)*180/pi(),sqrt(Cx(xhat_count,xhat_count))*180/pi(),width,decimals);
        xhat_count = xhat_count + 1;
    end
    % p
    if Estimate_p
        printEOP(fileID,'Phi',xhat(xhat_count)*180/pi(),sqrt(Cx(xhat_count,xhat_count))*180/pi(),width,decimals);
        xhat_count = xhat_count + 1;
    end
    % k
    if Estimate_k
        printEOP(fileID,'Kappa',xhat(xhat_count)*180/pi(),sqrt(Cx(xhat_count,xhat_count))*180/pi(),width,decimals);
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
    if Estimate_xp
        printEOP(fileID,'xp',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        PAR = [PAR; [{'xp'} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
        xhat_count = xhat_count + 1;
    end
    if Estimate_yp
        printEOP(fileID,'yp',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        PAR = [PAR; [{'yp'} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
        xhat_count = xhat_count + 1;
    end
    if Estimate_c
        printEOP(fileID,'c',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        PAR = [PAR; [{'c'} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
        xhat_count = xhat_count + 1;
    end
    if Estimate_radial
        for j = 1:Num_Radial_Distortions
            printDist(fileID,strcat('k',num2str(j)),xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
            PAR = [PAR; [{strcat('k',num2str(j))} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
            xhat_count = xhat_count + 1;
        end
    end
    if Estimate_decent
        for j = 1:2 
            printDist(fileID,strcat('p',num2str(j)),xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
            PAR = [PAR; [{strcat('p',num2str(j))} {xhat(xhat_count)} {sqrt(Cx(xhat_count,xhat_count))}];];
            xhat_count = xhat_count + 1;
        end
    end
    xhat_count_end = xhat_count - 1;
    fprintf(fileID,'\nCovariance sub-matrix\n-------------------------------\n');
    Corr_sub = Correlation(xhat_count_start:xhat_count_end,xhat_count_start:xhat_count_end); % get correlation submatrix
    names = PAR(i*9-4:i*9-4+(xhat_count_end-xhat_count_start),1); % get IOP list from PAR
    names = [{''}; names;]; % add whitespace
    fprintf(fileID,strcat('%-6.2s'),names{:}); % print cell names at the top
    fprintf(fileID,'\n');
    
    % print submatrix
    for j = 1:size(Corr_sub,1)
        fprintf(fileID,strcat('%-6.2s'),names{j+1}); % print name
        for k=1:j
            fprintf(fileID,strcat('%-+6.2f'),Corr_sub(j,k)); % print only lower trianglular part of matrix
        end
        fprintf(fileID,'\n'); % new line
    end
end

% Estimated Ground Coordinates
fprintf(fileID,['\n' line '\n\nEstimated Ground Coordinates of targets\nTargetID\tnumImages\tX\tY\tZ\tstdX\tstdY\tstdZ\n\n']);
for i = 1:size(TIE,1) % for each tie point/target
    targetID = TIE(i); % get target name
    numImages = countTargetImages(targetID,PHO); % get number of images for target
    XYZ = xhat(xhat_count:xhat_count+2); % get estimated XYZ form xhat
    stdxyz = zeros(3,1);
    for j = 1:3 % get estimated standard deviations of XYZ from Cxhat
        stdxyz(j) = sqrt(Cx(xhat_count+j-1,xhat_count+j-1));
    end
    printTIE(fileID,targetID,numImages,XYZ,stdxyz,width,decimals); % print to output file
    xhat_count = xhat_count + 3;    
end

% corrected image measurements
fprintf(fileID,['\n' line '\n\nCorrected Image Measurements\nPointID\tImageID\tCorrected x\tCorrected y\n\n']);
for i = 1:size(PHO_corr,1) % for each point
    fprintf(fileID, strcat('%1$-',width,'s%2$-',width,'s%3$-',width,'.',decimals,'f%4$-',width,'.',decimals,'f\n'),PHO_corr{i,1},PHO_corr{i,2},PHO_corr{i,3},PHO_corr{i,4});
end


if xhat_count ~= size(A,2) + 1
    disp("warning: xhat_count didn't end on it's expected value (unknowns + 1)");
end

% close file
fclose(fileID);

% output PAR and RSD files
[~,name,~] = fileparts(Output_Filename);
writecell(RSD,strcat(name,'.rsd'),'Delimiter','tab','FileType','text');
writecell(PAR,strcat(name,'.par'),'Delimiter','tab','FileType','text');

% change back to project directory if in batch mode
if batch
    cd(projectDir)
end

disp('Done!');

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
function count = countImagePoints(imageID,PHO)
count = 0;
for function_i=1:size(PHO,1)
    if strcmp(PHO{function_i,2},imageID)
        count = count + 1;
    end
end
end
function count = countTargetImages(targetID,PHO)
count = 0;
for function_i=1:size(PHO,1)
    if strcmp(PHO{function_i,1},targetID)
        count = count + 1;
    end
end
end

end