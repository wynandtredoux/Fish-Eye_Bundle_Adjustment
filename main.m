%% Bundle adjustment script using equidistant fish-eye model
% Wynand Tredoux, May 2020
clear all
close all
format longg
clc

tic %start time
date = char(datetime); %date
%% read in files
[filereaderror, files] = ReadFiles({'.pho','.ext','.cnt','.int','.cfg'});
if filereaderror == 1
    disp ('Error reading files')
    return
else
    disp('Files read successfully!')
end

PHO = files{1}; % image measurements [pointID image xmm ymm]
%n = length(PHO)*2; % number of measurements, (x,y) for each line of PHO
EXT = files{2}; % EOPs [imageID CamreaID Xc Yc Zc omega phi kappa]
%u = size(EXT,1)*6; % number of unknowns, (Xc, Yc, Zc, omega, phi, kappa) for each image in EXT

CNT = files{3}; % object coordinates [TargetID X Y Z]
INT = files{4}; % IOPs [CameraID yaxis_dir xmin ymin xmax ymax]
%                      [xp yp c]
TIE = []; % list of tie point target IDs (read in later if needed)

%% Get settings from cfg file
CFG = files{5};

% get output filename (is allowed to error without exiting the program)
cfg_errors = 0;
[Output_Filename,cfg_errors] = findSetting(CFG,'Output_Filename',cfg_errors);
if cfg_errors>0 % default filename if none is provided
    Output_Filename = 'output.out';
end

% get measurement standard deviation (is allowed to error without exiting the program)
cfg_errors = 0;
[Meas_std,cfg_errors] = findSetting(CFG,'Meas_std',cfg_errors);
if cfg_errors>0 % if no std is provided, set to 1
    Meas_std = 1;
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
    return
end

%% Read in tie point list if needed
if Estimate_tie == 1 && Estimate_AllGCP == 0
    [filereaderror, files] = ReadFiles({'.tie'});
    if filereaderror == 1
        disp ('Error reading files')
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
    for j = 1:(5+Num_Radial_Distortions)
        tmp(i+1,j) = {str2double(INT(i+1,j))};
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
    return
end


% build P weight matrix from Meas_std
Cl = diag(repmat(Meas_std,size(PHO,1)*2,1));
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
        k = delta_k(size(u,1)+1:end,:);
        % seperate Cx into only the unknowns part
        Cx = Cx(1:size(A,2),1:size(A,2));
    else
        Cx = N^-1;
        delta = -Cx*u;
    end
    
    % scale distortion parameters in delta with values from dist_scaling
    % for each camera
    for i = 1:size(dist_scaling,1)
        % radial distortion
        if Estimate_radial
            % get radial distortion index
            radial_index = dist_scaling(i,1);
            % scale paramaters
            for j = 1:Num_Radial_Distortions
                delta(radial_index+j-1) = delta(radial_index+j-1)/dist_scaling(i,j+2);
            end
        end
        % decentering distortion
        if Estimate_decent
            % get decentering distortion index
            decent_index = dist_scaling(i,2);
            % scale P1
            delta(decent_index) = delta(decent_index)/dist_scaling(i,3);
            % scale P2
            delta(decent_index+1) = delta(decent_index+1)/dist_scaling(i,3);
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
figure;
hold on
plot(1:length(deltasumarr),deltasumarr)
title('normal of \delta over time')
xlabel('Iteration')
ylabel('normal value')

%% Residuals
% calculate x,y residuals for all image measurements
v = A*delta + w;
% build RSD with the format point_id image _id x_obs y_obs radial_dist v_x v_y v_r v_t
RSD = BuildRSD(v,PHO,EXT,INT,xhat,... %data
    Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c, Estimate_radial, Num_Radial_Distortions, Estimate_decent); % settings

% create plot of radial component of the residuals - RSD(:,8) - as a function of radial distance - RSD(:,5)
%tmp = [[RSD{:,5}]' [RSD{:,8}]']; % get r and v_r from RSD
figure;
hold on
scatter([RSD{:,5}],[RSD{:,8}]);
title('radial component of the residuals v_r as a function of radial distance r')
xlabel('radial distance r')
ylabel('radial component of the residuals v_r')

%%  variance factor
sigma0 = v'*P*v/(size(A,2)-size(A,1))

%% Create output file
padding = 4;
disp('Writing output file...');
line = '*************************************************************************************************************';
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
    [gite, mfiles] = system('git ls-files -m');
end

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
tmp = [{'Iteration_Cap'} 	{num2str(Iteration_Cap)}
{'Threshold_Value'}	{num2str(threshold)}
{'Inner_Constraints'}	{num2str(Inner_Constraints)}
{'Estimate_Xc'}	{num2str(Estimate_Xc)}
{'Estimate_Yc'}	{num2str(Estimate_Yc)}
{'Estimate_Zc'}	{num2str(Estimate_Zc)}
{'Estimate_Omega'}	{num2str(Estimate_w)}
{'Estimate_Phi'}	{num2str(Estimate_p)}
{'Estimate_Kappa'}	{num2str(Estimate_k)}
{'Estimate_c'}	{num2str(Estimate_c)}
{'Estimate_xp'}	{num2str(Estimate_xp)}
{'Estimate_yp'}	{num2str(Estimate_yp)}
{'Estimate_Radial_Distortions'}	{num2str(Estimate_radial)}
{'Num_Radial_Distortions'}	{num2str(Num_Radial_Distortions)}
{'Estimate_tie'}	{num2str(Estimate_tie)}
{'Estimate_AllGCP'} {num2str(Estimate_AllGCP)}];

printCell(fileID, tmp, '\t\t', padding);
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
{'\n'}	{''}];

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
        printEOP(fileID,'Omeaga',xhat(xhat_count)*180/pi(),sqrt(Cx(xhat_count,xhat_count))*180/pi(),width,decimals);
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

for i = 1:size(INT,1)/2 % for each camera
    cameraID = INT{i*2-1,1};
    tmp = [{'Camera'} {cameraID}
        {'\line'} {''}];
    printCell(fileID, tmp, '', padding);
    if Estimate_xp
        printEOP(fileID,'xp',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        xhat_count = xhat_count + 1;
    end
    if Estimate_yp
        printEOP(fileID,'yp',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        xhat_count = xhat_count + 1;
    end
    if Estimate_c
        printEOP(fileID,'c',xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
        xhat_count = xhat_count + 1;
    end
    if Estimate_radial
        for j = 1:Num_Radial_Distortions
            printDist(fileID,strcat('k',num2str(j)),xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
            xhat_count = xhat_count + 1;
        end
    end
    if Estimate_decent
        for j = 1:2
            printDist(fileID,strcat('p',num2str(j)),xhat(xhat_count),sqrt(Cx(xhat_count,xhat_count)),width,decimals);
            xhat_count = xhat_count + 1;
        end
    end
end

% Estimated Ground Coordinates
fprintf(fileID,['\n' line '\n\nEstimated Ground Coordinates of targets\nTargetID\tX\tY\tZ\tstdX\tstdY\tstdZ\n\n']);
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
    
% close file
fclose(fileID);
disp('Done!');

% small functions not worth putting in their own files
function printEOP(fileID,name,value,std,width,decimals)
fprintf(fileID, strcat('%1$-',width,'.',decimals,'s%2$-',width,'.',decimals,'f%3$-',width,'.',decimals,'f\n'),name,value,std);
end
function printDist(fileID,name,value,std,width,decimals)
fprintf(fileID, strcat('%1$-',width,'.',decimals,'s%2$-',width,'.',decimals,'e%3$-',width,'.',decimals,'e\n'),name,value,std);
end
function printTIE(fileID,targetID,numImages,XYZ,stdxyz,width,decimals)
fprintf(fileID, strcat('%1$-',width,'.',decimals,'s%2$-',width,'.','0','i%3$-',width,'.',decimals,'f%4$-',width,'.',decimals,'f%5$-',width,'.',decimals,'f%6$-',width,'.',decimals,'f%7$-',width,'.',decimals,'f%8$-',width,'.',decimals,'f\n'),targetID,numImages,XYZ(1),XYZ(2),XYZ(3),stdxyz(1),stdxyz(2),stdxyz(3));
end
function count = countImagePoints(imageID,PHO)
count = 0;
for i=1:size(PHO,1)
    if strcmp(PHO{i,2},imageID)
        count = count + 1;
    end
end
end
function count = countTargetImages(targetID,PHO)
count = 0;
for i=1:size(PHO,1)
    if strcmp(PHO{i,1},targetID)
        count = count + 1;
    end
end
end