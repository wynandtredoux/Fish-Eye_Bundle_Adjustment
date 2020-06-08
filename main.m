%% Bundle adjustment script using equidistant fish-eye model
% Wynand Tredoux, May 2020
clear all
close all
format longg
clc
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

if cfg_errors>0
    disp ('Error getting settings')
    return
end

%% Read in tie point list if needed
if Estimate_tie
    [filereaderror, files] = ReadFiles({'.tie'});
    if filereaderror == 1
        disp ('Error reading files')
        return
    end
    TIE = files{1};
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
    Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_c, Estimate_xp, Estimate_yp, Estimate_radial, Num_Radial_Distortions, Estimate_decent); % settings
if xhaterror == 1
    disp('Error building xhat');
    return
end

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
    u = A'*w;
    N = A'*A;
    
    if Inner_Constraints
        NG = [N G;
            G' zeros(size(G,2));];
        uG = [u; zeros(size(G,2),1);];

        delta_k = -NG^-1*uG;
        delta = delta_k(1:size(u,1),:);
        k = delta_k(size(u,1)+1:end,:);
    else  
        delta = -(N)^-1*u;
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

disp('Done!');
%xhat

figure;
hold on
plot(1:length(deltasumarr),deltasumarr)
title('normal of \delta over time')
xlabel('Iteration')
ylabel('normal value')

% figure;
% hold on
% title('(X_c,Y_c,Z_c) over time')
% xlabel('Iteration')
% ylabel('value')
% X = repmat([1:size(xhat_arr,1)]',[1, 3]);
% plot(X,xhat_arr(:,1:3))
% legend
% figure;
% hold on
% title('(\omega,\phi,\kappa) over time')
% xlabel('Iteration')
% ylabel('value')
% plot(X,xhat_arr(:,4:end))
% legend