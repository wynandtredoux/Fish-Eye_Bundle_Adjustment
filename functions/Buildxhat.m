% create the ux1 xhat vector from the .ext file
function [error, xhat, xhatnames] = Buildxhat(data, EXT, INT, TIE, CNT)
error = 0;

% figure out how many unknowns based on settings
u = 0;
% Xc, Yc, Zc, omega, phi, and kappa add 1 unknown per image each
u = u + data.settings.Estimate_Xc*data.numImg + data.settings.Estimate_Yc*data.numImg + data.settings.Estimate_Zc*data.numImg;
u = u + data.settings.Estimate_w*data.numImg + data.settings.Estimate_p*data.numImg + data.settings.Estimate_k*data.numImg;
% c, xp, and yp add 1 unknown per camera each
u = u + data.settings.Estimate_c*data.numCam + data.settings.Estimate_xp*data.numCam + data.settings.Estimate_yp*data.numCam;
% radial distortion adds Num_Radial_Distortions unknowns per camera, and decentering distortion adds 2 per camera
u = u + data.settings.Estimate_radial*data.settings.Num_Radial_Distortions*data.numCam + data.settings.Estimate_decent*2*data.numCam;
% tie points add 3 unknowns per point (X, Y, Z)
u = u + size(TIE,1)*3;

xhat = zeros(u,1);
xhatnames = cell(u,1);
count = 1;

% image EOPs
for i = 1:data.numImg
    image = EXT{i,1};
    cam = EXT{i,2};
    Xc = EXT{i,3};
    Yc = EXT{i,4};
    Zc = EXT{i,5};
    w = EXT{i,6};
    p = EXT{i,7};
    k = EXT{i,8};
    
    if data.settings.Estimate_Xc
        xhat(count) = Xc;
        xhatnames(count) = {strcat('Xc_',image,'_',cam)};
        count = count + 1;
    end
    if data.settings.Estimate_Yc
        xhat(count) = Yc;
        xhatnames(count) = {strcat('Yc_',image,'_',cam)};
        count = count + 1;
    end
    if data.settings.Estimate_Zc
        xhat(count) = Zc;
        xhatnames(count) = {strcat('Zc_',image,'_',cam)};
        count = count + 1;
    end
    if data.settings.Estimate_w
        xhat(count) = w;
        xhatnames(count) = {strcat('w_',image,'_',cam)};
        count = count + 1;
    end
    if data.settings.Estimate_p
        xhat(count) = p;
        xhatnames(count) = {strcat('p_',image,'_',cam)};
        count = count + 1;
    end
    if data.settings.Estimate_k
        xhat(count) = k;
        xhatnames(count) = {strcat('k_',image,'_',cam)};
        count = count + 1;
    end
end

% camera IOPs and distortions
for i=1:data.numCam
    xp = INT{i*2,1};
    yp = INT{i*2,2};
    c = INT{i*2,3};
    cam = INT{i*2-1,1};
    
    k = [INT{i*2,4:4+data.settings.Num_Radial_Distortions-1}]';
    p = [INT{i*2,4+data.settings.Num_Radial_Distortions:5+data.settings.Num_Radial_Distortions}]';
    
    % IOPs
    if data.settings.Estimate_xp
        xhat(count) = xp;
        xhatnames(count) = {strcat('xp_','_',cam)};
        count = count + 1;
    end
    if data.settings.Estimate_yp
        xhat(count) = yp;
        xhatnames(count) = {strcat('yp_','_',cam)};
        count = count + 1;
    end
    if data.settings.Estimate_c
        xhat(count) = c;
        xhatnames(count) = {strcat('c_','_',cam)};
        count = count + 1;
    end
    % Distortions
    if data.settings.Estimate_radial
        for j=1:size(k,1)
            xhat(count) = k(j);
            xhatnames(count) = {strcat('k',num2str(j),'_',cam)};
            count = count + 1;
        end
    end
    if data.settings.Estimate_decent
        for j=1:size(p,1)
            xhat(count) = p(j);
            xhatnames(count) = {strcat('p',num2str(j),'_',cam)};
            count = count + 1;
        end
    end    
end

% Tie points
if size(TIE,1) > 0 % if TIE is not empty
    for i=1:size(TIE,1)
        % get tie point name from TIE
        targetID = TIE(i);
        % find initial coordinates from CNT
        X = [];
        Y = [];
        Z = [];
        for j=1:size(CNT,1)
            if strcmp(targetID, CNT{j,1})
                X = CNT{j,2};
                Y = CNT{j,3};
                Z = CNT{j,4};
                break;
            end
        end
        % if targetID's coordinates can't be found in CNT
        if size(X,1) == 0 || size(Y,1) == 0 || size(Z,1) == 0
            disp(['Error Buildxhat(): can''t find ' targetID ' from .tie in .cnt']);
            error = 1;
            return;
        end
        % add coordinates to xhat
        xhat(count:count+2) = [X;Y;Z;];
        xhatnames(count:count+2) = [{strcat('X_',targetID)}; {strcat('Y_',targetID)}; {strcat('Z_',targetID)};];
        count = count + 3;
    end
end
end