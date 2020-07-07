%% Build design matrix A, misclosure matrix w and inner constrains design matrix G
function [error, A, misclosure, G, dist_scaling] = BuildAwG(PHO,EXT,CNT,INT,TIE,xhat,... % data
    BuildG, Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c, Estimate_radial, Num_Radial_Distortions, Estimate_decent) % settings
%% Setup
error = 0;
Use_Collinearity = 0; % Functionality was added to use collinearity equations instead of fish-eye. This variable enables/disables that code
% Num_Radial_Distortions must always be at least 1
if Num_Radial_Distortions<1
    Num_Radial_Distortions = 1;
end

% A should be an nxu matrix, w is nx1
n = size(PHO,1)*2; % number of measurements
numimg = size(EXT,1);
numcam = size(INT,1)/2;
u = size(xhat,1); % number of unknowns
u_perimage = Estimate_Xc + Estimate_Yc + Estimate_Zc + Estimate_w + Estimate_p + Estimate_k; % unknowns per image
u_percam = Estimate_c + Estimate_xp + Estimate_yp + Estimate_radial*Num_Radial_Distortions + Estimate_decent*2; % unknowns per camera

% matrix dist_scaling is created to keep track of all the scale factors for the distortions of each camera
% dist_scaling = [(xhat index for radial)(xhat index for decentering)(scale1)(scale2)(scale3)...] <- Cam0
%                [(xhat index for radial)(xhat index for decentering)(scale1)(scale2)(scale3)...] <- Cam1
%                ...
% Where scale = rmax^2*i
dist_scaling = zeros(numcam,2+Num_Radial_Distortions);

if BuildG
    d = 7; %number of inner constraints (X, Y, Z, w, p, k, and scale)
    G = zeros(u,d);
else
    G = 0;
end

A = zeros(n,u);
misclosure = zeros(n,1);
images = {};

%% Main loop
% For every row of the image in PHO (which correlates to 2 rows of A and w)
for i=1:length(PHO)   
    % get variables for current rows
    targetID = PHO{i,1}; % target name/id
    imageID = PHO{i,2}; % image ID
    x_px = PHO{i,3}; % x coordinate on the image in pixels
    y_px = PHO{i,4}; % y coordinate on the image in pixels
    
    %% get EOPs from xhat or EXT depending on settings
    % find correct image in EXT
    ext_index = -1;
    for j = 1:length(EXT)
        if strcmp(EXT{j,1},imageID)
            ext_index = j;
            break;
        end        
    end
    if ext_index < 0
        errordlg(['Could not find image ' imageID ' from .pho in .ext. Check that the image ID exists in both files'],['Error on image ' imageID]);
        error = 1;
        return
    end
    % get cameraID in EXT
    cameraID = EXT{ext_index,2};
    
    % get index of xhat for this image
    xhat_index = ext_index*u_perimage-(u_perimage-1);
    
    xhat_count = 0;
    % XYZc
    if Estimate_Xc == 1
        Xc = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        Xc = EXT{ext_index,3};
    end    
    if Estimate_Yc == 1
        Yc = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        Yc = EXT{ext_index,4};
    end    
    if Estimate_Zc == 1
        Zc = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        Zc = EXT{ext_index,5};
    end
    
    % omega phi kappa
    if Estimate_w == 1
        w = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        w = EXT{ext_index,6};
    end    
    if Estimate_p == 1
        p = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        p = EXT{ext_index,7};
    end    
    if Estimate_k == 1
        k = xhat(xhat_index+xhat_count);
        %xhat_count = xhat_count + 1;
    else
        k = EXT{ext_index,8};
    end
    
    %% Get ground coordinates of object point
    % first check if targetID is a tie point
    tieIndex = -1;
    % find current targetID in TIE
    for j = 1:size(TIE,1)
        if strcmp(TIE(j),targetID)
            tieIndex = j;
            break;
        end
    end
    if tieIndex ~= -1 % if targetID was found in TIE, then it is a TIE point and must be estimated
        xhat_index_TIE = u_perimage*numimg + u_percam*numcam + tieIndex*3 - 2;
        % get object coordinates
        X = xhat(xhat_index_TIE);
        Y = xhat(xhat_index_TIE + 1);
        Z = xhat(xhat_index_TIE + 2);
    else % if targetID is not a tie point    
        % find correct target in CNT
        cnt_index = -1;
        for j = 1:length(CNT)
            if strcmp(CNT(j,1),targetID)
                cnt_index = j;
                break;
            end        
        end
        if cnt_index < 0
            errordlg(['Could not find target ' targetID ' from .pho in .cnt. Check that the target ID exists in both files'],['Error on target ' targetID]);
            error = 1;
            return
        end
        % get object coordinates
        X = CNT{cnt_index,2};
        Y = CNT{cnt_index,3};
        Z = CNT{cnt_index,4};
    end
    
    %% get IOPs and distortions from xhat or INT depending on settings
    % find correct camera in INT
    int_index = -1;
    for j = 1:2:length(INT)
        if strcmp(INT(j,1),cameraID)
            int_index = j;
            break;
        end        
    end
    if int_index < 0
        errordlg(['Could not find camera ' cameraID ' from .ext in .int. Check that the camera ID exists in both files'],['Error on camera ' cameraID]);
        error = 1;
        return
    end
    
    cam_num = (int_index + 1)/2; % camera number
    xhat_index_IOP = u_perimage*numimg + cam_num*u_percam-(u_percam-1); % rows in xhat where IOPs may be located
    xhat_count = 0;
    % IOPs
    if Estimate_xp == 1
        xp = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        xp = INT{int_index + 1,1};
    end
    if Estimate_yp == 1
        yp = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        yp = INT{int_index + 1,2};
    end
    if Estimate_c == 1
        c = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        c = INT{int_index + 1,3};
    end
    
    % distortions
    if Estimate_radial == 1
        % find radial distortions in xhat
        xhat_index_radial = xhat_index_IOP + xhat_count;
        K = [xhat(xhat_index_radial:xhat_index_radial+Num_Radial_Distortions-1)];
        % add xhat index to dist_scaling matrix
        dist_scaling(cam_num,1) = xhat_index_radial;
        
        xhat_count = xhat_count + Num_Radial_Distortions;
    else
        % get radial distortion from INT
        K = [INT{int_index + 1,4:4+Num_Radial_Distortions-1}]';
    end
    if Estimate_decent == 1
        % find decentering distortions in xhat
        xhat_index_decent = xhat_index_IOP + xhat_count;
        P = [xhat(xhat_index_decent:xhat_index_decent+1)];
        % add xhat index to dist_scaling matrix
        dist_scaling(cam_num,2) = xhat_index_decent;
        
        %xhat_count = xhat_count + 2;
    else    
        P = [INT{int_index + 1,4+Num_Radial_Distortions:5+Num_Radial_Distortions}]';
    end
    
    % get y_axis_dir
    y_dir = INT{int_index,2};
    % Check that y_dir is +-1
    if y_dir ~= 1 && y_dir ~= -1
        disp('y_dir should be +-1 only')
        error = 1;
        return
    end
       
    %% build rows of A for EOPs and distortions
    
    % UVW
    U = (Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc);
    V = (Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc);
    W = sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc);
    R = sqrt(U^2+V^2);
    % radial distortions
    x_bar = x_px-xp;
    y_bar = y_px-yp;
    r = sqrt((x_bar)^2 + (y_bar)^2);
    delta_r = 0;
    for j = 1:length(K)
        delta_r = delta_r + K(j)*r^(2*j);
    end
    %delta_r = K(1)*r^2 + K(2)*r^4 + K(3)*r^6 + K(4)*r^8 + K(5)*r^10;
    % decentering distortions
    if length(P) < 2 % if no decentering distortions are given, P is just 0
        P = [0 0]; % make it [0 0] to avoid errors
    end
    decentering_x = P(1)*(y_bar^2 + 3*x_bar^2) + 2*P(2)*x_bar*y_bar;
    decentering_y = P(2)*(x_bar^2 + 3*y_bar^2) + 2*P(1)*x_bar*y_bar;
    
    % fx fy
    if Use_Collinearity == 0
        % fish-eye model
        fx = -c*U/R*atan(R/W) + xp + delta_r*(x_bar) + decentering_x;
        fy = -c*y_dir*V/R*atan(R/W) + yp + delta_r*(y_bar) + decentering_y;
    elseif Use_Collinearity == 1
        % collinearity model
        fx = -c*U/W + xp + delta_r*(x_bar) + decentering_x;
        fy = -c*y_dir*V/W + yp + delta_r*(y_bar) + decentering_y;
    else
        % if Use_Collinearity has an invalid value
        errordlg('Use_Collinearity must be either 0 or 1');
        error = 1;
        return
    end
    
    % EOP partial derivatives (different for Collinearity vs Fisheye models)
    Ablock = zeros(2,u_perimage);
    count_A = 1;
    
    if Estimate_Xc
        if Use_Collinearity == 0% fish-eye model
            % dfx/dXc
            A14 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*cos(k)*cos(p))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*((sin(p)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
            % dfy/dXc
            A24 = - (c*y_dir*((sin(p)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*cos(p)*sin(k))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
        else % collinearity model
            A14 = (c*cos(k)*cos(p))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) - (c*sin(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
            A24 = - (c*y_dir*sin(p)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (c*y_dir*cos(p)*sin(k))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
        end            
        Ablock(:,count_A) = [A14; A24;];
        count_A = count_A + 1;
    end
    if Estimate_Yc
        if Use_Collinearity == 0% fish-eye model
            % dfx/dYc
            A15 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) + (c*((2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + (cos(p)*sin(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2));
            % dfy/dYc
            A25 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*y_dir*((2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + (cos(p)*sin(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        else % collinearity model
            A15 = (c*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) + (c*cos(p)*sin(w)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
            A25 = (c*y_dir*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) + (c*y_dir*cos(p)*sin(w)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
        end
        Ablock(:,count_A) = [A15; A25;];
        count_A = count_A + 1;
    end
    if Estimate_Zc
        if Use_Collinearity == 0% fish-eye model
            % dfx/dZc
            A16 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*((2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) - (cos(p)*cos(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
            % dfy/dZc
            A26 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*y_dir*((2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) - (cos(p)*cos(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        else % collinearity model
            A16 = (c*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) - (c*cos(p)*cos(w)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
            A26 = (c*y_dir*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) - (c*y_dir*cos(p)*cos(w)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;            
        end
        Ablock(:,count_A) = [A16; A26;];
        count_A = count_A + 1;
    end
    if Estimate_w
        if Use_Collinearity == 0% fish-eye model
            % dfx/dw
            A11 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*(((cos(p)*cos(w)*(Y - Yc) + cos(p)*sin(w)*(Z - Zc))*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
            % dfy/dw
            A21 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*y_dir*(((cos(p)*cos(w)*(Y - Yc) + cos(p)*sin(w)*(Z - Zc))*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        else % collinearity model
            A11 = (c*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) - (c*(cos(p)*cos(w)*(Y - Yc) + cos(p)*sin(w)*(Z - Zc))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
            A21 = (c*y_dir*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) - (c*y_dir*(cos(p)*cos(w)*(Y - Yc) + cos(p)*sin(w)*(Z - Zc))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
        end
        Ablock(:,count_A) = [A11; A21;];
        count_A = count_A + 1;
    end
    if Estimate_p
        if Use_Collinearity == 0% fish-eye model
            % dfx/dp
            A12 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*((2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + ((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(cos(p)*(X - Xc) - cos(w)*sin(p)*(Z - Zc) + sin(p)*sin(w)*(Y - Yc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
            % dfy/dp
            A22 = (c*y_dir*((2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + ((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(cos(p)*(X - Xc) - cos(w)*sin(p)*(Z - Zc) + sin(p)*sin(w)*(Y - Yc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2));
        else % collinearity model
            A12 = (c*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) + (c*(cos(p)*(X - Xc) - cos(w)*sin(p)*(Z - Zc) + sin(p)*sin(w)*(Y - Yc))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
            A22 = (c*y_dir*(cos(p)*(X - Xc) - cos(w)*sin(p)*(Z - Zc) + sin(p)*sin(w)*(Y - Yc))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (c*y_dir*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
        end
        Ablock(:,count_A) = [A12; A22;];
        count_A = count_A + 1;
    end
    if Estimate_k
        if Use_Collinearity == 0% fish-eye model
            % dfx/dk
            A13 = -(c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
            % dfy/dk
            A23 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
        else % collinearity model
            A13 = -(c*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
            A23 = (c*y_dir*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
        end
        Ablock(:,count_A) = [A13; A23;];
        %count_A = count_A + 1;
    end
           
    % determine where in the A matrix to put this Ablock
    % the row is determined by i, the column is determined by the image # and we need to keep track of which images have already been encountered
    Arow = 2*i-1; % PHO has x,y measurements on the same row, but A requires x and y on different rows 
    
    if isempty(images) % if images is empty
        images = {imageID};
    end
    Acol = -1;
    % find image in images vector
    for j = 1:length(images)
        if strcmp(images{j},imageID)
            Acol = size(Ablock,2)*j-size(Ablock,2)+1;
            break;
        end
    end
    % if image isn't already in images, then add it
    if Acol<0 
        images = [images; {imageID};];
        Acol = u_perimage*length(images) - u_perimage + 1;
    end
        
    if Acol<0 || Arow<0
        disp('this should never happen');
        error = 1;
        return;
    end
    % place Ablock so that the first element is at (Arow,Acol) in A
    A(Arow:Arow+1,Acol:Acol+u_perimage-1) = Ablock;
    
    %% Build rows of A for IOPs and distortions
    % xp, yp, and c for each camera
    Ablock_IOPs = zeros(2,u_percam);
    count_A = 1;
    
    % partial derivatives 
    if Estimate_xp % same for Fish-eye vs Collinearity
        Ablock_IOPs(:,count_A) = [1;0;];
        count_A = count_A + 1;
    end
    if Estimate_yp % same for Fish-eye vs Collinearity
        Ablock_IOPs(:,count_A) = [0;1;];
        count_A = count_A + 1;
    end
    if Estimate_c % different for Fish-eye vs Collinearity
        if Use_Collinearity == 0% fish-eye model
            Ax_c = -(atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
            Ay_c = -(y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
        else % collinearity model
            Ax_c = -((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
            Ay_c = -(y_dir*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
        end
        Ablock_IOPs(:,count_A) = [Ax_c;Ay_c;];
        count_A = count_A + 1;
    end
    
    % Distortion partial derivatives must be scaled to avoid poor conditioning of the normal equations matrix
    % calculate max possible radial distortion rmax
    xmin = INT{int_index,3};
    ymin = INT{int_index,4};
    xmax = INT{int_index,5};
    ymax = INT{int_index,6};
    rmax = sqrt(((xmax-xmin)*0.5)^2+((ymax-ymin)*0.5)^2);
    % calculate scale facors
    for j = 1:Num_Radial_Distortions
        dist_scaling(cam_num,2+j) = rmax^(2*j);
    end
    
    if Estimate_radial == 1 % same for Fish-eye vs Collinearity
        Ax_K = zeros(1,Num_Radial_Distortions);
        Ay_K = zeros(1,Num_Radial_Distortions);
        for j=1:Num_Radial_Distortions
            % calculate partial derivative with scale factor
            Ax_K(1,j) = r^(2*j)*x_bar/dist_scaling(cam_num,2+j);
            Ay_K(1,j) = r^(2*j)*y_bar/dist_scaling(cam_num,2+j);
        end        
        Ablock_IOPs(:,count_A:count_A+length(Ax_K)-1) = [Ax_K; Ay_K;];
        count_A = count_A + length(Ax_K);
    end
    if Estimate_decent == 1 % same for Fish-eye vs Collinearity
        % calculate partial derivative with scale factor
        Ax_P = [(y_bar^2 + 3*x_bar^2) 2*x_bar*y_bar]./dist_scaling(cam_num,3);
        Ay_P = [2*x_bar*y_bar (x_bar^2 + 3*y_bar^2)]./dist_scaling(cam_num,3);
        Ablock_IOPs(:,count_A:count_A+length(Ax_P)-1) = [Ax_P; Ay_P;];
        %count_A = count_A + length(Ax_P);
    end
    
    % Calculate Acol_IOP
    Acol_IOP = u_perimage*numimg + cam_num*u_percam - (u_percam - 1);
    
    % Place Ablock_IOP
    A(Arow:Arow+1,Acol_IOP:Acol_IOP+u_percam-1) = Ablock_IOPs;        
    
    %% Build Rows of A for TIE points (if needed)
    if tieIndex ~= -1 % if targetID was found in TIE, then it is a TIE point and must be estimated
        % partial derivatives for tie points (different for Fish-eye vs Collinearity)
        if Use_Collinearity == 0% fish-eye model
            dx_dX = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*cos(k)*cos(p))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) + (c*((sin(p)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
            dx_dY = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*((2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + (cos(p)*sin(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
            dx_dZ = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*((2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) - (cos(p)*cos(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));

            dy_dX = (c*y_dir*((sin(p)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) + (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*cos(p)*sin(k))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
            dy_dY = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*((2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + (cos(p)*sin(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
            dy_dZ = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*((2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) - (cos(p)*cos(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        else % collinearity model
            dx_dX = (c*sin(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (c*cos(k)*cos(p))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
            dx_dY = - (c*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) - (c*cos(p)*sin(w)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
            dx_dZ = (c*cos(p)*cos(w)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (c*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
            
            dy_dX = (c*y_dir*sin(p)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + (c*y_dir*cos(p)*sin(k))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
            dy_dY =  - (c*y_dir*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)) - (c*y_dir*cos(p)*sin(w)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2;
            dy_dZ = (c*y_dir*cos(p)*cos(w)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (c*y_dir*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc));
        end
        % build rows of A
        Ablock_TIE = [dx_dX dx_dY dx_dZ;
                    dy_dX dy_dY dy_dZ;];
        % place rows of A for tie points
        Acol_TIE = xhat_index_TIE;
        A(Arow:Arow+1,Acol_TIE:Acol_TIE+2) = Ablock_TIE;
    end
    
    %% Build w rows
    % w11 is fx - x_px
    w11 = fx - x_px;
    % w21 is fy - y_px
    w21 = fy - y_px;    
    % W is a vector, so the col is always 1, and the row is the same as Arow
     misclosure(Arow:Arow+1,1) = [w11;
                        w21;];
    %% Build G rows (same for Fish-eye vs Collinearity)
    if BuildG % if inner constraints are set in the settings
        if G(xhat_index,1) == 0 % if G hasn't already been built for these rows
            Gblock = [
                1 0 0 0 -Zc Yc Xc;
                0 1 0 Zc 0 -Xc Yc;
                0 0 1 -Yc Xc 0 Zc;
                0 0 0 -1 -sin(w)*tan(p) cos(w)*tan(p) 0;
                0 0 0 0 -cos(w) -sin(w) 0;
                0 0 0 0 sin(w)*sec(p) -cos(w)*sec(p) 0;
                ];
            % place Gblock in the correct rows of G
            G(xhat_index:xhat_index+5,:) = Gblock; % the rows of G correspond to the rows in xhat
        end
    end
end
end