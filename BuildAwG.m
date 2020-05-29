%% Build design matrix A, misclosure matrix w and inner constrains design matrix G
function [error, A, misclosure, G] = BuildAwG(PHO,EXT,CNT,INT,xhat,... % data
    BuildG, Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c) % settings
error = 0;
% A should be an nxu matrix, w is nx1
n = size(PHO,1)*2; % number of measurements
numimg = size(EXT,1);
numcam = size(INT,1)/2;
u = size(xhat,1); % number of unknowns
u_perimage = Estimate_Xc + Estimate_Yc + Estimate_Zc + Estimate_w + Estimate_p + Estimate_k;
u_percam = Estimate_c + Estimate_xp + Estimate_yp;

if BuildG
    d = 7; %number of inner constraints (X, Y, Z, w, p, k, and scale)
    G = zeros(u,d);
else
    G = 0;
end

A = zeros(n,u);
misclosure = zeros(n,1);
images = {};
%% EOPs
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
    
    %% get IOPs from xhat or INT depending on settings
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
        %xhat_count = xhat_count + 1;
    else    
        c = INT{int_index + 1,3};
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
    % fx fy
    fx = -c*U/R*atan(R/W) + xp;
    fy = -c*y_dir*V/R*atan(R/W) + yp;
    
    % partial derivatives
    Ablock = zeros(2,u_perimage);
    count_A = 1;
    
    % El-Sheimy's notes have the A matrix in the order [w p k Xc Yc Zc] but we need it in [Xc Yc Zc w p k], so the order is changed
%     Ablock = [A14 A15 A16 A11 A12 A13;
%                 A24 A25 A26 A21 A22 A23;];
    
    if Estimate_Xc
        % dfx/dXc
        A14 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*cos(k)*cos(p))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*((sin(p)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        % dfy/dXc
        A24 = - (c*y_dir*((sin(p)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*cos(k)*cos(p)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) - 2*cos(p)*sin(k)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*cos(p)*sin(k))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
        Ablock(:,count_A) = [A14; A24;];
        count_A = count_A + 1;
    end
    if Estimate_Yc
        % dfx/dYc
        A15 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) + (c*((2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + (cos(p)*sin(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2));
        % dfy/dYc
        A25 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*y_dir*((2*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + (cos(p)*sin(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        Ablock(:,count_A) = [A15; A25;];
        count_A = count_A + 1;
    end
    if Estimate_Zc
        % dfx/dZc
        A16 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*((2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) - (cos(p)*cos(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        % dfy/dZc
        A26 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*y_dir*((2*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) - (cos(p)*cos(w)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        Ablock(:,count_A) = [A16; A26;];
        count_A = count_A + 1;
    end
    if Estimate_w
        % dfx/dw
        A11 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w))))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*(((cos(p)*cos(w)*(Y - Yc) + cos(p)*sin(w)*(Z - Zc))*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        % dfy/dw
        A21 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w))))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) - (c*y_dir*(((cos(p)*cos(w)*(Y - Yc) + cos(p)*sin(w)*(Z - Zc))*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 - (2*((Y - Yc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) - (Z - Zc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)) + 2*((Y - Yc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - (Z - Zc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        Ablock(:,count_A) = [A11; A21;];
        count_A = count_A + 1;
    end
    if Estimate_p
        % dfx/dp
        A12 = (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2)) + (c*((2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + ((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(cos(p)*(X - Xc) - cos(w)*sin(p)*(Z - Zc) + sin(p)*sin(w)*(Y - Yc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2));
        % dfy/dp
        A22 = (c*y_dir*((2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))) + ((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)*(cos(p)*(X - Xc) - cos(w)*sin(p)*(Z - Zc) + sin(p)*sin(w)*(Y - Yc)))/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2)*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc))^2 + 1)*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2) - (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*(2*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))*(cos(k)*sin(p)*(X - Xc) + cos(k)*cos(p)*cos(w)*(Z - Zc) - cos(k)*cos(p)*sin(w)*(Y - Yc)) - 2*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))*(sin(k)*sin(p)*(X - Xc) + cos(p)*cos(w)*sin(k)*(Z - Zc) - cos(p)*sin(k)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(2*(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(3/2));
        Ablock(:,count_A) = [A12; A22;];
        count_A = count_A + 1;
    end
    if Estimate_k
        % dfx/dk
        A13 = -(c*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
        % dfy/dk
        A23 = (c*y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
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
    
    %% Build rows of A for IOPs
    % xp, yp, and c for each camera
    Ablock_IOPs = zeros(2,u_percam);
    count_A = 1;
    
    % partial derivatives
    if Estimate_xp
        Ablock_IOPs(:,count_A) = [1;0;];
        count_A = count_A + 1;
    end
    if Estimate_yp
        Ablock_IOPs(:,count_A) = [0;1;];
        count_A = count_A + 1;
    end
    if Estimate_c
        Ax_c = -(atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
        Ay_c = -(y_dir*atan((((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2)/(sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc)))*((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc)))/(((Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc))^2 + ((Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc))^2)^(1/2);
        Ablock_IOPs(:,count_A) = [Ax_c;Ay_c;];
        %count_A = count_A + 1;
    end
    
    % Calculate Acol_IOP
    Acol_IOP = u_perimage*numimg + cam_num*u_percam - (u_percam - 1);
    
    % Place Ablock_IOP
    A(Arow:Arow+1,Acol_IOP:Acol_IOP+u_percam-1) = Ablock_IOPs;        
    
    %% Build w rows
    % w11 is fx - x_px
    w11 = fx - x_px;
    % w21 is fy - y_px
    w21 = fy - y_px;    
    % W is a vector, so the col is always 1, and the row is the same as Arow
     misclosure(Arow:Arow+1,1) = [w11;
                        w21;];
    %% Build G rows
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