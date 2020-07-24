%% Build design matrix A, misclosure matrix w and inner constrains design matrix G
% data - contains points and settings structures
%       points - structure containing all info about each image point:
%           [x y targetID imageID ext_index cameraID Xc Yc Zc omeaga phi kappa int_index cam_num xp yp c K(radial distortion) P(decentering distortion) y_dir]
%       settings - contains all values from .cfg file
%           BuildG - build interior orientation constraints matrix
%           Estimate... - Toggle estimation of Xc,Yc,Zc,omeaga,phi,kappa,xp,yp,c,radial distortion parameters, and decentering distortion parameters
%           Num_Radial_Distortions - number of radial distortion parameters
%           type - model to be used, valid values are 'Fisheye', 'Pinhole'
% CNT - Object Point Coordinates [PointID  X  Y  Z]
% TIE - List of tie points [PointID] (can be empty if there are no tie points)
% xhat - estimated unknowns [Unknown_value] has length of u

function [error, A, misclosure, G, dist_scaling] = BuildAwG(data, xhat) % settings
%% Setup
error = 0;
% Num_Radial_Distortions must always be at least 1
if data.settings.Num_Radial_Distortions<1
    data.settings.Num_Radial_Distortions = 1;
end

% A should be an nxu matrix, w is nx1
u = size(xhat,1); % number of unknowns
u_perimage = data.settings.Estimate_Xc + data.settings.Estimate_Yc + data.settings.Estimate_Zc + data.settings.Estimate_w + data.settings.Estimate_p + data.settings.Estimate_k; % unknowns per image
u_percam = data.settings.Estimate_c + data.settings.Estimate_xp + data.settings.Estimate_yp + data.settings.Estimate_radial*data.settings.Num_Radial_Distortions + data.settings.Estimate_decent*2; % unknowns per camera

% matrix dist_scaling is created to keep track of all the scale factors for the distortions of each camera
% dist_scaling = [(xhat index for radial)(xhat index for decentering)(scale1)(scale2)(scale3)...] <- Cam0
%                [(xhat index for radial)(xhat index for decentering)(scale1)(scale2)(scale3)...] <- Cam1
%                ...
% Where scale = rmax^2*i
dist_scaling = zeros(data.numCam,2+data.settings.Num_Radial_Distortions);

if data.settings.Inner_Constraints
    d = 7; %number of inner constraints (X, Y, Z, w, p, k, and scale)
    G = zeros(u,d);
else
    G = 0;
end

A = zeros(data.n,u);
misclosure = zeros(data.n,1);
images = {};

%% Main loop
% For every row of the image in PHO (which correlates to 2 rows of A and w)
for i=1:data.n/2
    x = data.points(i).x; % x coordinate on the image
    y = data.points(i).y; % y coordinate on the image
    
    %% get EOPs from xhat or data.points depending on settings 
    % get index of xhat for this image
    xhat_index = data.points(i).ext_index*u_perimage-(u_perimage-1);
    
    xhat_count = 0;
    % XYZc
    if data.settings.Estimate_Xc == 1
        Xc = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        Xc = data.points(i).Xc;
    end    
    if data.settings.Estimate_Yc == 1
        Yc = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        Yc = data.points(i).Yc;
    end    
    if data.settings.Estimate_Zc == 1
        Zc = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        Zc = data.points(i).Zc;
    end
    
    % omega phi kappa
    if data.settings.Estimate_w == 1
        w = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        w = data.points(i).w;
    end    
    if data.settings.Estimate_p == 1
        p = xhat(xhat_index+xhat_count);
        xhat_count = xhat_count + 1;
    else
        p = data.points(i).p;
    end    
    if data.settings.Estimate_k == 1
        k = xhat(xhat_index+xhat_count);
        %xhat_count = xhat_count + 1;
    else
        k = data.points(i).k;
    end
    
    %% Get ground coordinates of point from xhat or data.points
    if data.points(i).isTie % if target is a tie point
        xhat_index_TIE = u_perimage*data.numImg + u_percam*data.numCam + data.points(i).tieIndex*3 - 2;
        % get object coordinates
        X = xhat(xhat_index_TIE);
        Y = xhat(xhat_index_TIE + 1);
        Z = xhat(xhat_index_TIE + 2);
    else % if targetID is not a tie point    
        % get object coordinates
        X = data.points(i).X;
        Y = data.points(i).Y;
        Z = data.points(i).Z;
    end
    
    %% get IOPs and distortions from xhat or data.points depending on settings
    xhat_index_IOP = u_perimage*data.numImg + data.points(i).cam_num*u_percam-(u_percam-1); % rows in xhat where IOPs may be located
    xhat_count = 0;
    % IOPs
    if data.settings.Estimate_xp == 1
        xp = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        xp = data.points(i).xp;
    end
    if data.settings.Estimate_yp == 1
        yp = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        yp = data.points(i).yp;
    end
    if data.settings.Estimate_c == 1
        c = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        c = data.points(i).c;
    end
    
    % distortions
    if data.settings.Estimate_radial == 1
        % find radial distortions in xhat
        xhat_index_radial = xhat_index_IOP + xhat_count;
        K = [xhat(xhat_index_radial:xhat_index_radial+data.settings.Num_Radial_Distortions-1)];
        % add xhat index to dist_scaling matrix
        dist_scaling(data.points(i).cam_num,1) = xhat_index_radial;
        
        xhat_count = xhat_count + data.settings.Num_Radial_Distortions;
    else
        % get radial distortion from INT
        K = data.points(i).K;
    end
    if data.settings.Estimate_decent == 1
        % find decentering distortions in xhat
        xhat_index_decent = xhat_index_IOP + xhat_count;
        P = [xhat(xhat_index_decent:xhat_index_decent+1)];
        % add xhat index to dist_scaling matrix
        dist_scaling(data.points(i).cam_num,2) = xhat_index_decent;
        
        %xhat_count = xhat_count + 2;
    else    
        P = data.points(i).P;
    end
    
    % get y_axis_dir
    y_dir = data.points(i).y_dir;
       
    %% build rows of A for EOPs and distortions
    
    % UVW
    U = (Y - Yc)*(cos(w)*sin(k) + cos(k)*sin(p)*sin(w)) + (Z - Zc)*(sin(k)*sin(w) - cos(k)*cos(w)*sin(p)) + cos(k)*cos(p)*(X - Xc);
    V = (Y - Yc)*(cos(k)*cos(w) - sin(k)*sin(p)*sin(w)) + (Z - Zc)*(cos(k)*sin(w) + cos(w)*sin(k)*sin(p)) - cos(p)*sin(k)*(X - Xc);
    W = sin(p)*(X - Xc) + cos(p)*cos(w)*(Z - Zc) - cos(p)*sin(w)*(Y - Yc);
    R = sqrt(U^2+V^2);
    % radial distortions
    x_bar = x-xp;
    y_bar = y-yp;
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
    if strcmp(data.settings.type,'fisheye')
        % fish-eye model
        fx = -c*U/R*atan(R/W) + xp + delta_r*(x_bar) + decentering_x;
        fy = -c*y_dir*V/R*atan(R/W) + yp + delta_r*(y_bar) + decentering_y;
        typeint = 0;
    elseif strcmp(data.settings.type,'pinhole')
        % collinearity model
        fx = -c*U/W + xp + delta_r*(x_bar) + decentering_x;
        fy = -c*y_dir*V/W + yp + delta_r*(y_bar) + decentering_y;
        typeint = 1;
    else
        % if type has an invalid value
        errordlg('BuildAwG, invalid type in data.settings.type');
        error = 1;
        return
    end
    
    % EOP partial derivatives (different for Collinearity vs Fisheye models)
    Ablock = zeros(2,u_perimage);
    count_A = 1;
    
    if data.settings.Estimate_Xc
        if typeint == 0% fish-eye model
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
    if data.settings.Estimate_Yc
        if typeint == 0% fish-eye model
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
    if data.settings.Estimate_Zc
        if typeint == 0% fish-eye model
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
    if data.settings.Estimate_w
        if typeint == 0% fish-eye model
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
    if data.settings.Estimate_p
        if typeint == 0% fish-eye model
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
    if data.settings.Estimate_k
        if typeint == 0% fish-eye model
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
    % the row is determined by i, the column is determined by the image #
    Arow = 2*i-1; % PHO has x,y measurements on the same row, but A requires x and y on different rows 
    Acol = size(Ablock,2)*data.points(i).ext_index-size(Ablock,2)+1;

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
    if data.settings.Estimate_xp % same for Fish-eye vs Collinearity
        Ablock_IOPs(:,count_A) = [1;0;];
        count_A = count_A + 1;
    end
    if data.settings.Estimate_yp % same for Fish-eye vs Collinearity
        Ablock_IOPs(:,count_A) = [0;1;];
        count_A = count_A + 1;
    end
    if data.settings.Estimate_c % different for Fish-eye vs Collinearity
        if typeint == 0% fish-eye model
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
    rmax = sqrt(((data.points(i).xmax-data.points(i).xmin)*0.5)^2+((data.points(i).ymax-data.points(i).ymin)*0.5)^2);
    % calculate scale facors
    for j = 1:data.settings.Num_Radial_Distortions
        dist_scaling(data.points(i).cam_num,2+j) = rmax^(2*j);
    end
    
    if data.settings.Estimate_radial == 1 % same for Fish-eye vs Collinearity
        Ax_K = zeros(1,data.settings.Num_Radial_Distortions);
        Ay_K = zeros(1,data.settings.Num_Radial_Distortions);
        for j=1:data.settings.Num_Radial_Distortions
            % calculate partial derivative with scale factor
            Ax_K(1,j) = r^(2*j)*x_bar/dist_scaling(data.points(i).cam_num,2+j);
            Ay_K(1,j) = r^(2*j)*y_bar/dist_scaling(data.points(i).cam_num,2+j);
        end        
        Ablock_IOPs(:,count_A:count_A+length(Ax_K)-1) = [Ax_K; Ay_K;];
        count_A = count_A + length(Ax_K);
    end
    if data.settings.Estimate_decent == 1 % same for Fish-eye vs Collinearity
        % calculate partial derivative with scale factor
        Ax_P = [(y_bar^2 + 3*x_bar^2) 2*x_bar*y_bar]./dist_scaling(data.points(i).cam_num,3);
        Ay_P = [2*x_bar*y_bar (x_bar^2 + 3*y_bar^2)]./dist_scaling(data.points(i).cam_num,3);
        Ablock_IOPs(:,count_A:count_A+length(Ax_P)-1) = [Ax_P; Ay_P;];
        %count_A = count_A + length(Ax_P);
    end
    
    % Calculate Acol_IOP
    Acol_IOP = u_perimage*data.numImg + data.points(i).cam_num*u_percam - (u_percam - 1);
    
    % Place Ablock_IOP
    A(Arow:Arow+1,Acol_IOP:Acol_IOP+u_percam-1) = Ablock_IOPs;        
    
    %% Build Rows of A for TIE points (if needed)
    if data.points(i).isTie % if targetID was found in TIE, then it is a TIE point and must be estimated
        % partial derivatives for tie points (different for Fish-eye vs Collinearity)
        if typeint == 0% fish-eye model
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
    % w11 is fx - x
    w11 = fx - x;
    % w21 is fy - y
    w21 = fy - y;    
    % W is a vector, so the col is always 1, and the row is the same as Arow
     misclosure(Arow:Arow+1,1) = [w11;
                        w21;];
    %% Build G rows (same for Fish-eye vs Collinearity)
    if data.settings.Inner_Constraints % if inner constraints are set in the settings
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