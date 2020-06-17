function RSD = BuildRSD(v,PHO,EXT,INT,xhat,... %data
    Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c, Estimate_radial, Num_Radial_Distortions, Estimate_decent) % settings
% calculate r,t residuals for all image measurements
RSD = cell(size(PHO,1),9);
RSD(:,1:4) = PHO;

for i = 1:size(RSD,1)
    % get vx,vy from v
    vx = v(2*i-1);
    vy = v(2*i);
    %% get (xp,yp) for this image (different for each camera)
    imageID = RSD{i,2};
    % find correct image in EXT
    ext_index = -1;
    for j = 1:length(EXT)
        if strcmp(EXT{j,1},imageID)
            ext_index = j;
            break;
        end        
    end
    % get cameraID in EXT
    cameraID = EXT{ext_index,2};

    int_index = -1;
    % find correct camera in INT
    for j = 1:2:length(INT)
        if strcmp(INT(j,1),cameraID)
            int_index = j;
            break;
        end        
    end
    cam_num = (int_index + 1)/2; % camera number
    
    % get xp,yp from xhat or INT depending on settings
    numimg = size(EXT,1); % total number of images
    u_perimage = Estimate_Xc + Estimate_Yc + Estimate_Zc + Estimate_w + Estimate_p + Estimate_k; % unknowns per image
    u_percam = Estimate_c + Estimate_xp + Estimate_yp + Estimate_radial*Num_Radial_Distortions + Estimate_decent*2; % unknowns per camera
    xhat_index_IOP = u_perimage*numimg + cam_num*u_percam-(u_percam-1); % rows in xhat where IOPs may be located
    xhat_count = 0; % xhat_count = 1 to skip c in xhat if Estimate_c == 1
    if Estimate_xp == 1
        xp = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        xp = INT{int_index + 1,1};
    end
    if Estimate_yp == 1
        yp = xhat(xhat_index_IOP+xhat_count);
        %xhat_count = xhat_count + 1;
    else    
        yp = INT{int_index + 1,2};
    end
    
    %% calculate vr and vt
    xbar = RSD{i,3} - xp;
    ybar = RSD{i,4} - yp;
    theta = atan2(ybar,xbar);
    Phi = atan2(vy,vx);
    v_dist = sqrt(vx^2 + vy^2);
    vr = v_dist*cos(theta-Phi);
    vt = v_dist*sin(theta-Phi);
    
    %% Fill RSD
    r = sqrt(xbar^2+ybar^2);
    RSD(i,5:9) = [{r} {vx} {vy} {vr} {vt}];    
    
end
end