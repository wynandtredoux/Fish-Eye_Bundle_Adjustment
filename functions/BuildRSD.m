function RSD = BuildRSD(v, data, xhat)
u_perimage = data.settings.Estimate_Xc + data.settings.Estimate_Yc + data.settings.Estimate_Zc + data.settings.Estimate_w + data.settings.Estimate_p + data.settings.Estimate_k; % unknowns per image
u_percam = data.settings.Estimate_c + data.settings.Estimate_xp + data.settings.Estimate_yp + data.settings.Estimate_radial*data.settings.Num_Radial_Distortions + data.settings.Estimate_decent*2; % unknowns per camera
% calculate r,t residuals for all image measurements
RSD = cell(data.n/2,9);
RSD(:,1:4) = [{data.points.targetID}' {data.points.imageID}' {data.points.x}' {data.points.y}'];


for i = 1:size(RSD,1)
    % get vx,vy from v
    vx = v(2*i-1);
    vy = v(2*i);
    %% get xp,yp from xhat or data depending on settings
    xhat_index_IOP = u_perimage*data.numImg + data.points(i).cam_num*u_percam-(u_percam-1); % rows in xhat where IOPs may be located
    xhat_count = 0; % xhat_count = 1 to skip c in xhat if Estimate_c == 1
    if data.settings.Estimate_xp == 1
        xp = xhat(xhat_index_IOP+xhat_count);
        xhat_count = xhat_count + 1;
    else    
        xp = data.points(i).xp;
    end
    if data.settings.Estimate_yp == 1
        yp = xhat(xhat_index_IOP+xhat_count);
        %xhat_count = xhat_count + 1;
    else    
        yp = data.points(i).yp;
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