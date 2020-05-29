% create the ux1 xhat vector from the .ext file
function [error, xhat] = Buildxhat(EXT, INT, ... % data
    Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_xp, Estimate_yp, Estimate_c, Estimate_radial, Estimate_decent) % settings
error = 0;
numImg = size(EXT,1); % number of images
numCam = size(INT,1)/2;% number of cameras

% figure out how many unknowns based on settings
u = 0;
% Xc, Yc, Zc, omega, phi, and kappa add 1 unknown per image each
u = u + Estimate_Xc*numImg + Estimate_Yc*numImg + Estimate_Zc*numImg;
u = u + Estimate_w*numImg + Estimate_p*numImg + Estimate_k*numImg;
% c, xp, and yp add 1 unknown per camera each
u = u + Estimate_c*numCam + Estimate_xp*numCam + Estimate_yp*numCam;
% radial distortion adds 5 unknowns per camera, and decentering distortion adds 2 per camera
u = u + Estimate_radial*5*numCam + Estimate_decent*2*numCam;

xhat = zeros(u,1);
count = 1;

% image EOPs
for i = 1:numImg
    Xc = EXT{i,3};
    Yc = EXT{i,4};
    Zc = EXT{i,5};
    w = EXT{i,6};
    p = EXT{i,7};
    k = EXT{i,8};
    
    if Estimate_Xc
        xhat(count) = Xc;
        count = count + 1;
    end
    if Estimate_Yc
        xhat(count) = Yc;
        count = count + 1;
    end
    if Estimate_Zc
        xhat(count) = Zc;
        count = count + 1;
    end
    if Estimate_w
        xhat(count) = w;
        count = count + 1;
    end
    if Estimate_p
        xhat(count) = p;
        count = count + 1;
    end
    if Estimate_k
        xhat(count) = k;
        count = count + 1;
    end
end

% camera IOPs and distortions
for i=1:numCam
    xp = INT{i*2,1};
    yp = INT{i*2,2};
    c = INT{i*2,3};
    
    k = [INT{i*2,4:8}]';
    p = [INT{i*2,9:10}]';
    
    % IOPs
    if Estimate_xp
        xhat(count) = xp;
        count = count + 1;
    end
    if Estimate_yp
        xhat(count) = yp;
        count = count + 1;
    end
    if Estimate_c
        xhat(count) = c;
        count = count + 1;
    end
    % Distortions
    if Estimate_radial
        for j=1:size(k,1)
            xhat(count) = k(j);
            count = count + 1;
        end
    end
    if Estimate_decent
        for j=1:size(p,1)
            xhat(count) = p(j);
            count = count + 1;
        end
    end    
end
end