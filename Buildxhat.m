% create the ux1 xhat vector from the .ext file
function [error, xhat] = Buildxhat(EXT, INT, ... % data
    Estimate_Xc, Estimate_Yc, Estimate_Zc, Estimate_w, Estimate_p, Estimate_k, Estimate_c, Estimate_xp, Estimate_yp) % settings
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

% count = 1;
% for i = 1:size(EXT,1)
%     for j = 3:8
%         xhat(count) = EXT{i,j};
%         count = count + 1;    
%     end
% end
end