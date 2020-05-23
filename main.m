%% Main bundle adjustment script
% Wynand Tredoux, May 2020
clear all
close all
format longg
clc
%% read in files
[filereaderror, files] = ReadFiles({'.pho','.ext','.cnt','.int'});
if filereaderror == 1
    disp ('Error reading files')
    return
else
    disp('Files read successfully!')
end

PHO = files{1}; % image measurements [pointID image xmm ymm]
n = length(PHO)*2; % number of measurements, (x,y) for each line of PHO
EXT = files{2}; % EOPs [imageID CamreaID Xc Yc Zc omega phi kappa]
u = size(EXT,1)*6; % number of unknowns, (Xc, Yc, Zc, omega, phi, kappa) for each image in EXT

CNT = files{3}; % object coordinates [TargetID X Y Z]
INT = files{4}; % IOPs [CameraID yaxis_dir xmin ymin xmax ymax]
%                      [xp yp c]

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
    tmp(i+1,1) = {str2double(INT(i+1,1))};
    tmp(i+1,2) = {str2double(INT(i+1,2))};
    tmp(i+1,3) = {str2double(INT(i+1,3))};
end
INT = tmp;
clear tmp

% initialize xhat (initial values for unknowns [Xc Yc Zc omeaga phi kappa], from EXT)
[xhaterror, xhat] = Buildxhat(EXT);
if xhaterror == 1
    disp('Error building xhat');
    return
end

threshold = 0.00000001;
deltasum = 100;
count = 0;
%xhat_arr = xhat';
deltasumarr = [];
% main loop
while deltasum > threshold
    count = count + 1;
    disp(['Itteration ' num2str(count) ':']);
    %% build A and w matricies
    [Awerror, A, w] = BuildAw(PHO, EXT, CNT, INT, xhat);
    if Awerror == 1
        disp('Error building A and w');
        return
    end
    %disp('Aw Built!')
    
    %% Calculate solution
    u = A'*w;
    N = A'*A;
    delta = -(N)^-1*u;
    xhat = xhat + delta;
    %xhat_arr = [xhat_arr; xhat';];
        
    deltasum = sumabs(delta)
    deltasumarr = [deltasumarr deltasum];
    % residuals
    %v = A*delta + w;
    % constrain loop to 100 itterations
    if count > 100
        break;
    end
end

disp('Done!');
%xhat

figure;
hold on
plot(1:length(deltasumarr),deltasumarr)
title('normal of \delta over time')
xlabel('Itteration')
ylabel('normal value')

% figure;
% hold on
% title('(X_c,Y_c,Z_c) over time')
% xlabel('Itteration')
% ylabel('value')
% X = repmat([1:size(xhat_arr,1)]',[1, 3]);
% plot(X,xhat_arr(:,1:3))
% legend
% figure;
% hold on
% title('(\omega,\phi,\kappa) over time')
% xlabel('Itteration')
% ylabel('value')
% plot(X,xhat_arr(:,4:end))
% legend