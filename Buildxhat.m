% create the ux1 xhat vector from the .ext file
function [error, xhat] = Buildxhat(EXT)
error = 0;
if size(EXT,2)~=8
    error = 1;
    errordlg('Error! .ext file is in the wrong format. It should be in the format [ImageID CameraID Xc Yc Zc omega phi cappa]','Formatting error');
    return
end

xhat = zeros(size(EXT,1)*6,1);
count = 1;
for i = 1:size(EXT,1)
    for j = 3:8
        xhat(count) = EXT{i,j};
        count = count + 1;    
    end
end
end