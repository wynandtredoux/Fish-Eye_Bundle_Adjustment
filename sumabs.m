% find the absolute sum of a vector
function sum = sumabs(vect)
sum = 0;
if size(vect,2)~=1
    disp('Error: sumabs(vect), vect must be a vector')
    return
end
if size(vect,1)==0
    return
end

for i = 1:length(vect)
    sum = sum + abs(vect(i));
end
end