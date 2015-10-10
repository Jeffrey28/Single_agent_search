% 7/11/14
% The function calculates the angle (in radian) of a vector
% 12/9/14
% the domain is [0,2*pi)
function angle = calAngle (vec)
len = norm(vec);
x = vec(1);
y = vec(2);
if x >= 0
    angle = asin(y/len);
else
    angle = pi - asin(y/len);
end
if angle < 0 
    if angle < -pi/2
        error('angle out of range, wrong design of the calAngel algorithm')
    end
    angle = angle+2*pi;
end
end