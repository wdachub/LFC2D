function vel = vel13RotationCompress(time, xy, T)
% return the velocity of a compressed rotation

b = 0.1;
x = xy(1,:);
y = xy(2,:);
vel = [-y-b*x; x-b*y]*(1+0.1*cos(time*pi));

