function vel = velRotation(time, xy, T)
% return the velocity of a pure rotation

x = xy(:,1);
y = xy(:,2);
vel =[-y , x];

