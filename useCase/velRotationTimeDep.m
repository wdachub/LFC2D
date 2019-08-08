function vel = velRotationTimeDep(time, xy, T)
% return the velocity of a pure rotation


x = xy(:,1);
y = xy(:,2);
vel = 1/(1+time.^2)*[-y , x];