function vel = velStrainTimeDep(time, xy, T)
% return the velocity of a pure strain

x = xy(:,1);
y = xy(:,2);
T = 1;
vel = [x, -y]*time/T;
