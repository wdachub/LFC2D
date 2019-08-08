function vel = velStrain(time, xy, T)
% return the velocity of a pure strain

x = xy(:,1);
y = xy(:,2);
vel = [x, -y];
