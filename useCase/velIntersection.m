function vel = velIntersection(time, xy,T)
x = xy(:,1);
y = xy(:,2);
% this velocity field produces winding numbers +1,0,-1,-2
a=1;
b=2*pi;
vel = [a*x+b*y  -b*x+a*y];

