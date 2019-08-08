function vel = velDeformation(time, xy, T)
% return the velocity of the deformation test.

if nargin<3
    T = 2;
end

n = 4;
x = n*pi*(xy(:,1)+0.5);
y = n*pi*(xy(:,2)+0.5);

vel= -cos(pi*time/T)*[sin(x).*sin(y) cos(x).*cos(y)];
