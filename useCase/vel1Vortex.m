function vel = vel1Vortex(time, xy, T, n)
% return the velocity of Rider-Kothe reversed single vortex test
%  as in PAM 2008 JCP.
% Input
%  xy is a 2-by-n matrix denoting n points.
%  n  is the number of points for plotting the velocity field.

if nargin<3
    T = 2;
end
if nargin<4
    n = 0;
end

x = xy(1,:);
y = xy(2,:);
nP = max(size(x));
vel = zeros(2,nP);
vel(1,:) = sin(pi*x).*sin(pi*x).*sin(2*pi*y);
vel(2,:) =-sin(pi*y).*sin(pi*y).*sin(2*pi*x);
vel = vel*cos(pi/T*time);

% add a plot to visualize the velocity
if n>0
    [X,Y] = meshgrid(0:1/n:1, 0:1/n:1);
    U = sin(pi*X).*sin(pi*X).*sin(2*pi*Y);
    V =-sin(pi*Y).*sin(pi*Y).*sin(2*pi*X);
    quiver(X,Y,U,V);
end
