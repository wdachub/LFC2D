function [vel,M] = vel2VortexTimeDep(time, xy, T, n, nFrames)
% return the velocity of Rider-Kothe reversed single vortex test
%  modified to couple time and space
% n is the number of points for plotting the velocity field.

if nargin<3
    T = 2;
end
if nargin<4
    n = 0;
end
if nargin<5
    nFrames = 0;
end

x = xy(1,:);
y = xy(2,:);
nP = max(size(x));
vel = zeros(2,nP);
xx = x + time/T;
vel(1,:) = sin(pi*xx).*sin(pi*xx).*sin(2*pi*y);
vel(2,:) =-sin(pi*y) .*sin(pi*y) .*sin(2*pi*xx);

% add a plot to visualize the velocity
if n>0
    S=1;
    if nFrames==0
        [X,Y] = meshgrid(0:1/n:1, 0:1/n:1);
        XX = X + time/T;
%         U = cos(pi*time/T)*sin(pi*XX).*sin(pi*XX).*sin(2*pi*Y);
%         V =-cos(pi*time/T)*sin(pi*Y).*sin(pi*Y).*sin(2*pi*XX);
        U = sin(pi*XX).*sin(pi*XX).*sin(2*pi*Y);
        V =-sin(pi*Y).*sin(pi*Y).*sin(2*pi*XX);
        quiver(X,Y,U,V,S);
%         uv = sqrt(U.*U+V.*V);
%         maxVel = norm(uv(:), inf);
        axis equal;
%        title(['Max vel = ', num2str(maxVel)]);
    else
        for i=1:nFrames+1
            t = time;
            if nFrames>0
                t = t + (i-1)/nFrames*T;
            end
            [X,Y] = meshgrid(0:1/n:1, 0:1/n:1);
            XX = X + t/T;
            U = sin(pi*XX).*sin(pi*XX).*sin(2*pi*Y);
            V =-sin(pi*Y).*sin(pi*Y).*sin(2*pi*XX);
            quiver(X,Y,U,V,S);
%             uv = sqrt(U.*U+V.*V);
%             maxVel = norm(uv(:), inf);
            axis equal;
%             title(['Max vel = ', num2str(maxVel)]);
            M(i) = getframe(gcf);
        end
        movie(M,1);
    end
end
