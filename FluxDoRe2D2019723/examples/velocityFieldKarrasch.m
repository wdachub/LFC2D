function u = velocityFieldKarrasch(t,Pos)
% Definition of the velocity field.
% The velocity field must be defined for evaluation of multiple Position
% vectors at the given time. The following form must be satisfied
%
%          [x1         [vx(t,x1,y1)
%           x2          vx(t,x2,y2)         
%           ...         ...
%   Pos =   xn    ->    vx(t,xn,yn)      
%           y1          vy(t,x1,y1)
%           y2          vy(t,x2,y2)
%           ...         ...
%           yn]         vy(t,xn,yn)] 

half = numel(Pos)/2;

A = 0.25;


u = [-A*pi*sin(pi*Pos(1:half)).*cos(pi*Pos(half+1:end)) - t*sin(pi*t);
      A*pi*cos(pi*Pos(1:half)).*sin(pi*Pos(half+1:end)) - t*cos(pi*t)];


end