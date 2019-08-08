function testVelocityField(v,C)
% The function tests if the velocity field u is defined properly, this is
% if the calculation at time t and the positions (x1,y1) ,... (xn,yn) is 
% defined as
%
%          [x1         [vx(t,x1,y1)
%           x2          vx(t,x2,y2)         
%           ...         ...
%           xn    ->    vx(t,xn,yn)      
%           y1          vy(t,x1,y1)
%           y2          vy(t,x2,y2)
%           ...         ...
%           yn]         vy(t,xn,yn)] 
%
% Therfore random points of the surface C are taken and the velocity at
% those points is calculated in both ways, point by point and all points at
% once, for a random time value. If the results do not match an error
% message is displayed
%
% List of Parameters:
%
% =========================================================================
%   INPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
%      C      | Defines the surface, stored as 
%             |    C = [x1 y1; 
%             |         x2 y2; 
%             |          ... 
%             |         xn yn]
%             | The number of points that are stored defines the spatial
%             | resolution (Each point is advected back and the
%             | back-advected points define the backadvected surface.
%      -------|------------------------------------------------------------
%      v      | Either a anonymous function or the handle to a m-Function
%             | in the current folder. It must be defined for evaluation 
%             | of multiple Position vectors at the given time. The 
%             | following form must be satisfied
%             |
%             |                [x1         [vx(t,x1,y1)
%             |                 x2          vx(t,x2,y2)         
%             |                 ...         ...
%             |         Pos =   xn    ->    vx(t,xn,yn)      
%             |                 y1          vy(t,x1,y1)
%             |                 y2          vy(t,x2,y2)
%             |                 ...         ...
%             |                 yn]         vy(t,xn,yn)] 


%% Choose random points
nOfRandPoints = 3;

randIndices = round((size(C,1)-1).*rand(nOfRandPoints,1))+1;
Points = C(randIndices,:);

%% Create random time value
tVal = rand(1);

%% Calc velocity values point by point
uPbyP = zeros(nOfRandPoints,2);
for i=1:nOfRandPoints
    uPbyP(i,:) = v(tVal,Points(i,:)')';
end
uPbyP = uPbyP(:);

%Calc velocity valures all at once
uAao = v(tVal,Points(:));

if ~isempty(find(uPbyP-uAao,1))
    error('The velocity field is not defined correctly!');
end

end