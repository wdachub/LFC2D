function [BoundingPolygon,BoundingPolygonComponents] = CalcBoundaryD(v,C,T,nTime,RelTol)
% Calculates the bounding polygon of the donating regions. The polygon is
% composed of the surface itself, the two streaklines through the start and
% the endpoint and the backward image of the section.
%
% List of Parameters:
%
% =========================================================================
%   IN
%       Name  |         Description
%-------------|------------------------------------------------------------
%       C     | Defines the surface, stored as 
%             |        [x1 y1 
%             |         x2 y2 
%             |         ... 
%             |         xn yn]
%             | the number of points that are stored defines the spatial
%             | resolution (Each point is advected back and the
%             | back-advected points define the backadvected surface.
%       ------|------------------------------------------------------------
%       v     | Either a anonymous function or the handle to a m-Function
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
%       ------|------------------------------------------------------------
%       T     | Defines the end of the time interval over which the flux
%             | integral shall be calculated. Time interval: [0,T]
%       ------|------------------------------------------------------------
%       nTime | Defines how many divisions of the time interval are made.
%             | Equals the number of points on the streak line
%       ------|------------------------------------------------------------
%       RelTol| Relative Tolerance used for the ODE solver (wich is ode45)
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
% BoundingPolygon 
%             | Stores the bounding Polygon of the subset D as
%             |        [x1 y1 
%             |         x2 y2 
%             |         ... 
%             |         xn yn]
%             |
%             | The order of the parts is
%             |        [Surface
%             |         Streak Line starting from the endpoint C
%             |         backward image of the section
%             |         Streak line ending at the startpoint of C]

nOfPointsC = size(C,1);

%% Calculate the streak lines through start and endpoint of C
[~,StreakLines] = calcStreakLinesODE45(v,C([1;nOfPointsC;nOfPointsC+1;end]),T,nTime,RelTol);
StreakLine1 = [StreakLines(1,:)',StreakLines(3,:)'];
StreakLine2 = [StreakLines(2,:)',StreakLines(4,:)'];

%% Perform backadvection of C
%Only the last position is needed => use StraklineFunction with 2 time
%steps
[~,SectionBackwards] = calcStreakLinesODE45(v,C(:),T,2,RelTol); 
BWimageSection = [SectionBackwards(1:nOfPointsC,2),SectionBackwards(nOfPointsC+1:end,2)];

%% Assemble polygon
BoundingPolygon = [C(1:end-1,:);
                     StreakLine2;
                     flipud(BWimageSection(2:end-1,:));
                     flipud(StreakLine1)];

% Store components of the polygon
BoundingPolygonComponents = struct('Surface',C,...
                                   'StreakLine1',StreakLine1,...
                                   'StreakLine2',StreakLine2,...
                                   'BWimageSurface',BWimageSection);
                                   
                 
end