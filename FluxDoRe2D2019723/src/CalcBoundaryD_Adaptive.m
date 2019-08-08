function [BoundingPolygon,BoundingPolygonComponents] = CalcBoundaryD_Adaptive(C,v,T,optLagrange)
% Calculates the bounding polygon of the donating regions. The polygon is
% composed of the surface itself, the two streaklines through the start and
% the endpoint of the section and the backward image of the section.
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
%             | the number of points that are stored defines the initial 
%             | spatial resolution. Each point is advected back and the
%             | back-advected points define the first candidate for the 
%             | backadvected surface which is then checked for the need of
%             | refinement.
%       ------|------------------------------------------------------------
%       v     | Velocity field. Either a anonymous function or the handle 
%             | to a m-Functionin the current folder. It must be defined 
%             | for evaluation of multiple Position vectors at the given 
%             | time. The following form must be satisfied
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
%       T     | Defines the time interval over which the flux integral
%             | shall be calculated. Time interval: T = [t1,t2]
%       ------|------------------------------------------------------------
% optLagrange | Struct that stores the options for the Lagrangian Method.
%             | Initialized by setOptLagrange() before calling
%             | fluxLagrangeSteadySurface2D_adaptive(...)
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%       ------|------------------------------------------------------------
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
%       ------|------------------------------------------------------------
% BoundingPolygonComponents
%             | 
%             | Struct that stores the parts of the bounding polygon
%             | including the refined surface, the streaklines and the
%             | corresponding time discretization and the backward-image of 
%             | the surface 


%% Calculate the streak lines through start and endpoint of C
[StreakLine1,StreakLine2,tau1,tau2] = CalcStreaklinesAdaptive(C(1,:),C(end,:),T,v,optLagrange);


%% Perform backadvection of C
[C_backadvec,C_refined] = BackadvectSectionAdaptive(C,v,T,optLagrange);

%% Assemble polygon
BoundingPolygon = [C_refined(1:end-1,:);
                     StreakLine2;
                     flipud(C_backadvec(2:end-1,:));
                     flipud(StreakLine1)];

% Store Components of the polygon
BoundingPolygonComponents = struct('RefinedSurface',C_refined,...
                                   'StreakLine1',StreakLine1,...
                                   'StreakLine2',StreakLine2,...
                                   'DiscrT_Streakline1',tau1,...
                                   'DiscrT_Streakline2',tau2,...
                                   'BWimageSurface',C_backadvec);
                                   
                 
end