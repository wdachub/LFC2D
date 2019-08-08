function [flux, addData, flag] = fluxLagrangeSteadySurface2D(C,v,T,nTime,RelTol,region)
% Relate to formula in BA!! 
%
%
%The function implements a calculation of formula (1) of the paper
%
%   Lagrangian Transport through Surfaces in Volume-Preserving Flows
%       Daniel Karrasch (2016)
%
% i.e it calculates the flux due to the velocity field u through the 
% surface C over the time T by calculating the right hand side of formula
% (1)
%
% f is assumed to be constantly 1. To change this the function
% integrateOverSubsets must be modified
%
% For the function to work properly it is required to add the ComplexVisual
% Toolbox of Christian Ludwig to the MATLAB functions
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
%      -------|------------------------------------------------------------
%      T      | Defines the end of the time interval over which the flux
%             | integral shall be calculated. Time interval: [0,T]
%      -------|------------------------------------------------------------
%      nTime  | Defines how many divisions of the time interval are made.
%             | Equals the number of points on the streak line
%      -------|------------------------------------------------------------
%      RelTol | Relative Tolerance used for the ODE solver (wich is ode45)
%      -------|------------------------------------------------------------
%      region | (optional) Only the flux corresponding to this region is
%             | calculated. Particles originated from outside this region
%             | are not considered for the flux. The region is must be
%             | connected without holes and is defined by its bounding 
%             | polygon in the following form:
%             |    region = [x1 y1; 
%             |              x2 y2; 
%             |               ... 
%             |              xn yn]    
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
%     flux    | Scalar value of the calculated flux
%     --------|------------------------------------------------------------
%     addData | Struct containing the following fields:
%             |     Bounding Polygon - bounding polygon of subset D
%             |     DividedPolygons -  table containing the sub loops of
%             |                        the bounding polygon
%             |     WindingNumbers   - Winding numbers corresponding to the 
%             |                        sub loops in 'DividedPolygons'  


%% Test the correctness of the implementation of the velocity field
testVelocityField(v,C);

%% Calculate the bounding polygon of the set D
% TolAngle = 30 * pi/180;
% TolLength = 5e-2;

[BoundingPolygon,BoundingPolygonComponents] = CalcBoundaryD(v,C,T,nTime,RelTol);
% [BoundingPolygon,BoundingPolygonComponents] = CalcBoundaryD_Adaptive(v,C,T,TolAngle,TolLength);

% if exist('region','var')
%     [regionX,regionY] = poly2cw(region(:,1),region(:,2));
%     [BoundingPolygonX,BoundingPolygonY] = poly2cw(BoundingPolygon(:,1),BoundingPolygon(:,2));
%     [IntersectionX, IntersectionY] = polybool('intersection',BoundingPolygonX,BoundingPolygonY,regionX,regionY);
%     BoundingPolygon = [IntersectionX, IntersectionY];
% end

%% Divide Polygon
[BoundingPolygon,simpleLoops,connectedComponents,intersections] = dividePolygon(BoundingPolygon(1:end-1,:));

%% Calculate Winding numbers
[wNo,~,~] = polygonCCwindings(BoundingPolygon,connectedComponents);

%% Filter for NaN's in the winding numbers
NaN_wNo = find(isnan(wNo));
wNo(NaN_wNo) = 0;
flag = false;
if ~isempty(NaN_wNo)
    flag = true;
    fprintf('\n%i winding numbers have been changed from NaN to 0\n\n',numel(NaN_wNo))
end

%% Calculate Flux
if exist('region','var')
    flux = integrateOverSubsets(simpleLoops,wNo, region);
else
    flux = integrateOverSubsets(simpleLoops,wNo);
end

%% Store Results in struct
addData = struct('BoundingPolygon',[real(BoundingPolygon)',imag(BoundingPolygon)'],...
                 'BoundingPolygonComponents',BoundingPolygonComponents,...
                 'Intersections',intersections,...
                 'DividedPolygons',cell2table(simpleLoops),...
                 'WindingNumbers',wNo);


end