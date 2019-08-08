function [intFlux, addData, flag] = fluxLagrangeSteadySurface2D_adaptive(C,v,T,optLagrange,region)
% The function implements an adaptive calculation of formula (1) of the 
% paper
%
%   Lagrangian Transport through Surfaces in Volume-Preserving Flows
%       Daniel Karrasch (2016)
%
% i.e it integrates the flux due to the velocity field v through the 
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
%      -------|------------------------------------------------------------
%      T      | Defines the time interval over which the flux integral
%             | shall be calculated. Time interval: T = [t1,t2]
%      -------|------------------------------------------------------------
% optLagrange | Struct that stores the options for the Lagrangian Method.
%             | Initialized by setOptLagrange() before calling
%             | fluxLagrangeSteadySurface2D_adaptive(...)
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
%     intFlux | Scalar value of the integrated flux
%     --------|------------------------------------------------------------
%     addData | Struct containing the following fields:
%             |     Bounding Polygon - bounding polygon of subset D
%             |     DividedPolygons -  table containing the sub loops of
%             |                        the bounding polygon
%             |     WindingNumbers   - Winding numbers corresponding to the 
%             |                        sub loops in 'DividedPolygons'
%     --------|------------------------------------------------------------
%     flag    | returns true if some of the winding numbers could not have
%             | been computed unambiguously


%% Test the correctness of the implementation of the velocity field
testVelocityField(v,C);

%% Calculate the bounding polygon of the set D
[BoundingPolygon,BoundingPolygonComponents] = CalcBoundaryD_Adaptive(C,v,T,optLagrange);

%% Divide Polygon
[BoundingPolygon,simpleLoops,connectedComponents,intersections] = dividePolygon(BoundingPolygon);

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
% the distribution must be evaluated at the time T(1)
if exist('region','var')
    intFlux = integrateOverSubsets(simpleLoops,wNo, region);
else
    intFlux = integrateOverSubsets(simpleLoops,wNo);
end

%% Store Results in struct
addData = struct('BoundingPolygon',[real(BoundingPolygon)',imag(BoundingPolygon)'],...
                 'BoundingPolygonComponents',BoundingPolygonComponents,...
                 'Intersections',intersections,...
                 'DividedPolygons',cell2table(simpleLoops),...
                 'WindingNumbers',wNo);


end