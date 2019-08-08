% This file demonstrates the use of the Code for the adaptive Lagrangian 
% Method for flux integration on the example of the cross-channel mixer
% described in the paper. 
% In the setup of the cross-channel mixer the flux across a defined section
% is integrated (possibly restricted to a region of interest)

%% Add required Paths
if ismac
    separator = '/';
elseif isunix
    separator = '/';
elseif ispc
    separator = '\';
else
    disp('Platform not supported')
end

addpath(genpath('..'))
% addpath(['..',separator,'SubFunctions'])
% addpath(['..',separator,'FunctionsVisualization'])
% addpath(['..',separator,'..',separator,'ComplexVisual'])
% addpath(['..',separator,'..',separator,'CrossChannelMixer'])
% addpath(['..',separator,'..',separator,'CrossChannelMixer',separator,'Functions'])

%% Step 1: Define Setting
% Cross-Channel mixer:
% The parameters of the cross-channel mixer are stored in a struct calles
% CrossChannelMixer. This includes the geometric properties as well as the
% Information to crate the velocity field. To modify the parameters open 
% the createCrossChannelMixer-file
createCrossChannelMixer

% Section:
% The section is defined as a polyline. Possibly the first definition
% stated here is refined during the calculation. Defined as (nPoints x 2)
% matrix, with the x Positions in the first column and the y Positions in
% the second column
% Section must oriented such that the normal vector field from the Eulerian
% flux integral points to the RIGHT.

nPointsC = 20;
% horizontal section
C = [linspace(CrossChannelMixer.p(end)+CrossChannelMixer.r(end) , CrossChannelMixer.p(1)-CrossChannelMixer.r(1) , nPointsC)' , ...
     zeros(nPointsC,1)];
% vertical section
% C = [(CrossChannelMixer.p(end)+CrossChannelMixer.r(end))*ones(nPointsC,1),...
%     linspace(0,CrossChannelMixer.f,nPointsC)'];
 
% Time Interval:
% The flux is integrated over a time interval T=[t0,t1]
T = [0,5.5];

% Velocity Field:
% The velocity field is created from the parameters stored in the
% CrossChannelMixer struct
v = @(t,Pos)VelocityFieldCrossChannelMixer(t,Pos,CrossChannelMixer);



%% Step 2: Define Region (optional)
% The flux calculation can be restricted to a material region of interest.
% Only paricles originating from within this region contribute to the flux
% integration. This region is defined by its closed bounding polygon. This 
% region is an optional parameter. If no region shall be used no variable 
% 'PolygonRegion' must exist in the workspace

% Example (Rectangular lower half)
yMax = 0;
yMin = -CrossChannelMixer.e;
xMax = CrossChannelMixer.p(end) + CrossChannelMixer.r(end);
maxVel = (CrossChannelMixer.e + CrossChannelMixer.f)^2 / 4 / CrossChannelMixer.e / CrossChannelMixer.f * CrossChannelMixer.U;   %Formula for velocity at vertex
xMin = CrossChannelMixer.p(1) - CrossChannelMixer.r(1) - T(end)*maxVel;

PolygonRegion = [xMin , yMin;
                 xMax , yMin;
                 xMax , yMax;
                 xMin , yMax;
                 xMin , yMin];
             
             
             
%% Step 3: Options for the Lagrangian Method
% Several options are available for the Lagrangian Method. Open
% setOptLagrange.m for more information and modification
setOptLagrange



%% Step 4: Flux integration using the Lagrangian Method
% Depending on whether a bounding Polygon for the region of interest is
% defined or not the Lagrangian method is executed as follows. 'IntFlux' is
% the scalar result of the flux integration, 'addData' is a struct
% containing additional Data
tic
if exist('PolygonRegion','var')
    [intFlux, addData,flag] = fluxLagrangeSteadySurface2D_adaptive(C,v,T,optLagrange,PolygonRegion);
else
    [intFlux, addData,flag] = fluxLagrangeSteadySurface2D_adaptive(C,v,T,optLagrange);
end
toc
disp(['Integrated (restricted) flux = ' num2str(intFlux)])

%% Step 5: Visualizing results

% Using the struct 'addData' the following visualizations can be performed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the donating regions colored depending on the corresponding 
% winding number. If existent plot of Region of interest
[hFigure1,ax1] = plotCrossChannelMixer(CrossChannelMixer,xMin-1,xMax+1);
if exist('PolygonRegion','var')
    plotDonatingRegions(addData,C,'RegionOfInterest',PolygonRegion,'hFigure',hFigure1);
else
    plotDonatingRegions(addData,C,'hFigure',hFigure1)
end
% axis([-5. 6. -1.85 1.55])
% axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the components of the bounding polygon of the donating region
[hFigure2,ax2] = plotCrossChannelMixer(CrossChannelMixer,xMin-1,xMax+1);
plotComponentsD(addData,hFigure2)
% axis([-5. 6. -1.85 1.55])