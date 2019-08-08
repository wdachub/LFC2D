% This file demonstrates the use of the Code for the adaptive Lagrangian 
% Method for flux integration. For a fixed 2D section and a time-
% dependent velocity field the transport of a conserved quantity across the 
% section is calculated over a given time interval. For the quantity, the initial 
% distribution must be known. In this code the distribution is assumed to 
% be constantly 1.

%% Add required paths
if ismac, separator = '/';
elseif isunix, separator = '/';
elseif ispc, separator = '\';
else, disp('Platform not supported')
end

addpath(genpath('..'))


%% Step 1: Define setting
% Section:
% The section is defined as a polyline. Possibly the first definition
% stated here is refined during the calculation. Defined as (nPoints x 2)
% matrix, with the x Positions in the first column and the y Positions in
% the second column
% Section must oriented such that the normal vector field from the Eulerian
% flux integral points to the RIGHT.
nPointsC = 300;
C = [0.75*ones(nPointsC,1), linspace(-0.2,1.2,nPointsC)'];

% Time interval:
% The flux is integrated over a time interval T=[t0,t1]
T = [0,2.5];

% Velocity field:
% The velocity field must be a function that allows the evaluation of
% multiple positions for one instant of time. For more information and
% modification of the velocity field open velocityField.m 
v = @(t,x) velocityFieldKarrasch(t,x);



%% Step 2: Define Region (optional)
% The flux calculation can be restricted to a material region of interest.
% Only paricles originating from within this region contribute to the flux
% integration. This region is defined by its bounding polygon. This region 
% is an optional parameter. If no region shall be used no variable 
% 'PolygonRegion' must exist in the workspace

% Example (nPointsRegion-sided-figure)
nPointsRegion = 99;
angRegion = linspace(0,1,nPointsRegion+1)'*2*pi;

PolygonRegion = [0.3+.4*cos(angRegion), .5+.6*sin(angRegion)];



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
if exist('PolygonRegion','var')
    plotDonatingRegions(addData,C,'RegionOfInterest',PolygonRegion);
else
    plotDonatingRegions(addData,C)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the components of the bounding polygon of the donating region
plotComponentsD(addData)