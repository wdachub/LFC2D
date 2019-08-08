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

%% Step 0: Eulerian flux integral for comparison
syms x y t real

t_init = 0.0;
t_end = 2.2;
x1 = -0.1;
x2 = 0.3;

rho = exp(-2*t - exp(-2*t)*(x^2 + y^2));
vsym = [x + 2*pi*y; -2*pi*x + y];
 
% Check conservation law
consLaw = simplify(diff(rho, t) + divergence(rho*vsym, [x;y]))
 
% Evaluate Euler Integral
% Evaluated for y = 0 and n = [0;-1], so we build a simplified integrand
integrand = subs(rho, y, 0) * subs(vsym, y, 0)' * [0; -1];

% Space integration
IntSpace = simplify( int(integrand, x, x1, x2) );
 
% Time integration
IntComplete = simplify( int(IntSpace, t, t_init, t_end) );

intEuler = double( IntComplete )
intEulerNum = integral2(matlabFunction(integrand), t_init, t_end, x1, x2, 'AbsTol', 0, 'RelTol', 1e-12)
abs(intEuler-intEulerNum) / abs(intEuler)
%% Step 1: Define setting
% Section:
% The section is defined as a polyline. Possibly the first definition
% stated here is refined during the calculation. Defined as (nPoints x 2)
% matrix, with the x Positions in the first column and the y Positions in
% the second column
% Section must oriented such that the normal vector field from the Eulerian
% flux integral points to the RIGHT.

nPointsC = 50;
C = [linspace(x1, x2, nPointsC)', zeros(nPointsC, 1)];

% Time interval:
% The flux is integrated over a time interval T=[t0,t1]
T = [t_init, t_end];

% Velocity field:
% The velocity field must be a function that allows the evaluation of
% multiple positions for one instant of time. For more information and
% modification of the velocity field open velocityFieldComp.m 
v = @(t, x) velocityFieldComp(t, x);

%% Step 3: Options for the Lagrangian Method
% Several options are available for the Lagrangian Method. Open
% setOptLagrange.m for more information and modification

hs = [0.05, 0.01, 0.005, 0.001];
% optionally: add point grid, TODO: replace by bounding box coords, adjust resolution
refine = false;
if refine
    [X, Y] = meshgrid(-1:optLagrange.TolLength:1, -1:optLagrange.TolLength:1);
    Pi = [X(:), Y(:)];
end
% optionally: visualize triangulation of each donating region
visualize = false;

setOptLagrange

% test different length tolerances
for j=1:length(hs)
    optLagrange.TolLength = hs(j);

    %% Step 4: Flux integration using the Lagrangian Method
    % Depending on whether a bounding Polygon for the region of interest is
    % defined or not the Lagrangian method is executed as follows. 'IntFlux' is
    % the scalar result of the flux integration, 'addData' is a struct
    % containing additional Data

    % POOR MAN'S INTEGRATION
    % manual integration, to be included in integrateOverSubsets.m
    density = matlabFunction(subs(rho, t, t_init));

    tic
    [intFlux, addData, flag] = fluxLagrangeSteadySurface2D_adaptive(C, v, T, optLagrange);
    % [intFlux, addData,flag] = fluxLagrangeSteadySurface2D(C, v, T, N, optLagrange.ODE_RelTol);

    DP = table2cell(addData.DividedPolygons);
    wNo = addData.WindingNumbers;
    [lam,w] = quadpts(4);
    nq = size(lam, 1);
    intLagrange = 0;
    for i=1:numel(addData.DividedPolygons) % for each simple loop with nonzero wNo
        if not(wNo(i)==0)
            DPi = DP{1,i};
            index = (sum(DPi(1:end-1,:) == DPi(2:end,:),2)<2); % exclude duplicates
            q = sum(index);
            DPi = DPi(index,:);
            if refine
                DT = delaunayTriangulation([DPi; Pi], [(1:(q-1))', (2:q)'; q 1]);
            else
                DT = delaunayTriangulation(DPi, [(1:(q-1))', (2:q)'; q 1]);
            end
            inside = isInterior(DT);
            ti = DT.ConnectivityList(inside, :);
            p = DT.Points;
            if visualize
                plotDonatingRegions(addData,C)
                triplot(DT.ConnectivityList(inside, :),DT.Points(:,1),DT.Points(:,2))
            end
            v1 = p(ti(:,3),:) - p(ti(:,2),:); v2 = p(ti(:,1),:) - p(ti(:,3),:); v3 = p(ti(:,2),:) - p(ti(:,1),:);
            area = 0.5*(-v3(:,1).*v2(:,2) + v3(:,2).*v2(:,1));       % areas of triangles
            intList = zeros(size(ti,1),1);
            for k = 1:nq
                x = lam(k,1)*p(ti(:,1),:) + lam(k,2)*p(ti(:,2),:) + lam(k,3)*p(ti(:,3),:);
                intList = intList + w(k)*density(x(:,1), x(:,2));
            end
            intLagrange = intLagrange + wNo(i) * sum(area .* intList);
        end
    end

    disp(['Length tolerance = ' num2str(hs(j)) '.'])
    disp(['Integration took ' num2str(toc) ' seconds.'])
    disp(['Integrated compressible (restricted) flux = ' num2str(intLagrange)])
    disp(['Relative error: ' num2str(abs(intLagrange-intEuler)/abs(intEuler))])
end
%% Step 5: Visualizing results
% Using the struct 'addData' the following visualizations can be performed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the donating regions colored depending on the corresponding 
% winding number. If existent plot of Region of interest
if exist('PolygonRegion', 'var')
    plotDonatingRegions(addData, C, 'RegionOfInterest', PolygonRegion);
else
    plotDonatingRegions(addData, C)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the components of the bounding polygon of the donating region
plotComponentsD(addData)
