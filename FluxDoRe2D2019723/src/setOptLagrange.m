% Set options used for the adaptive Lagrangian flux integration and stores
% them in the struct optLagrange


%% Options ODE solver
% The standard MATLAB ode45 solver is used. Here absolute and relative
% tolerances can be given as .ODE_AbsTol and .ODE_RelTol. Optional
% parameter
optLagrange.ODE_RelTol = 1e-6;
optLagrange.ODE_AbsTol = 1e-8;


%% Adaptive refinement
% Both the streaklines and the backadvected sections are adaptively refined
% to meet some prescribed tolerances. The procedure is as follows:
% An initial discretization of the time interval in case of the streaklines
% and the discretization of the section given by the user are advected
% back, resulting in polylines that are candidates for the streaklines and 
% the backadvected section. Angles between adjacent segments and the 
% lengths of the segments are measured. If an angle exceeds a prescribed 
% tolerance both adjacent segments need refinement. If a length exceeds the
% prescribed tolerance, this segement needs refinement. A fixed number of 
% new points is then inserted in the corresponding time segment or 
% segment of the section, the new points are advected back and 
% the resulting refined polylines are again checked for the need of
% refinement.

% DIFFERENT NUMBER OF POINTS/ANGLES ETC. FOR STREAKLINE AND SECTION??


% Number of new Points that are inserted in one segment per refinement
optLagrange.nNewPoints = 2;

% Maximum number of refinements. If reached, the last result is accepted
optLagrange.MaxNumberRefinement = 20;

% Number of intervals for the initial discretization of the time interval.
% Linearly spaced in [t1 t2]
optLagrange.nIniIntervalT = 100;

% Tolerance for the angle. Angles that exceed this value are to be refined.
% Angle given in rad
optLagrange.TolAngle = 20 * pi/180;

% Tolerance for the length. Lengths that exceed this value are to be
% refined
optLagrange.TolLength = 1e-1;

% Factor of TolLength below wich edges are not refined. If an edge is to be
% refined due to angular tolerance but its length is shorter than 
% factorStopRefinement*TolLength then it is not refined. (Short edges
% contribute to the error very little)
% Must be from [0,1]
optLagrange.factorStopRefinement = 0.2;