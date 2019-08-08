function [streakline1,streakline2,tau1,tau2] = CalcStreaklinesAdaptive(x1,x2,T,v,optLagrange)
% Calculates the straklines through x1 and x2 for time interval T = [t1,t2]
% The discretisation of the time interval is adaptively refined to reach
% the maximal lengths and angles prescribed by the settings in optLagrange
%
% List of Parameters:
%
% =========================================================================
%   IN
%       Name  |         Description
%-------------|------------------------------------------------------------
%       x1/x2 | Point through which the streaklines shall be calculated.
%             | Must be given as a row vector
%       ------|------------------------------------------------------------
%       T     | Defines the time interval for which the streaklines shall
%             | be calculated. Time interval: T = [t1,t2]
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
% optLagrange | Struct that stores the options for the Lagrangian Method.
%             | Initialized by setOptLagrange() before calling
%             | fluxLagrangeSteadySurface2D_adaptive(...)
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
% streakline1 | Space points of the streaklines. Given as
%      /      |                      [x1 y1
% streakline2 |     streakline1/2 =    ...
%             |                       xn yn]
%      -------|------------------------------------------------------------
%     tau1/2  | Refined time points corresponding to the space points of
%             | the streaklines

% Set options for ODE solver
optODE = odeset;
if isfield(optLagrange,'ODE_AbsTol')
    optODE = odeset(optODE,'AbsTol',optLagrange.ODE_AbsTol);
end
if isfield(optLagrange,'ODE_RelTol')
    optODE = odeset(optODE,'RelTol',optLagrange.ODE_RelTol);
end

% Calculate Streaklines
[streakline1, tau1] = StreaklineSinglePoint(x1);
[streakline2, tau2] = StreaklineSinglePoint(x2);




%% Nested function for the calculation of the streaklines
function [streakline, tau] = StreaklineSinglePoint(x)
% Calculates the streakline trough x in the time interval T = [t2 t2]. Uses
% adaptivity. 
% Nested function in CalcStreaklines

% Initial discretization of the time interval
tau = linspace(T(1),T(2),ceil(optLagrange.nIniIntervalT)+1)';

% Initial backadvection
streakline = backadvectionPoint(x,tau,T(1));

% Find indices of edges that are to be refined
indEdgesToRefine = findEdgesToRefine(streakline,optLagrange);

counter = 0;

while (~isempty(indEdgesToRefine) && counter<optLagrange.MaxNumberRefinement)
    % Create new discretization points for the time intervals corresponding
    % to the edges that need to be refined (+2 as first and last are
    % removed)
    newPointsTime = linspaceV(tau(indEdgesToRefine),tau(indEdgesToRefine+1),optLagrange.nNewPoints)';
    
    % Calculate corresponding new points on the streakline
    newPoints_StreakLine = backadvectionPoint(x,newPointsTime(:),T(1));
    
    % Add new taus and new points to result 
    tau = [tau;newPointsTime(:)];
    streakline = [streakline;newPoints_StreakLine];
    
    % Sort new points according to time stamps
    [tau,I] = sort(tau);
    streakline = streakline(I,:);
    
    % Find indices of edges that are to be refined
    indEdgesToRefine = findEdgesToRefine(streakline,optLagrange);
    counter = counter+1;

end
fprintf('Refinement of streakline through (%.2f,%.2f) required %i steps.\n',x(1),x(2),counter)
end




%% Nested function for the 
function Pos = backadvectionPoint(x,tau,t0)
% Advects the point x = [x,y] from the time points stored in the vector tau 
% back to time t0. Results are stored in Pos.
% Nested function in CalcStreaklines

% Pre-allocate storage
Pos = zeros(numel(tau),2);



% Loop over all time points in tau. Use parallelization
parfor i=1:numel(tau)
    if tau(i) == t0     
        Pos(i,:) = x';  
        
    else    %Advect x0 back. Time defined by current value of tau 
        [~, Positions] = ode45(v,[tau(i),mean([tau(i),t0]),t0],x,optODE);   
        % Only the last postion is needed. The [tau(i),mean(),t0] 
        % construction reduces the size of the output vector
        
        %Store last Positions as points on the streak lines 
        Pos(i,:) = Positions(end,:); 
    end
    
end

end
end




%% Function to determine the indices of segemnts that need refinement
function indEdgesToRefine = findEdgesToRefine(Points,optLagrange)
% Measures the angles between the segments of the poly-line defined by the
% points. If the absolute value of the angle is bigger than the TolAngle,
% both adjacent edges need refinement. Also the lengths of the segments are
% calculated and sements with a length that exceeds the tolerance must also
% be refinded. The function returns the indices of  those lines that need 
% refinement. 

% Calculate connecting vectors
v = diff(Points);

% Calculate scalar products of adjacent vectors
sp = sum(v(2:end,:).*(v(1:end-1,:)),2);

% Calculate norms of vectors
norms = sqrt(sum(v.^2,2));

% Calculate Angles
angles = acos(sp./norms(1:end-1)./norms(2:end));

% Find indices of angles that are too big
indAngTooBig = find(abs(angles) > optLagrange.TolAngle & norms(1:end-1) > optLagrange.factorStopRefinement*optLagrange.TolLength...
    & norms(2:end) > optLagrange.factorStopRefinement*optLagrange.TolLength);

% Find indices of segments that are too long
indVecTooLong = find(norms>optLagrange.TolLength);

% Return indices of the edges adjacent to the angles that are too big
indEdgesToRefine = unique([indAngTooBig;indAngTooBig+1;indVecTooLong]);

end




%% Function for linearized point insertion 
function linSpaces = linspaceV(u,v,n)
% Creates n equidistant points between each pair u(i) v(i), respectively. u
% and v must be column vectors of the same length. The result is a matrix 
% containing the new points in each row. 
% 
% Example: u = [1;2;3] , v = [4;5;7] , n = 4
%
% linSpaces = [ 1.6   2.2   2.8   3.4
%               2.6   3.2   3.8   4.4
%               3.8   4.6   5.4   6.2]

linSpaces = u*ones(1,n+2) + (v-u)*linspace(0,1,n+2);

linSpaces = linSpaces(:,2:end-1);
end