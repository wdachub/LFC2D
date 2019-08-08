function [tSteps,StreakLines] = calcStreakLinesODE45(v,x0,T,nTime,RelTol)
% Calculates the streaklines trough the points stored in x0 over the 
% time interval [0,T]. The time interval is divided in nTime steps and 
% RelTol is used as tolerance for the ode45 routine
%
% List of Parameters:
%
% =========================================================================
%   INPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
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
%       x0    | Points through which the strak lines are calculated.
%             | The following data structure must be satisfied
%             |
%             |                [x1
%             |                 x2
%             |                 ...
%             |         x0 =    xn
%             |                 y1
%             |                 y2
%             |                 ...
%             |                 yn]
%       ------|------------------------------------------------------------
%       T     | Defines the end of the time interval for which the streak 
%             | lines shall be calculated
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
%       tSteps| Time parameters for the streak lines. (see below for 
%             | connection to streak lines
%             |     tSteps = [0,t1, t2, ... T]
%       ------|------------------------------------------------------------
%  StreakLines| Stores the points of the streak lines in the following form
%             |
%             |                [x1(0)   x1(t1)  ... x1(T)] 
%             |                |x2(0)   x2(t1)  ... x2(T)|
%             |                |...                 ...  |
%             |StreakLines =   |xn(0)   xn(t1)  ... xn(T)|
%             |                |y1(0)   y1(t1)  ... y1(T)|
%             |                |y2(0)   y2(t1)  ... y2(T)| 
%             |                |...                 ...  |
%             |                [yn(0)   yn(t1)  ... yn(T)]
%             |  
%             | For example the points of the first streak line are 
%             | stored in the first row (x-coordinates) and in the n+1st 
%             | row (y-coordinates) 
%             | The connection between the points of a streak line and the
%             | time steps is  the following: If a point is advected
%             | forward as long as the corresponding time parameter 
%             | indicates, the advection ends at the initial position in x0
% -------------------------------------------------------------------------

% Define time steps
tSteps = linspace(T(1),T(2),nTime);

% Initialize StreakLines-matrix
StreakLines = zeros(size(x0,1),nTime);

%Store the initial step
StreakLines(:,1) = x0;

% Set options for the ODE solver
opt = odeset('RelTol',RelTol,'AbsTol',1e-6);

% Loop over the time Steps
for i=2:nTime
    %Advect x0 back. Time defined by current value of tSteps
    [~, Positions] = ode45(v,[tSteps(i),mean([tSteps(i),0]),0],x0,opt);
    
    %Store last Positions as points on the streak lines 
    StreakLines(:,i) = Positions(end,:)'; 
end

end
