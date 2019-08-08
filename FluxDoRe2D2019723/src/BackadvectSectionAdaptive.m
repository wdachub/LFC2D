function [C_backadvec,C_refined] = BackadvectSectionAdaptive(C,v,T,optLagrange)
% Calculates the backadvected surface. Backadvection from T(2) to T(1). The
% discretization of the surface C is adaptively refined to reach the
% maximal lengths and angles prescribed by the settings in optLagrange.
% Refinement of C by linear interpolation.
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
%       T     | Defines the time interval for which the section shall be
%             | advected back. Time interval: T = [t1,t2]
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
% C_backadvec | Polyline of the backadvected section
%      -------|------------------------------------------------------------
%   C_refined | Refined section C

% Set options for ODE solver
optODE = odeset;
if isfield(optLagrange,'ODE_AbsTol')
    optODE = odeset(optODE,'AbsTol',optLagrange.ODE_AbsTol);
end
if isfield(optLagrange,'ODE_RelTol')
    optODE = odeset(optODE,'RelTol',optLagrange.ODE_RelTol);
end 


%First backadvection with given discretization of C
[~, Points] = ode45(v,[T(2),mean(T),T(1)],C(:),optODE);
% Only the last postion is needed. The [tau(i),mean(),t0] 
% construction reduces the size of the output vector
C_backadvec = reshape(Points(end,:),size(Points,2)/2,2);

% Calculate Angles and identify edges that need refinement
indEdgesToRefine = findEdgesToRefine(C_backadvec,optLagrange);

counter = 0;

while (~isempty(indEdgesToRefine) && counter<optLagrange.MaxNumberRefinement)
    % Create new discretization points for C corresponding
    % to the edges that need to be refined (+2 as first and last are
    % removed)
    newPointsCx = linspaceV(C(indEdgesToRefine,1),C(indEdgesToRefine+1,1),optLagrange.nNewPoints)';
    newPointsCy = linspaceV(C(indEdgesToRefine,2),C(indEdgesToRefine+1,2),optLagrange.nNewPoints)';
    
    % Advect new points back
    [~,newPoints] = ode45(v,[T(2),mean(T),T(1)],[newPointsCx(:);newPointsCy(:)],optODE);
    
    % Create auxiliary entries to index-vector for sorting
    auxInd = linspaceV(indEdgesToRefine,indEdgesToRefine+1,optLagrange.nNewPoints)';
    auxInd = auxInd(:);
    
    % Sort auxiliary entries to index vector to create new index vector
    [~,I] = sort([(1:size(C,1))';auxInd]);
    
    % Add new Points to C_backadvec and C
    C = [C;[newPointsCx(:),newPointsCy(:)]];
    C_backadvec = [C_backadvec;reshape(newPoints(end,:),size(newPoints,2)/2,2)];
    
    % Sort C and C_backadvec according to I
    C = C(I,:); 
    C_backadvec = C_backadvec(I,:);
    
    % Calculate Angles and identify edges that need refinement
    indEdgesToRefine = findEdgesToRefine(C_backadvec,optLagrange);
    counter = counter+1;
    
end

C_refined = C;
fprintf('\nBackadvected section refinement required %i steps.\n',counter)
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
