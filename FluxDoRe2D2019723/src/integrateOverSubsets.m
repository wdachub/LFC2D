function integral = integrateOverSubsets(simpleLoops,wNo, region)
% Calculates the right hand side of formula (1) on page 3 of the paper. It
% is assumed that f=1 everywhere an thus the integrals are in fact just the
% areas of the polygon loops
%
% List of Parameters:
%
%   IN
%       Name  |         Description
%-------------|------------------------------------------------------------
% dividedPolygons 
%             | Cell that stores the closed Loops in its fields. The closed
%             | loops must be stored as 2D coordinates
%             |        [x1 y1 
%             |         x2 y2 
%             |         ... 
%             |         xn yn]
%     --------|------------------------------------------------------------
%     wNo     | matrix that stores for each closed loop the corresponding
%             | winding number
%     region  | (optional) Only the flux corresponding to this region is
%             | calculated. Particles originated from outside this region
%             | are not considered for the flux. The region is must be
%             | connected and is defined by its bounding polygon in the
%             | following form:
%             |    region = [x1 y1; 
%             |              x2 y2; 
%             |               ... 
%             |              xn yn]
% ------------|------------------------------------------------------------
%                               OUTPUT
%-------------|------------------------------------------------------------
%     integral| Result of the sum over the integrals of formula (1)

integral = 0;

% If a region is specified, intersect each loop with the bonding polygon of
% the region before calculating the area
if exist('region','var') 
    %Convert region to clockwise order
    [regionX,regionY] = poly2cw(region(:,1),region(:,2));
    
    %Loop over all polygon loops
    for i=1:numel(simpleLoops)
        if wNo(i) ~=0 
            %Convert current loop to clockwise order
            [curLoopX,curLoopY] = poly2cw(simpleLoops{i}(:,1),simpleLoops{i}(:,2));
            
            %Intersect current loop and region
            [curPolygonX, curPolygonY] = polybool('intersection',curLoopX,curLoopY,regionX,regionY);
            
            % Intersection can produce non connected loops, separated by 
            % NaN's. 
            %Detect indices of NaNs
            indNaN = find(isnan(curPolygonX));
            indNaN = [0 ; indNaN ; numel(curPolygonX)+1];     %Add 0 and "end"+1 for easy indexing
            
            for j=1:numel(indNaN)-1
                curIndices = indNaN(j)+1 : indNaN(j+1)-1;
                integral = integral + wNo(i)*polyarea(curPolygonX(curIndices) , curPolygonY(curIndices));
            end
        end
    end
    
% If no region is specified, just calculate the area and weight it with the
% winding number
else
    for i=1:numel(simpleLoops)
        if wNo(i) ~= 0
            integral = integral + wNo(i)*polyarea(simpleLoops{i}(:,1),simpleLoops{i}(:,2));
        end
    end

end