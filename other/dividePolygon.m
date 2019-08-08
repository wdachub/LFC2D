function [Polygon,simpleLoops,connectedComponents,intersections] = dividePolygon(Polygon)
% Divides 'Polygon' in simple, closed loops. The coordinates of the polygon
% are re-interpreted as points in the complex plane. That way the
% ComplexVisual Toolbox by Christian Ludwig can be used
%
% List of Parameters:
%
% =========================================================================
%   INPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
%     Polygon | Stores the bounding Polygon of the donating region D as
%             |        [x1 y1 
%             |         x2 y2 
%             |         ... 
%             |         xn yn]
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
%-------------|------------------------------------------------------------
%     Polygon | Stores the initial Polygon as a complex polygon augmented 
%             | by the crossing points. Each crossing point is stored twice 
%             | i.e. both segments that are intersecting are augmented by 
%             | the coordinates of the crossing point. The Polygon is
%             | stored as a complex polygon because this is required to use
%             | the ComplexVisual Toolbox to calculate the winding numbers.
%     --------|------------------------------------------------------------
% simpleLoops | Cell that stores the closed Loops in its fields. The closed
%             | loops are stored as coordinates in the x-y-plane
%             |        [x1 y1 
%             |         x2 y2 
%             |         ... 
%             |         xn yn]
%     --------|------------------------------------------------------------
% connectedComponents
%             | Stores the connected Components required for the
%             | calculation of the winding number. Stored in a cell
% -------------------------------------------------------------------------
%
% The following functions are provided by the ComplexVisual Toolbox by 
% Christian Ludwig
%
%   intersectPolylines
%   augmentPolylinePoI
%   calcConnectedComps
%   getPolygonCCpaths
%
% Lizenz:
%
% The ComplexVisual Matlab-Toolbox
% Copyright (c) 2011, C. Ludwig <ludwig@ma.tum.de>
% All rights reserved.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


%% Find intersections and add Crossing Points to the Polygon
% ALTERNATIVE 1: should be used when the curve doesn't have too many points
% Polygon = Polygon';
% intersections = InterX(Polygon);
% % complexification for later use
% Polygon(2,:) = Polygon(2,:)*1i;
% Polygon = sum(Polygon,1);
% % END OF ALTERNATIVE 1

% % ALTERNATIVE 2: is faster on curve with many points (>5000), but may
% miss intersection points, fixable by increasing variable "block" in line
% 242 of intersectPolylines.m
Polygon = Polygon';
Polygon(2,:) = Polygon(2,:)*1i;
Polygon = sum(Polygon,1);
intersections = intersectPolylines(Polygon,[]);
% % END OF ALTERNATIVE 2

[Polygon,~,intersections] = augmentPolylinePoI(Polygon,[],intersections);

%% Calculate Connected Components and Divide polygon 
connectedComponents = calcConnectedComps(intersections);

simpleLoops = getPolygonCCpaths(Polygon,connectedComponents);

%% Re-Format Loops as coordinates in the x-y-plane
for i=1:numel(simpleLoops)
    simpleLoops{i} = [real(simpleLoops{i})',imag(simpleLoops{i})'];
end

end