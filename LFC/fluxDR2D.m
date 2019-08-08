function [flux] = fluxDR2D...
    (curve, vel, pfun,t0,te ,nSeg,order,option)


% =========================================================================
% INPUTS
% =========================================================================
%       Name  |         Description
%-------------|------------------------------------------------------------
% curve       | data format:
%             | 1  function handle
%             | 2  vector which record the coordinate of points on curve.
%             | Notice:
%             | 1 the orientation of curve is very important which
%             | will affect the sign of result.  n is normal vector, v is the
%             | tangent vector of curve. They form the local coordinate system
%             | (n,v).
%             | 2 Close curve has faster algorithm. Equal first and last
%             | component of vector to indicate the curve is closed.
%             |
%             |
% ------------|------------------------------------------------------------
%   vel       |  the velocity field.
%             |  data format£º
%             |  1  function handle
%             |  2  cell, cell{i} is matrix. (vel{1},vel{2}) record the
%             |  coordinates of nodes,(vel{3},vel{4}) record the velocity
%             |  on each nodes.
% ------------|------------------------------------------------------------
%   order     |  the intended order-of-convergence, which is the minimum of
%             |  the order of spline and runge kutta
%             |  data format£º
%             |  1. [order] is a number, that means four order are same.
%             |  2. [curve_spline_order, DR_spline_order,RKorder£¬cubatureOrder]
% ------------|------------------------------------------------------------
%   nSeg      | the number of segments to divide the curve, so there is nSeg+1
%             | nodes on the curve.
%             | In order to implement spline interpolation, nSeg must be
%             | at least order-1;
%             |
% ------------|------------------------------------------------------------
%   pfun      | the passive function to be integrated on the donating
%             | region.
%             |
% ------------|------------------------------------------------------------
%   te, t0    | ending and beginning times of the test.
%             | t0 will be zero if missing.
%             |
% ------------|------------------------------------------------------------
%             | option(1) If we use spline approximating given curve and sampling the spline when construct DR,
%             | adding input point into sampling point set will increase the accuracy.
%             | option(1)=1 adding input point;
%             | option(1)=0 not adding input point;
%             | option(2) refinement
%             | option(2)=1 the max tolerance distance in refinement is set to the
%             | average of each segment.
%  option     | option(2)=0, don't do refinement
%             | default setting: not adding, not refine.
% =========================================================================
% OUTPUTS
% =========================================================================
%       Name  |         Description
% ------------|------------------------------------------------------------
%   flux      | flux pass through the given during [t0, te].
%             |
%             |
% ------------|------------------------------------------------------------
%
%%
%input preconditioning
if exist('t0','var')==0
    t0=0;
end

if exist('option','var')==0
    option=[0,0];% default setting: do not add points, perform the refinement.
end


if exist('order','var')==0
    order=2;
end

if size(order,2)~=1
    order=order';%The order suppose to be a vertical vector.
end
if numel(order)<4
    order=[order ;order(end)*ones(4-numel(order),1)];
    % fill the missing parameter automatically
    %     curve_spline_order=order(1);
    %     DR_spline_order=order(2);
    %     RKorder = order(3);
    %     cubatureOrder=order(4);
end




%%
%construct generating curve.
DonatingRegion = donatingRegion(curve, vel,t0,te ,nSeg,order(1:3),option);

%%
%integrate pfun over the generating curve.
splType = 'not-a-knot';%4th order spline type
cubature_type=4;% guass legendre.
[xNodes, yNodes, weights] = splinegauss(order(4), DonatingRegion.DR,...
    DonatingRegion.spline_order_vett,  splType,cubature_type);
if isempty(weights)
    flux=0;
    return
end

fNodes = pfun(xNodes, yNodes,t0);
flux = weights'*fNodes;

end




