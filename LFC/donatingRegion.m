function [DonatingRegion, TimeCurve]= donatingRegion(curve, vel,t0,te ,nSeg,order,option)
% construct the donating region
%  associated to the curve with endpoint p1,p2 (or L, N).
%  within time interval [t0, te] under the velocity field vel.
% =========================================================================
% INPUTS
% =========================================================================
%       Name  |         Description
%-------------|------------------------------------------------------------
% SAcure      | data format:
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
%             |  1. [order] only a number, that means the order of spline,
%             |  runge kutta are same.
%             |  2. [curve_spline_order, RKorder]
% ------------|------------------------------------------------------------
%   nSeg      | the number of segments to divide the curve, so there is nSeg+1
%             | nodes on the curve.
%             | In order to implement spline interpolation, nSeg must be
%             | at least order-1;
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
% ------------|------------------------------------------------------------
% =========================================================================
% OUTPUTS
% =========================================================================
%       Name  |         Description
% ------------|------------------------------------------------------------
%   pts       | points on the donating region boundary.
%             | format: LN, streakline2, timeline, streakline1,
%             |
%             |
% ------------|------------------------------------------------------------
%  vertIDs    | indices marking the boundary of LN, its preimage and streaklines.
%             |
% ------------|------------------------------------------------------------





% enforce preconditions.
if (t0>te)
    error('t0 must less than te!');
end

if numel(order)<3
    if size(order,2)~=1
        order=order';
    end
    % if input order less than 2, increase them to 2; It is meaningless if
    %spline_order£¬RKorder is less than 2.
    order(order<2)=2;
    % fill the missing parameter automatically
    order=[order ;order(end)*ones(3-numel(order),1)];
end
curve_spline_order=order(1);
DR_spline_order=order(2);
RKorder = order(3);



if exist('option','var')==0
    option=[1,0]; %will not add points
    %will not use the adaptive mesh
end
if numel(option)<2
    option=[option ;option(end)*zeros(2-numel(option),1)];%add zero at end so that the length of option becomes 2.
end
add_point=option(1);
refine=option(2)>0;
if size(curve,1)==2
    add_point=0;%It is not necessary to add point when curve is a line segment. 
end
%%
global RK
switch RKorder %RKorder
    case 2       % Modified Euler
        RK.a = zeros(2,2);
        RK.a(2,1) = 1/2;
        RK.b = [0 1];
        RK.c = [0 1/2];
    case 3   % Heun's formula
        RK.a = zeros(3,3);
        RK.a(2,1) = 1/3;
        RK.a(3,2) = 2/3;
        RK.b = [1 0 3]/4;
        RK.c = [0 1 2]/3;
    case 4
        RK.a = zeros(4,4);
        RK.a(2,1) = 1/2;
        RK.a(3,2)=1/2;
        RK.a(4,3)=1;
        RK.b = [1/6, 1/3, 1/3, 1/6];
        RK.c = [0 1/2 1/2 1];
    case 5   % Verner 1978 (DVERK)
        RK.a = zeros(6,6);
        RK.a(2,1)   = [1/6];
        RK.a(3,1:2) = [4/75 16/75];
        RK.a(4,1:3) = [5/6  -8/3  5/2];
        RK.a(5,1:4) = [-165/64 55/6 -425/64 85/96];
        RK.a(6,1:5) = [12/5 -8 4015/612 -11/36 88/255];
        RK.b = [13/160 0 2375/5984 5/16 12/85 3/44];
        RK.c = [0 1/6 4/15 2/3 5/6 1];
    case 6   % Verner 1978 (DVERK)
        RK.a = zeros(8,8);
        RK.a(2,1)   = [1/6];
        RK.a(3,1:2) = [4/75 16/75];
        RK.a(4,1:3) = [5/6  -8/3  5/2];
        RK.a(5,1:4) = [-165/64 55/6 -425/64 85/96];
        RK.a(6,1:5) = [12/5 -8 4015/612 -11/36 88/255];
        RK.a(7,1:6) = [-8263/15000 124/75 -643/680 -81/250 2484/10625 0];
        RK.a(8,1:7) = [3501/1720 -300/43 297275/52632 -319/2322 24068/84065 0 3850/26703];
        RK.b = [3/40 0 875/2244 23/72 264/1955 0 125/11592 43/616];
        RK.c = [0 1/6 4/15 2/3 5/6 1 1/15 1];
    case 7   % Prince and Dormand 1981
        RK.a = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0;
            29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
            16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
            39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
            246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
            -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
            185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
            403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0];
        RK.b = [ 13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321, 656045339/265891186,   -3867574721/1518517206,   465885868/322736535,  53011238/667516719,                  2/45,    0];
        RK.c = [ 0, 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1];
    case 8   % Prince and Dormand 1981
        RK.a = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0;
            29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
            16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
            39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
            246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
            -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
            185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
            403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0];
        RK.b = [ 14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731,   561292985/797845732,   -1041891430/1371343529,  760417239/1151165299, 118820643/751138087, -528747749/2220607170,  1/4];
        RK.c = [ 0, 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1];
    otherwise
        error('invalid order!');
end

%%
%spline approximation

if isa(curve, 'function_handle')
    add_point=0;% we don't need adding points into sampling set, since curve function can give exact value.
    p1=curve(0,t0);%L
    p2=curve(1,t0);%N
    %     SAcurve=curve;
else
    if size(curve,2)~=2
        curve=curve';
        warning('Data format of curve is not correct, which should be N*2 matrix. Input has been transposed automatically.')
    end
    p1=curve(1,:);%L
    p2=curve(end,:);%N
    SAcurve.ct=[0;cumsum(edgeLength(curve))];
    curveLength=SAcurve.ct(end);
    SAcurve.ct=SAcurve.ct/curveLength;%SAcurve.ct
    
    switch curve_spline_order
        % case 2
        % this case is hard to write, so i include this into the 'otherwise' case.
        % polygon approximation will be implement by the build-in function spapi.
        case 4
            % csape is faster than spapi but only work when spline order is 4.
            SAcurve.ppx=csape(SAcurve.ct,curve(:,1),'not-a-knot');
            SAcurve.ppy=csape(SAcurve.ct,curve(:,2),'not-a-knot');
        otherwise
            SAcurve.ppx=spapi(order(1),SAcurve.ct,curve(:,1));
            SAcurve.ppy=spapi(order(1),SAcurve.ct,curve(:,2));
    end
end






%%
%The generating curve is a closed curve and DonatingRegion.p is parameter of it.
%DonatingRegion.p=0 <=> L
%DonatingRegion.p=1 <=> N
%DonatingRegion.p=2 <=> preimage of N
%DonatingRegion.p=3 <=> preimage of L
%DonatingRegion.p=4 <=> L
%Hence we have,
%DonatingRegion.p=0~1 <=> LN
%DonatingRegion.p=1~2 <=> streakline of N
%DonatingRegion.p=2~3 <=> preimage of LN
%DonatingRegion.p=3~4 <=> streakline of L
%There are two advantage of using DonatingRegion.p:
%We can locate point on the point on the generating curve easily by
%checking DonatingRegion.p. The corresponding rule have shown above.
%We can apply this program to moving curve case without changing anything.



%TimeCurve_pvertex and TimeCurve_fvertex represent the L, N and preimage of
%them. Those two variables will not change during the program.
%p is independent variable, and f is dependent variable.
% They are a bridge between generating curve and DonatingRegion.p.
% DonatingRegion.p->TimeCurve.f->generating curve
% where TimeCurve.f is computed by interpolation of TimeCurve_pvertex and
% TimeCurve_fvertex.
global TimeCurve_pvertex TimeCurve_fvertex
TimeCurve_pvertex =[0,1,2,3,4];
% if we don't normalize SAcurve.ct, we have
% TimeCurve_fvertex =[0,t0;
%     SAcurve.ct(end),t0;
%     SAcurve.ct(end),te;
%     0,te;
%     0,t0];
TimeCurve_fvertex =[0,t0;
    1,t0;
    1,te;
    0,te;
    0,t0];


%If the curve is close, we don't need calculate the streakline.
if nnz(p1==p2)==2
    DonatingRegion.p=[linspace(0,1,nSeg+1),linspace(2,3,nSeg+1),4]';
else
    DonatingRegion.p=linspace(0,4,4*nSeg+1)';
end


if isa(curve, 'numeric')
    if add_point==1
        DonatingRegion.p=sort(unique([DonatingRegion.p;SAcurve.ct;3-SAcurve.ct]));
        %They are useful when construct spline_order_vett.
        indTCp1=find((DonatingRegion.p==1));%index of N
        indTCp2=find((DonatingRegion.p==2));%index of preimage of N
    else
        indTCp1=nSeg+1;
        indTCp2=nSeg+2;
    end
end

%the basic parameter
%TimeCurve.f(:,1) space or say curve parameter,
%TimeCurve.f(:,2) time.
TimeCurve.f=[
    interp1(TimeCurve_pvertex , TimeCurve_fvertex (:,1),DonatingRegion.p), ...
    interp1(TimeCurve_pvertex ,TimeCurve_fvertex (:,2),DonatingRegion.p)];


if isa(curve, 'function_handle')
    DonatingRegion.space=curve(TimeCurve.f(:,1),TimeCurve.f(:,2) ) ;
else
    switch curve_spline_order
        case 4
            DonatingRegion.space=[ppval(SAcurve.ppx, TimeCurve.f(:,1)) ppval(SAcurve.ppy, TimeCurve.f(:,1))] ;
            %the curve is not moving, so the position is not depend on time.
        otherwise
            DonatingRegion.space=[fnval(SAcurve.ppx, TimeCurve.f(:,1)) fnval(SAcurve.ppy, TimeCurve.f(:,1))] ;
    end
end




dt = (te-t0)/nSeg;
%compute the DR.
DonatingRegion.DR=flowmap(DonatingRegion.space, TimeCurve.f(:,2),t0 , vel, -dt);
%Delete same points, which will happen when velocity is zero. Recall the
%orginal point in the intersection test.
% [DonatingRegion.DR,DonatingRegion.p] = DeleteSamePts(DonatingRegion.DR,DonatingRegion.p,SAcurve.ct,curveLength/nSeg/100);
[DonatingRegion.DR,DonatingRegion.p] = DeleteSamePts(DonatingRegion.DR,DonatingRegion.p,SAcurve.ct,eps);
%%
%refine DR
if refine==1
    if isa(curve, 'function_handle')
        indTCp1=find((DonatingRegion.p==1));
        indTCp2=find((DonatingRegion.p==2));
        curveLength=edgeLength(DonatingRegion.DR(1:indTCp1,:));%length of curve
        curveLength=curveLength(end);
    end
    
    if option(2)==1
       maxDist=curveLength/nSeg;
    else
        maxDist=option(2);
    end
    % maxDist=1110;
    % the length of all edges
    
    if nnz(p1==p2)==2
        if curve_spline_order==2%program optimization
            %when curve is line segment, it is not neccessary to refine it.
            DRtemp1=DonatingRegion.DR(1:indTCp1,:);
            TCtemp1=DonatingRegion.p(1:indTCp1,:);
        else
            %get more sampling points of LN.
            [DRtemp1,TCtemp1]=AdaptiveBackAdectPts(DonatingRegion.DR(1:indTCp1,:),DonatingRegion.p(1:indTCp1,:),t0, vel,dt,maxDist,order(1),SAcurve);
        end
        [DRtemp2,TCtemp2]=AdaptiveBackAdectPts(DonatingRegion.DR(indTCp2:end-1,:),DonatingRegion.p(indTCp2:end-1,:),t0, vel,dt,maxDist,order(1),SAcurve);%timeline
        DonatingRegion.DR=[DRtemp1;DRtemp2;p1];
        DonatingRegion.p=[TCtemp1;TCtemp2;4];
    else
        
        [DonatingRegion.DR,DonatingRegion.p]=AdaptiveBackAdectPts(DonatingRegion.DR,DonatingRegion.p,t0, vel,dt,maxDist,order(1),SAcurve);
    end
end

TimeCurve.f=[
    interp1(TimeCurve_pvertex , TimeCurve_fvertex (:,1),DonatingRegion.p), ...
    interp1(TimeCurve_pvertex ,TimeCurve_fvertex (:,2),DonatingRegion.p)];

%%
%construct spline_order_vett, which provide the infomation of spline order
%to integrator.
if curve_spline_order==2&& DR_spline_order~=2 && add_point==1
    % when curve_spline_order is 2 and DR_spline_order is larger than 2,
    % the preimage of ploygon vetex will become non-differentiable point on
    % the preiamge of curve. 
    [~,vertIDs]=ismember([1;3-flipud(SAcurve.ct);4],DonatingRegion.p);
    DonatingRegion.spline_order_vett=[order(2)*ones(size(vertIDs)) vertIDs];
else
    %if curve and DR order is 2
    %if curve order is large than 2, so the curve is smooth.
    [~,vertIDs]=ismember([1;2;3;4],DonatingRegion.p);
    DonatingRegion.spline_order_vett=[order(2)*ones(size(vertIDs)) vertIDs];
end
%the spline order of curve

DonatingRegion.spline_order_vett(1,1)=curve_spline_order;

% DonatingRegion.spline_order_vett(3,1)=order(1);
%we make an assumption here, the spline order of preimage of curve is equal
%to the spline order of curve.
%questionable, is that best assumption?
%It is useful when the curve is not smooth, since the flow mapping will not
%increase the regularity.






%if the points on a block is less than the given order of spline,
%set the order of spline on this block to be the number of points.
%so the order is not less than 1.When order is 1, the integrator will skip
%this block.
tmp=diff([1;DonatingRegion.spline_order_vett(:,2)])+1;
DonatingRegion.spline_order_vett(tmp<DR_spline_order,1)=tmp(tmp<DR_spline_order);

if nnz(p1==p2)==2
    DonatingRegion.spline_order_vett([2,end],1)=[1;1];
    %the integral of spline is offset. we don't need compute them.
end

if nnz(DonatingRegion.spline_order_vett<0)
    
end 
end

function [pts, TimeCurve]=AdaptiveBackAdectPts(pts,TimeCurve,t0, vel,dt,maxDist,spline_order,SAcurve)
global TimeCurve_pvertex TimeCurve_fvertex
if maxDist==0
    %maxDist==0 will case endless loop, so stop the function immediately.
    return;
end

elTL = edgeLength(pts);
% find those edges with length bigger than maxDist
% it may spend many time on distribution memory if n is large.
idx = find(elTL>maxDist);
while ~isempty(idx)
    newTimeLn = pts(1:idx(1),:);
    newTimeSpace = TimeCurve(1:idx(1),:);
    for i=1:length(idx)
        eid = idx(i);  % the index into elTL
        nSubEdges = ceil(elTL(eid)/maxDist);
        toAdd_parameter=TimeCurve(eid,:)+bsxfun(@times,(1:nSubEdges-1)'/nSubEdges,(TimeCurve(eid+1,:)-TimeCurve(eid,:)));
        newTimeSpace = [newTimeSpace; toAdd_parameter];
        if i==length(idx)
            endId = length(TimeCurve);
        else
            endId = idx(i+1);
        end
        newTimeSpace = [newTimeSpace; TimeCurve(eid+1:endId,:)];
        toAdd_ts=[interp1(TimeCurve_pvertex , TimeCurve_fvertex (:,1),toAdd_parameter), ...
            interp1(TimeCurve_pvertex ,TimeCurve_fvertex (:,2),toAdd_parameter)];
        
        
        
        
        if isa(SAcurve, 'function_handle')
            toAdd=SAcurve(toAdd_ts(:,1),toAdd_ts(:,2) ) ;
        else
            switch spline_order
                case 4
                    toAdd=[ppval(SAcurve.ppx, toAdd_ts(:,1)) ppval(SAcurve.ppy, toAdd_ts(:,1))] ;
                    %the curve is not moving, so the position is not depend on time.
                otherwise
                    toAdd=[fnval(SAcurve.ppx, toAdd_ts(:,1)) fnval(SAcurve.ppy, toAdd_ts(:,1))] ;
            end
        end
        
        newTLpts =  flowmap(toAdd, toAdd_ts(:,2),t0, vel, -dt);
        newTimeLn = [newTimeLn; newTLpts; pts(eid+1:endId,:)];
    end
    % update p12static and timeline for the next loop.
    TimeCurve = newTimeSpace;
    pts = newTimeLn;
    elTL = edgeLength(pts);
    idx = find(elTL>maxDist);
end
end


function [pts,TimeCurve] = DeleteSamePts(pts,TimeCurve,SAcurvect,minDist)
% delete the adjacent same points in the polygonal representation of curve.
% =========================================================================
% INPUTS
% =========================================================================
%       Name  |         Description
%-------------|------------------------------------------------------------
% pts         | vector which record the coordinate of points on curve.
%             | data format:
%             | [ x1 x2 ...xn
%             |   y1 y2 ...yn]
%             | if x1=x2 and y1=y2, [x2, y2] will be dropped.
% ------------|------------------------------------------------------------
%   minDist   |  When the distance between adjacent points is less the
%             |  minDist, they will be consider as same point.
%             |
%             |  minDist will be the machine number if input is missing.
% ------------|------------------------------------------------------------
% =========================================================================
% OUTPUTS
% =========================================================================
%       Name  |         Description
% ------------|------------------------------------------------------------
%   pts       | polygonal representation of curve with no
%             | adjacent same points.
%             |
% ------------|------------------------------------------------------------
%%
%input preconditioning
if exist('minDist','var')==0
    minDist=eps;%let minDist be the machine number if input is missing.
end
%%
i=1;
while i<size(pts,1) %execute the loop if i less than the number of points.
    if norm(pts(i,:)-pts(i+1,:))<minDist %the adjacent points are same.
        %ie. the distance between them is less the minDist.
        %if sum((pts(:,i)-pts(:,i+1)).^2)<eps
%         if ~isempty(ismember([1,2,3,4],TimeCurve(i+1)))
        if nnz(ismember([1;3-flipud(SAcurvect);4],TimeCurve(i+1)))

            i=i+1; %those value is indicator of vectex of DR.
        else
            pts(i+1,:)=[];%delete the latter one.
            TimeCurve(i+1)=[];%delete the coresponding parameters.
        end
    else
        i=i+1;
    end
end

end
