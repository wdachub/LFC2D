
function [nodes_x,nodes_y,weights]=splinegauss(N,control_points,...
    spline_order_vett,SPLtypestring,cubature_type)

%--------------------------------------------------------------------------
% INPUTS.
%--------
%
% [N]      : DEGREE OF PRECISION OF THE ALGEBRAIC CUBATURE RULE.
%
% [control_points]: IF THE POLYGON HAS "L" SIDES, "control_points" IS A
%           VARIABLE CONTAINING ITS VERTICES, ORDERED COUNTERCLOCKWISE.
%           AS LAST ROW MUST HAVE THE COMPONENTS OF THE FIRST VERTEX.
%           IN OTHER WORDS, THE FIRST ROW AND LAST ROW ARE EQUAL.
%           "control_points" IS A "L+1 x 2" MATRIX.
%           IN CASE OF CURVILINEAR POLYGONS, "control_points" ARE
%           VERTICES OF THE CURVILINEAR POLYGON.
%
% [rotation]:SUPPOSE x_min, x_max, y_min, y_max ARE THE MINIMUM AND MAXIMUM
%           VALUES IN x AND y REACHED BY THE POLYGON, I.E.
%           "R=[x_min,x_max] x [y_min,y_max]" IS THE SMALLEST RECTANGLE
%           WITH SIDES PARALLEL TO THE AXIS x AND y, CONTAINING THE POLYGON.
%
%           [rotation=0]: GENERAL CASE, BUT OBSERVE THAT THE FUNCTION MUST
%                         BE DEFINED IN THE RECTANGLE "R" DESCRIBED ABOVE.
%           [rotation=1]: GOOD FOR CONVEX POLYGONS. WITH THIS CHOICE, FOR
%                         CONVEX POLYGONS IT SUFFICES THAT THE FUNCTION IS
%                         DEFINED IN THE POLYGON. FOR CURVILINEAR POLYGONS
%                         IT PROVIDES A GOOD REFERENCE SEGMENT BUT DOES NOT
%                         INSURE THAT ALL THE CUBATURE NODES ARE INSIDE THE
%                         DOMAIN.
%           [rotation=2]: THE USER CHOOSES A SPECIAL REFERENCE SEGMENT "PQ".
%                         CHECK [1] FOR FURTHER DETAILS. SEE THE VARIABLES
%                         "P", "Q" BELOW.
%
% [P, Q]   : IF [rotation=2] THEN THE ALGORITHM CHOOSES A PREFERRED SEGMENT
%            "PQ" HAVING "P", "Q" AS EXTREMA, AS SET BY THE USER.
%            "P" AND "Q" ARE "1 x 2" ROW VECTORS. (SEE [1] FOR DETAILS).
%
% [spline_order_vett]: SPLINE ORDER AND BLOCK (VECTOR).
%
%            [spline_order_vett(i,1)=2] : PIECEWISE LINEAR SPLINES.
%            [spline_order_vett(i,1)=4] : CUBIC SPLINES
%                         (DEPENDING ON "SPLtypestring" IN INPUT).
%            [spline_order_vett(i,1)=k] : IN THE CASE k IS NOT 2 OR 4 IT
%                          CHOOSES THE k-TH ORDER SPLINES (FOR ADDITIONAL
%                          HELP DIGIT "help spapi" IN MATLAB SHELL).
%
%            [spline_order_vett(:,2)]: VECTOR OF FINAL COMPONENTS OF A BLOCK.
%
%            EXAMPLE:
%
%            "spline_order_vett=[2 31; 4 47; 8 67]" MEANS THAT FROM THE 1st
%             VERTEX TO THE 31th VERTEX WE HAVE AN ORDER 2 SPLINE (piecewise
%             linear), FROM THE 32th VERTEX TO THE 47th WE USE A 4th ORDER
%             SPLINE (i.e. A CUBIC AND PERIODIC SPLINE BY DEFAULT), FROM
%             THE 48th TO THE 67th (AND FINAL!) WE USE AN 8th ORDER SPLINE.
%
% [cumulative]: IT CHOOSES THE PARAMETRIZATION.
%            [cumulative=0]: 1:N.
%            [cumulative=1]: CUMULATIVE.
%
% [SPLtypestring]: IF [spline_order_vett=4] IT DECIDES THE TYPE OF END
%             CONDITIONS OF THE CUBIC SPLINE. IT IS A STRING. IT CAN BE:
%
%             'complete'   : match endslopes (as given in VALCONDS, with
%                     default as under *default*).
%             'not-a-knot' : make spline C^3 across first and last interior
%                     break (ignoring VALCONDS if given).
%             'periodic'   : match first and second derivatives at first
%                     data point with those at last data point (ignoring
%                     VALCONDS if given).
%             'second'     : match end second derivatives (as given in
%                    VALCONDS, with default [0 0], i.e., as in variational).
%             'variational': set end second derivatives equal to zero
%                     (ignoring VALCONDS if given).
%
% [cubature_type]:
%
%            [cubature_type=0]: PADUA POINTS.
%            [cubature_type=1]: FEJER 1 TYPE (TENSORIAL).
%            [cubature_type=2]: FEJER 2 TYPE (TENSORIAL).
%            [cubature_type=3]: CLENSHAW-CURTIS (TENSORIAL).
%            [cubature_type=4]: GAUSS-LEGENDRE (TENSORIAL).
%            [cubature_type=5]: GAUSS-LEGENDRE-LOBATTO (TENSORIAL).
%
%----------
% OUTPUTS.
%----------
%
% [nodes_x, nodes_y]: THE CUBATURE RULE PRODUCES THE NODES
%                          "(nodes_x,nodes_y)".
%            "nodes_x" AND "nodes_y" ARE COLUMN VECTORS.
%
%
% [weights]: CUBATURE WEIGHT. COLUMN VECTOR.
%
%--------------------------------------------------------------------------
% REFERENCE PAPERS.
%-------------------
%
% [1]. A. SOMMARIVA and M. VIANELLO "Gauss-like and triangulation-free
%      cubature over polygons". BIT Numerical Methematics 47 (2007),
%      441-453.
%
% [2]. A. SOMMARIVA and M. VIANELLO "Gauss-Green cubature an moment
%      computation over arbitrary geometries". JCAM.
%
%--------------------------------------------------------------------------


% Notation: "xi" as in the paper.
% xi=min(control_points(:,1));
xi=median(control_points(:,1));


%--------------------------------------------------------------------------
% COMPUTE NODES.
%--------------------------------------------------------------------------

nodes_x=[];
nodes_y=[];
weights=[];

% Number of Blocks.
L=size(spline_order_vett,1);

for block_index=1:L
    spline_block_order=spline_order_vett(block_index,1);
    if spline_block_order==1
        continue
    end
    spline_block_degree=spline_block_order-1;
    
    % Initial and final indices of "control points" in the block.
    if (block_index ==1)
        initial_point_index=1;
    else
        initial_point_index=spline_order_vett(block_index-1,2);
    end
    final_point_index=spline_order_vett(block_index,2);
    
    % Spline order in the block.
    
    % Control points (x_loc,y_loc) in the block.
    pts_loc=control_points(initial_point_index:final_point_index,:);
    
    % Parametrical description of the block.
    
    s_loc=[0;cumsum(edgeLength(pts_loc))];
    
    
    
    
    
    % Computing the spline parametrical description of the block.
    % "ppx", "ppy" describe S_i1, S_i2 in the paper, while "ppy1"
    % describe S'_i2.
    
    switch spline_block_order
        case 2
        case 4
            
            % CUBIC SPLINES BY CSAPE. AS DEFAULT WE USE PERIODIC CUBIC SPLINES.
            % Derivatives parameters are computed as well.
            ppx=csape(s_loc,pts_loc(:,1),SPLtypestring);
            ppy=csape(s_loc,pts_loc(:,2),SPLtypestring);
            [breaks_y,coeffs_y]=unmkpp(ppy);
            N_y=size(coeffs_y,1);
            dcoeffs_y=[zeros(N_y,1) 3*coeffs_y(:,1) 2*coeffs_y(:,2) ...
                coeffs_y(:,3)];
            ppy1=mkpp(breaks_y,dcoeffs_y);
            
        otherwise
            
            ppx=spapi(spline_block_order,s_loc,pts_loc(:,1));
            ppy=spapi(spline_block_order,s_loc,pts_loc(:,2));
            ppy1=fnder(ppy,1);
            
    end
    
    
    
    % Every block is subdivided in "number_of_subblocks" curves determined
    % by successive control points.
    number_of_subblocks=final_point_index-initial_point_index;
    
    % Cubature rule on the square [-1,1] x [-1,1].
    % Padua Points: 0, Gauss-Legendre: 4.
    [x_pts, y_pts, wpd]=cubature_manager(N,spline_block_degree,...
        cubature_type);
    
    
    % Computing quadrature points from a general sub-block. The cases in
    % which the order is 2 is a little different from other spline orders.
    % Consequently, we distinguish between them.
    for index_control_point=1:number_of_subblocks
        
        if (spline_block_order == 2)
            
            x1=pts_loc(index_control_point,1);
            x2=pts_loc(index_control_point+1,1);
            
            y1=pts_loc(index_control_point,2);
            y2=pts_loc(index_control_point+1,2);
            
            if ~(x2 == xi && x1 == xi)
                if ~( (y2-y1) == 0)
                    
                    % Computing nodes.
                    
                    qij_u= (y_pts +1)/2;
                    
                    si1_ij=x1+(x2-x1)*qij_u;
                    si2_ij=y1+(y2-y1)*qij_u;
                    
                    x=((si1_ij-xi)/2).*x_pts + (si1_ij+xi)/2;
                    y=si2_ij;
                    
                    
                    nodes_x=[nodes_x; x];
                    nodes_y=[nodes_y; y];
                    
                    % Computing weights via formula (18).
                    diff_s=s_loc(index_control_point+1)-...
                        s_loc(index_control_point);
                    diff_y=y2-y1;
                    S1_i2_q_ij_u=diff_y/diff_s;
                    
                    initial_scaling_factor=(s_loc(index_control_point+1)...
                        -s_loc(index_control_point))/4;
                    scaling_fact_minus=si1_ij-xi;
                    local_weights=initial_scaling_factor.*...
                        scaling_fact_minus.*S1_i2_q_ij_u.*wpd;
                    
                    weights=[weights; local_weights];
                    
                end
            end
            
        else
            % spline_block_order ~= 2.
            
            %-----------------------------------------
            % COMPUTING NODES.
            %-----------------------------------------
            
            % Computing "q_ij_u" (see formula (19)).
            
            half_pt_s=(s_loc(index_control_point+1)+...
                s_loc(index_control_point))/2;
            half_length_s=(s_loc(index_control_point+1)-...
                s_loc(index_control_point))/2;
            q_ij_u=half_pt_s+half_length_s*y_pts;
            
            % Evaluating S_i1(q_ij(u)), S_i2(q_ij(u)).
            switch spline_block_order
                case 4
                    S_i1_q_ij_u=ppval(ppx,q_ij_u);
                    S_i2_q_ij_u=ppval(ppy,q_ij_u);
                    
                otherwise
                    S_i1_q_ij_u=fnval(ppx,q_ij_u);
                    S_i2_q_ij_u=fnval(ppy,q_ij_u);
            end
            
            
            % See formula (18). Terms involving the spline in the arguments
            % of "f".
            scaling_fact_plus=(S_i1_q_ij_u+xi)/2;
            scaling_fact_minus=(S_i1_q_ij_u-xi)/2;
            
            x=scaling_fact_minus.*x_pts+scaling_fact_plus;
            y=S_i2_q_ij_u;
            
            nodes_x=[nodes_x; x];
            nodes_y=[nodes_y; y];
            
            %-----------------------------------------
            % COMPUTING WEIGHTS.
            %-----------------------------------------
            
            switch spline_block_order
                case 4
                    
                    S1_i2_q_ij_u=ppval(ppy1,q_ij_u);
                    % SCALING FACTOR: WE PASS FROM INDEXES [i,i+1].
                    initial_scaling_factor=(s_loc(index_control_point+1)-...
                        s_loc(index_control_point))/4;
                    
                otherwise
                    
                    S1_i2_q_ij_u=fnval( fnder( ppy, 1), q_ij_u );
                    % SCALING FACTOR: WE PASS FROM INDEXES [i,i+1].
                    initial_scaling_factor=(s_loc(index_control_point+1)-...
                        s_loc(index_control_point))/4;
                    
            end
            
            local_weights=initial_scaling_factor.*...
                (2*scaling_fact_minus).*S1_i2_q_ij_u.*wpd;
            
            weights=[weights; local_weights];
        end
    end
    
end


zeroweights=weights==0;
nodes_x(zeroweights)=[];
nodes_y(zeroweights)=[];
weights(zeroweights)=[];

end