
%--------------------------------------------------------------------------
% cubature_manager.
%--------------------------------------------------------------------------

function [nodes_x, nodes_y, weights]=cubature_manager(M,p_i,cubature_type)


%------------------ --------------------------------------------------------
% OBJECT:
%-----------
%
% THE TARGET OF THIS ROUTINE IS TO PROVIDE NODES AND WEIGHTS OF A
% QUADRATURE ROUTINE ON [-1,1]. THE CODES ARE DESCRIBED BY L.N. TREFETHEN 
% IN HIS CLENSHAW-CURTIS PAPER AND BY WALDVOGEL (PUBLISHED BY BIT).
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% M : DEGREE OF PRECISION OF 2D CUBATURE RULE.
%
% p_i: SPLINE DEGREE
%
% cubature_type: IT POINTS A CUBATURE ON THE SQUARE [-1,1] x [-1,1].
%
%     [cubature_type=0]: PADUA POINTS.
%     [cubature_type=1]: FEJER 1 (TENSORIAL).
%     [cubature_type=2]: FEJER 2 (TENSORIAL).
%     [cubature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS) (TENSORIAL).
%     [cubature_type=4]: GAUSS-LEGENDRE (TENSORIAL).
%     [cubature_type=5]: GAUSS-LEGENDRE-LOBATTO (TENSORIAL).
%
%----------
% OUTPUTS:
%----------
%
% nodes : M x 2 MATRIX OF NODES, IN THE INTERVAL [-1,1].
% 
% weights: M x 1 COLUMN VECTOR OF WEIGHTS.
%
%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES :
%----------------------
%
% 1. quadrature_rules_1D
% 2. pdcub
%
%--------------------------------------------------------------------------


    
    % lower degree of a Gauss-Legendre rule with ADE equal to M.
    Ntau=ceil((M+1)/2); 
    % In the quadrature rule, for M as input we have a rule with M+1 pts.
    Ntau_cub=Ntau-1; 
    [nodes_tau,weights_tau]=quadrature_rules_1D(Ntau_cub,cubature_type);
    
    % lower degree of a Gauss-Legendre rule with ADE equal to Nu_deg.
    Nu=ceil((M+2)*p_i/2); 
    Nu_cub=Nu-1;
    [nodes_u,weights_u]=quadrature_rules_1D(Nu_cub,cubature_type);
    
    [nodes_x,nodes_y]=meshgrid(nodes_tau,nodes_u);
    nodes_x=nodes_x(:);
    nodes_y=nodes_y(:);
    
    [weights1,weights2]=meshgrid(weights_tau,weights_u);
    weights=weights1(:).*weights2(:);
    

    
end





function [nodes,weights]=quadrature_rules_1D(n,quadrature_type)

%--------------------------------------------------------------------------
% OBJECT:
%-----------
% THE TARGET OF THIS ROUTINE IS TO PROVIDE NODES AND WEIGHTS OF A
% QUADRATURE ROUTINE ON [-1,1]. THE CODES ARE DESCRIBED BY L.N. TREFETHEN 
% IN HIS CLENSHAW-CURTIS PAPER AND BY WALDVOGEL (PUBLISHED BY BIT).
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% n : NUMBER OF NODES OF THE QUADRATURE RULE (NOT THE DEGREE!!).
%
% quadrature_type: IT POINTS A QUADRATURE RULE
%           [quadrature_type=1]: FEJER 1.
%           [quadrature_type=2]: FEJER 2.
%           [quadrature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS).
%           [quadrature_type=4]: GAUSS-LEGENDRE.
%           [quadrature_type=5]: GAUSS-LEGENDRE-LOBATTO.
%           [quadrature_type=6]: COMPOSITE TRAPEZOIDAL RULE.
%
%----------
% OUTPUTS:
%----------
%
% nodes : M x 2 MATRIX OF NODES, IN THE INTERVAL [-1,1].
% 
% weights: M x 1 COLUMN VECTOR OF WEIGHTS.
%
%--------------------------------------------------------------------------
% ADDITIONAL ROUTINES :
%----------------------
%
% 1. r_jacobi
% 2. gauss
% 3. lobatto_jacobi
%
% THESE ROUTINES ARE WRITTEN BY D. LAURIE AND W. GAUTSCHI, AND CAN BE FOUND
% IN W. GAUTSCHI HOMEPAGE. THEY ARE ATTACHED IN THIS FILE (SEE THE BOTTOM
% OF THE FILE). NO EXTERNAL FILE IS NEEDED.
%
%--------------------------------------------------------------------------
% EXAMPLE:
%----------
%
% >> [nodes,weights]=quadrature_rules_1D(5,5)
% 
% nodes =
%
%   -1.0000
%   -0.8302
%   -0.4688
%    0.0000
%    0.4688
%    0.8302
%    1.0000
%
%
% weights =
%
%    0.0476
%    0.2768
%    0.4317
%    0.4876
%    0.4317
%    0.2768
%    0.0476
%
%--------------------------------------------------------------------------
% RELATED PAPERS:
%-----------------
%
% [1] W. GAUTSCHI, ORTHOGONAL POLYNOMIALS AND QUADRATURE.
%
% [2] LLOYD N. TREFETHEN, IS GAUSS QUADRATURE BETTER THAN CLENSHAW?CURTIS?
%
% [3] J. WALDVOGEL, FAST CONSTRUCTION OF THE FEJER AND CLENSHAW-CURTIS
%     QUADRATURE RULES.
%
%--------------------------------------------------------------------------

% GAUSS LEGENDRE.
        beta=0.5./sqrt(1-(2*(1:n)).^(-2));
        T=diag(beta,1)+diag(beta,-1);
        [V,D]=eig(T);
        nodes=diag(D); 
%         [nodes,index]=sort(nodes); 
%         weights=2*V(1,index).^2;
        weights=2*V(1,:).^2;
end