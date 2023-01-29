function x = chebyshev_nodes(n)


%====================================================================================================
% Panagiotis Veneris, U of Liverpool, November 2022

% Function that: creates Chebyshev nodes (gridpoints) at which we will
% approximate V(k). The x's are the roots of the Chebyshev polynomial, and
% they are called "collocation points"

% Input:  n, number on nodes; corresponds to the order of the Chebyshev
%         polynomial specified in main m.file

% Output: x: (m x 1) vector of nodes {x_i}_i=1,...m in [-1,1] at which we will approximate V(k),
%         i.e nodes (gridpoints) at which we will evaluate our Cheby polynomials T_j
%====================================================================================================

m = n+1; % # nodes at which we will approximate V(k): since m = n+1 then Chebyshev interpolation
x = -cos( (pi/(2*m)) .* (2.*(1:m)' -1));

end

