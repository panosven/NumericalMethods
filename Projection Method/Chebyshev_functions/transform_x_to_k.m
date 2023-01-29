function k = transform_x_to_k(x,m,p) % x,m,p : inputs, k : output

%=========================================================================================================================
% Panagiotis Veneris, U of Liverpool, November 2022

% Function that: maps the nodes x_i (x_i in [-1,1]) to the domain of the gridpoints for capital (k_i in [p.kmin,p.kmax])
% In other words, it rescales Chebyshev nodes from [-1,1] to [p.kmin,p.kmax]
% In other words, it projects collocation points x in the k space (this
% explanation comes from Villaverde - Caldara code)

% Input:  x:  nx1 vector of nodes {x_i}_i=1,...m at which we will approximate V(k),
%             i.e nodes (gridpoints) at which we will evaluate our Chebyshev polynomials T_x
%         m:  # nodes at which we will approximate V(k)
%         p:  struct array of model parameters 

% Output: k: {k_i}_i=1...m in [p.kmin, p.kmax]
%==========================================================================================================================

% map the nodes x to k
k = (1 + x .* sec(pi/(2*m))) * ((p.kmax - p.kmin)./2) + p.kmin;  % k(x): function that maps the nodes x (x in [-1,1]) to the domain of the gridpoints for capital (k in [p.kmin,p.kmax])

end
