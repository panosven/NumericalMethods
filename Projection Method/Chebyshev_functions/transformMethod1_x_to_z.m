function z = transformMethod1_x_to_z(x,p) % x,m,p : inputs, k : output

%=========================================================================================================================
% Panagiotis Veneris, U of Liverpool, November 2022

% Map Chebyshev nodes x_i from [-1,1] to the domain of the gridpoints for productivity [p.zmin,p.zmax]

% In other words, it rescales Chebyshev nodes from [-1,1] to [p.zmin,p.zmax]
% In other words, it projects collocation points x in the z space (this
% explanation comes from Villaverde - Caldara code)

% Input:  x:  nx1 vector of nodes {x_i}_i=1,...m at which we will approximate V(z),
%             i.e nodes (gridpoints) at which we will evaluate our Chebyshevcpolynomials T_x_z
%         m:  # nodes at which we will approximate V(k)
%         p:  struct array of model parameters 

% Output: z: {z_i}_i=1...m in [p.zmin,p.zmax]

% Note: in this case, p.zmin, p.zmax will NOT be one of the gridpoints
%==========================================================================================================================

% map the nodes x to z
z = (1 + x) * ((p.zmax-p.zmin)/2) + p.zmin; % z(x): function that maps the nodes x (x in [-1,1]) to the domain of the gridpoints for productivity (z in [p.zmin,p.zmax])

end
