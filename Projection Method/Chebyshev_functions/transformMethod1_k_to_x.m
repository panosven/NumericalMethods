function x =  transformMethod1_k_to_x(k,p)

%==================================================================================================================
% Panagiotis Veneris, U of Liverpool, November 2022

% Maps capital nodes (gridpoints) k (k_i in [p.kmin,p.kmax]) to the domain of chebyshev nodes x (x_i in [-1,1])

% Input:  k:  {k_i}_i=1...m; 
%         m:  # nodes at which we will approximate V(k)
%         p:  struct array of model parameters 

% Output: x: {x_i}_i=1...m in [-1,1]

% Notes: in this case p.kmin, p.kmax (p.zmin, p.zmax if we have a 2nd state) will NOT be one of the gridpoints
%==================================================================================================================

 x = (2*(k-p.kmin))/(p.kmax-p.kmin) - 1;  % x(k) = {x_i}_i=1...m: function that maps capital gridpoints (k) to the domain [-1,1]

end 

