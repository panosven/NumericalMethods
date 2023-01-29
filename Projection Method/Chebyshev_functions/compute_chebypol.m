function T = compute_chebypol(x,n)

%==================================================================================================================================================
% Panagiotis Veneris, U of Liverpool, November 2022

% Function that: computes Chebyshev polynomials up to order n / evaluates each Chebyshev polynomial j at all gridpoints (Chebyshev nodes) x_i:x_1,.....,x_m 
% Input:  n: number on nodes; corresponds to the order of the Chebyshev
%         polynomial specified in main m.file         
%         x: (m x 1) vector of Chebyshev nodes defined on [-1,1] (i.e nodes at which we evaluate the Chebyshev
%         polynomials) 

% Output: T: m x (n+1);  T_j(x_i) with (i,j): i=1,....,m and j=0,....,n (or j= 1,....,n+1)
%===================================================================================================================================================

m = length(x);
T = zeros(m,n+1);  % initialize matrix of chebyshev polynomial (we evaluate the 29th-order cheby polynomial at 30 gridpoints x), T: m x (n+1) matrix

% Option 1
% for j = 1:n
%     T(:,j+1) = cos(j.*arccos(x));
% end

% Option 2
T(:,1) = 1;   % T_0(x) = cos(0) = 1
T(:,2) = x;   % T_1(x) = cos(arccos(x)) = x

for j = 2:n
    T(:,j+1) = 2*x .* T(:,j) - T(:,j-1); % T_j(x_i): column j is the j-1 polynomial evaluated at each x_i:x_1,.....,x_n (rows)
end

end
