function poly = polynomial_comb(vX1, vX2)

%%======================================================================================
% Panagiotis Veneris, U of Liverpool, November 2022

% Since we want to compute a multidimensional (2d) approximation of a function,
% we will do so using Complete Polynomials Basis (see Collard's notes pg.45)

% Complete polynomials bases take products of order lower than a priori given 'kappa' 
% into account, ignoring higher terms of higher degrees.

% Creates polynomial combinations of inputs
%           - vX1: 1 x nz vector, row of Cheb. polynomials for productivity
%           - vX2: 1 x nk vector, row of Cheb. polynomials for capital
% Output:   - poly
% vector

% This code is a direct implementation of Mohammed Ait Lahcen Python code 
% on Quant Econ notes, Chebyshev approximation
%%=====================================================================================

% Retrieve dimensions
nz = size(vX1,2);  % return the  length of the second dimension (2) of matrix vX1 (order of Chebyshev polynomial+1)
nk = size(vX2,2);  % return the  length of the second dimension (2) of matrix vX2 (order of Chebyshev polynomial+1)      

kappa = max(nz,nk); % find the maximum length
poly = [];

for ik = 1:nk
    for iz = 1:nz
        if (ik+iz-2)<= (kappa-1) % We subtract 2 is because nz, nk are the # of Cheby nodes for each of productivity 
                                 % and capital, while we want the Chebyshev polynomial order for each of them
                                 % don't get why we subtract 1 from kappa ???
            poly = [poly kron(vX2(:,ik), vX1(:,iz))]; % complete set of polynomials basis that will be used for Chebyshev approximation
        end
    end

end

end
 
