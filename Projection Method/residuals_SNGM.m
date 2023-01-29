function res = residuals_SNGM(theta,param)

%=====================================================================
% Find vector of residuals
% Inputs:
%          - theta: vector of initial coefficients for polynomial
%          - param: vector of parameters
% Output:  - res

% This function is a replica of Collard's code for the 
% stochastic growth model solve with the Collocation Method and 
% Gauss-Hermite quadrature to approximate the integral
%=====================================================================

global param 
global kpolynomial_nodesnum kpolynomial_order k_i T_x_k
global zpolynomial_nodesnum zpolynomial_order z_i T_x_z
global gh_weights gh_nodes

% Recover parameters
alppha = param(1);  % capital share
betta  = param(2);  % discount factor
dellta = param(3);  % capital depreciation rate
kminn  = param(4);  % lower bound for capital
kmaxx  = param(5);  % upper bound for capital
zminn  = param(6);  % lower bound for productivity
zmaxx  = param(7);  % upper bound for productivity
rho    = param(8);  % persistence of stochastic process
mu_z   = param(10); % mean of stochastic process
siggma = param(11); % inverse EIS

% Initialize auxiliary matrices
euler_rhs = []; % rhs of Euler equation
euler_lhs = []; % lhs of Euler equation


for ik = 1:kpolynomial_nodesnum     % loop over today's capital stock (capital nodes)
    for iz = 1:zpolynomial_nodesnum % loop over today's productivity (productivity nodes)
        
        % get polynomial combinations, complete polynomials basis
        T_x_kz = polynomial_comb(T_x_z(iz,:),T_x_k(ik,:)); 
        
        % approximate today's consumption function c_t = c(k_t,z_t) with Chebyshev polynomials
        c = exp(T_x_kz * theta);  % take exponent to avoid negative values on consumption
        
        % The idea is to approximate today's consumption function c_t = c(k_t,z_t) with
        % Chebyshev polynomials, and back out tomorrow's capital k_t+1 = k(k_t,z_t) using the
        % aggregate resource constraint.
              
        % back out k_t+1 from aggregate resource constraint
        kprime          = exp(z_i(iz)) * k_i(ik).^(alppha) + (1-dellta)*k_i(ik) - c; % k_t+1
        kprime_original =  (2*(log(kprime)-kminn))/(kmaxx-kminn) - 1;  % map log(k_t+1) in [p.kmin,p.kmax] to the domain [-1,1]
        T_x_kprime      = compute_chebypol(kprime_original,kpolynomial_order); % chebyshev polynomials for log(k_t+1)        
        
        % compute z_t+1 using GH quadrature nodes and weights
        zprime          = rho * z_i(iz) + (1 - rho) * mu_z + gh_nodes; % Z_t+1, using GH-grid
        zprime_original = (2*(zprime-zminn))/(zmaxx-zminn) - 1; % map productivity gridpoints (z') to the domain [-1,1]
        T_x_zprime      = compute_chebypol(zprime_original,zpolynomial_order); % chebyshev polynomials for z_t+1
             
        % compute complete chebyshev polynomials basis matrix for k',z'
        T_x_kprime_zprime  = polynomial_comb(T_x_zprime, T_x_kprime); 
             
        % compute next period's consumption
        cprime = exp(T_x_kprime_zprime * theta); % again, exponenent to get non-negative consumption
              
        % approximate the expectation using Gauss-Hermite quadrature (integration)
        integral = gh_weights'*(betta*(alppha*exp(zprime)*kprime.^(alppha-1) + 1-dellta).* cprime.^(-siggma)) /sqrt(pi);
        
        % reads as: E [ beta * (c')^(-sigma) * (a * e^(z') * (k')^(a-1) +
        % (1-delta) ] but the expectation (integral) is approximated by a
        % Gauss-Hermite polynomial, hence we multiply (*) by gh_weights and divide (/) by sqrt(pi)
       
        % compute both sides of Euler equation
        euler_rhs = [euler_rhs; log(integral)];
        euler_lhs = [euler_lhs; log(c.^(-siggma))];

    end
end


% Compute residuals and transform into vector

res = (euler_lhs - euler_rhs); % residuals
res = res(:);                  % vector of residuals

% Defining the approximation error as e(k_t,z_t) = - (1/c_t) + E[(beta*alpha*z_t+1*(k_t+1)^(alpha-1)+(1-delta))/c_t+1],
% at the true solutions c(k_t,z_t) and k(k_t,z_t) it must be that
% e(k_t,z_t) = 0 at each gridpoint k_t = k_i, z_t = z_i

% We seek to solve for the coefficients ('thetas') of the Chebyshev polynomial for which the 
% approximation error is zero. Since k_i and z_i are knowns, the only unknown are the Chebyshev
% coefficients (the 'thetas'). Hence we need to solve for the 'thetas' such that the error 
% terms are as small as possible

end
