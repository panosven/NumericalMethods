function [x,w] = GHernodes(n)


%=======================================================================
% Panagiotis Veneris, U of Liverpool, November 2022

% Function that: computes the nodes and weights of the Gauss-Hermite
% polynomials

% Input: # of nodes to compute (n)

% Output:
% - x: GH nodes 
% - w: associated weights

% This codes is an exact translation of Numerical Recipes'
% Fortran code (3rd edition, pg 185)

%======================================================================

maxit = 30; % max number of iterations
PIM4  = 0.7511255444; 
tol   = sqrt(eps);

% Initialize matrices
x    = zeros(n,1);
w    = zeros(n,1);
m    = (n+1)/2;

% Initial guesses for the first four roots
for i = 1:m   % loop over desired roots 
    if i == 1
        z = sqrt(2*n+1) - 1.85575 * ((2*n+1)^(-0.16667)); % initial guess for the largest root
    elseif i == 2
        z = z - 1.14 * (n^0.426) / z; % initial guess for the second largest root
    elseif i == 3
        z = 1.86 * z - 0.86 * x(1);   % initial guess for the third largest root
    elseif i == 4
        z = 1.91 * z - 0.91 * x(2);   % Initial guess for the fourth largest root
    else
        z = 2 * z - x(i-2);           % initial guess for the other roots
    end
    
    
    for it = 1:maxit % refinement by Newton’s method.
        p1 = PIM4;
        p2 = 0;
        for j = 1:n % Loop up the recurrence relation to get the Hermite polynomial evaluated at z
            p3 = p2;
            p2 = p1;
            p1 = z * sqrt(2/j) * p2 - sqrt((j-1)/j) * p3;
        end
        % p1 is now the desired Hermite polynomial. We next compute pp, its derivative, by the relation (4.6.21) using p2, 
        % the polynomial of one lower order
        pp = sqrt(2*n) * p2;
        z1 = z;
        z = z1 - p1/pp; % Newton’s formula
        if abs(z1-z) < tol
            break;
        end
    end
    x(i) = z;             % store the root (roots of each Hermite polynomial)
    x(n+1-i) = -z;        % and its symmetric counterpart
    w(i)	= 2/(pp*pp);  % Compute the weight (weights of each Hermite polynomial)
	w(n+1-i)= w(i);       % and its symmetric counterpart
end
        
end
