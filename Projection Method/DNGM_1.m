% Deterministic NGM with VFI + Chebyshev Approximation (collocation) of V(k)

% Panagiotis Veneris, U of Liverpool
% 07/11/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some Theory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chebyshev approximation is part of what we call "Spectral methods", as it uses spectral bases, 
% such as orthogonal (Chebyshev) polynomials, to approximate the function of interest.

% Spectral (polynomial) methods use the same polynomial over the whole
% domain of x (the function we want to approximate is f(x)).

% In contrast, "Finite element methods" split the domain of x into non-intersecting subdomains
% (elements) and fit a different polynomial for each element (i.e fit a different approximating
% function). For example, one region can be where a borrowing constraint binds and another region 
% where the borrowing constraint does not bind. So, we can simply fit two polynomials, one in 
% each region.


% Note: With full depreciation (delta = 1) and log utility we can get
% analytical solutions for the policy and value functions, and compare them
% with the approximated to eyeball the accuracy of our numerical algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
clear;
close all;
clc;

% addpath('C:\Users\pvene\Dropbox\Codes_Library\Numerical Methods\Chebyshev Collocation\Chebyshev_functions')

addpath('Chebyshev_functions')
global param k_init theta0 


% Parameters
p.A     = 10;      % productivity
p.alpha = 0.3;     % capital share
p.beta  = 0.7;     % discount factor
p.delta = 1;       % capital depreciation rate

tol     = 1e-6;    % tolerance for VFI convergence
maxit   = 1000;    % max number of iterations
dif     = 0.1; 
it      = 0;       % iterations counter


% Grid for capital
p.kmin  = 0.1;   % min level of capital
p.kmax  = 15;    % max level of capital
p.Ngrid = 30;    % # gridpoints for capital (if we want to increse the # gridpoints, make sure we also increse the order of the polynomial: ex: p.Ngrid = 40, n=39)
p.kgrid = linspace(p.kmin,p.kmax,p.Ngrid)';


%-----------------------------
% Chebyshev Approximation
%-----------------------------

% Step 1: Set the order of the polynomial approximation (n)
n = 29;    % order of Chebyshev polynomial


% Step 2: Pick the number of nodes m (i.e # of gridpoints) at which we will approximate V(k) + compute the nodes x [-1,1]
% Note: if m = n+1 is called Chebyshev collocation method or Chebyshev interpolation. For m>n+1 we have Chebyshev regression.
m   = n+1;                                % # Chebyshev nodes
x_i = chebyshev_nodes(n);                 % (m x 1) vector of Chebyshev nodes in [-1,1], at which we will approximate V(k),
                                          % i.e nodes (gridpoints) at which
                                          % we will evaluate each Chebyshev polynomial {T_j}_j=0,...,n
                                          % x_i = x(k) in lecture notes                                      

% Since x_i's are the roots of the n-th order Chebyshev polynomial, x_i's
% are called "collocation points". In other words, x_i's are the values for
% which the basis functions T_j are zero. 

% For the orthogonal Chebyshev polynomials, the Chebyshev nodes x_i are
% also the "Quadrature nodes".

param = [p.A p.alpha p.beta p.delta p.kmin p.kmax m n]';



% Step 3: Evaluate each Chebyshev polynomial j:0,...n at all Chebyshev nodes {x_i}_i=1,..,m in [-1,1]. This gives T(x).
m   = length(x_i);
T_x = zeros(m,n+1);             % m x (n+1);
                                % row i: nodes of x_i:x_1,....,x_m, 
                                % column j: is the j-1 polynomial evaluated at each {x_i}_i=1,..,m
T_x = compute_chebypol(x_i,n);  % m x (n+1) vector T_x in [-1,1]
                                % T_x(i,j) = T_x(k) in lecture notes
                                
% The set of orthogonal Chebyshev polynomials T are called basis functions,
% as they form an orthogonal basis wrt the weight function w = (1-x^2)^(-1/2)
% Remember: Chebyshev polynomials arise when we orthogonalize ( a la
% Gram-Schmidt) monomials 1,x^1,x^2... wrt a weight function w = (1-x^2)^(-1/2) 
% in [-1,1]


% Step 4: Mappings (Option 2 in class notes)

% 4.A) map each x_i to space [p.kmin,p.kmax]
% We project each collocation point {x_i}_i=1,..,m in [-1,1] to the corresponding value {k_i}_i=1,..,m in [p.kmin,p.kmax]

k_i = transform_x_to_k(x_i,m,p); % map nodes {x_i}_i=1...m in [-1,1] to gridpoints over capital {k_i}_i=1...m in [p.kmin,p.kmax]
                                 % k_i = k(x) in lecture notes

% 4.B) map each k_i to space [-1,1]
% We project each gridpoint {k_i}_i=1,..,m in [p.kmin,p.kmax] to the corresponding value {x_i}_i=1,..,m in [-1,1]

%x_i1 = transform_k_to_x(k_i,m,p); % coincides with x_i defined above


% Auxiliary matrices (pre-allocation)
k_init  = zeros(n+1,1);
theta   = zeros(n+1,1);  % initialize matrix for chebyshev coefficients
K11     = zeros(n+1,1);  % initialize matrix for kprime (i.e savings policy today)
V_tilde = zeros(n+1,1);  % initialize matrix for value function


% Step 5: Initial guesses 
theta0 = 0.1*ones(n+1,1);  % initial guess for the Chebyshev coefficients (coefficients that weight the various polynomials)
                           % (n+1) x 1 vector of chebyshev coefficients
                        
V_tilde0 = zeros(m,1);  % initial guess for the approximated value function at each collocation point x_i
                        % m x 1 vector; 

k_length = length(k_i);

tic
while dif>tol && it<maxit         
    it = it+1;

    for h = 1:k_length         % loop over grid for k_i's
        k_init = k_i(h,1);     % for each h, get the initial value of k_i (i.e h=1 get k_1,..., h=m get k_m)
        kn_min = param(5);
        kn_max = min( param(6), param(1)*(k_init^param(2)) + (1-param(4))*k_init );
        K1     = fminbnd(@value_cheb,kn_min,kn_max); % for each value of k_i (k_1,..,k_m) find k' that max the Bellman equation  
        V_tilde(h,1) = -value_cheb(K1); % collect the optimized value into a new value function called V_tilde
        K11(h,1) = K1;         % savings policy function
    end
        
    theta(1)     = (1/m) * sum(V_tilde);
    theta(2:end) = (2/m) * (T_x(:,2:end)' * V_tilde);
            
    % check approximation
    dif1 = max(abs(theta - theta0));
    dif2 = max(abs(V_tilde - V_tilde0));
    dif  = max([dif1,dif2]);
    fprintf('Iteration: %d\tConv. tol_pol.: %g\n',it,dif); 
       
    % update
    theta0   = theta;
    V_tilde0 = V_tilde;  
    
end
toc

% Actual (true) policy and value functions (formulas are valid only for full depreciation, i.e delta = 1, and log utility).
P = p.alpha/(1-p.alpha*p.beta);
Q = (log(p.A)+(1-p.alpha*p.beta)*log(1-p.alpha*p.beta)+...
    p.alpha*p.beta*log(p.alpha*p.beta))/((1-p.beta)*(1-p.alpha*p.beta));

V_actual = P.*log(k_i)+Q*ones(size(k_i));        % actual (true) VF
kprime_actual = p.alpha*p.beta*p.A*k_i.^p.alpha; % actual (true) savings policy function



% Euler Equation Errors
e0 = (K11 - kprime_actual)./kprime_actual;% percentage error between the true and the approximated policy functions

% Plots
figure
subplot(2,1,1)
hold on
plot(k_i,V_tilde,'bo','Linewidth',2);         % chebyshev
plot(k_i, V_actual,'r:', 'Linewidth',2);  % true
hold off
xlabel('k');
ylabel('V(k)');
legend('Cheb.','True','Location','best');
title('Value Funtion');
grid on;
set(gca,'FontSize',12,'fontname','times')


%figure;
subplot(2,1,2)
hold on
plot(k_i, K11, 'bo', 'LineWidth', 2);            % chebyshev
plot(k_i, kprime_actual,'r:','Linewidth',2); % true
plot(k_i, k_i, 'k:', 'LineWidth', 2) % 45 degree line
hold off
xlabel('k')
ylabel('g(k)')
legend1 = legend('Cheb.','True','$45^{\circ}$ line','Location', 'best');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Savings Policy')
ylim([0 10])
grid on
      

figure
hold on
plot(k_i, e0, 'bs-', 'LineWidth', 2)
hold off
grid on
box on
title('Euler Equation Errors')
xlabel('k')
ylabel('percent error')
legend('Cheb.', 'Location','best')


function valuef = value_cheb(kprime)

global param k_init theta0

% Recover parameters
Alpha  = param(1);
alppha = param(2);
betta  = param(3);
dellta = param(4);
kminn  = param(5);
kmaxx  = param(6);
mm     = param(7);
nn     = param(8);

% Consumption (returns a scalar)
c = Alpha*(k_init^alppha) + (1-dellta)*k_init - kprime; % consumption function (1st argument of RHS in Bellman)

x = cos(pi/(2*mm)) * (2* ((kprime-kminn)./(kmaxx - kminn)) -1  );  % x: function that maps tomorrow's capital gridpoints (k') to the domain [-1,1]
T = compute_chebypol(x,nn); 
  
valuef = - ( log(c) + betta * T*theta0 ); % approximation of value function V(k)_tilde
end
