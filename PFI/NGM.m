% Problem Set 4: Advanced Macroeconomics
% Solve the NGM with PFI and Linear Interpolation
% Panagiotis Veneris, U of Liverpool
% 18/11/2020 

% Notes: In PFI we give an initial guess about the policy function and then iterate on
% the Euler equation until we find the true policy function

clear; 
close all;
clc;

% Denote global variables to be used inside functions
global    kGrid alpha   k_prime beta A  k0 

% Define parameters
alpha   = 0.3;
beta    = 0.7;
A       = 10;
tol_pol = 1e-6;
maxit   = 1000; %max # of iterations
dif     = 0.1;
it      = 0; % iteration counter


% Grid for state variable (k
k_ss     = (beta*A*alpha)^(1/(1-alpha));  % ss level of capital
kGridnum = 30;        % # nodes for capital grid
kmin     = 0.25*k_ss; % lower bound for capital grid
kmax     = 1.75*k_ss; % upper bound for capital grid
kGrid    = linspace(kmin,kmax,kGridnum)'; % equally-spaced grid for k, 30x1 vector


% Initial guess for k' = g_0(k) 
k_prime = zeros(kGridnum,1);

x0 = [kmin 0.99*kmax]; % bounds for fzero

% Auxiliary matrices
k0    = zeros(kGridnum,1);
k_new = ones(kGridnum,1);


tic
while dif > tol_pol && it < maxit
  it=it+1;
  
    for i = 1:kGridnum                       % loop over grid for capital
         k0 =kGrid(i,1);                     % for each i in the state space get the initial value of k_t (i=1 k1, i=2 k2 etc...)
         k1 = fzero(@(kp) euler(kp),x0);     % k1 is a vector of k' found by fzero --> here we iterate on the Euler equation
         k_new(i,1) = k1;                    % collect k' in the matrix k_new, this matrix is the policy function for capital (k')
    end
    
    dif = norm(k_new-k_prime);
    k_prime = k_new;
     
    %fprintf('Iteration: %d\tConv. tol_pol.: %g\n',it,dif); 
    fprintf('iteration = %d, dif = %e\n', it, dif);   
end
toc

% Consumption and output policy functions
c     = A*kGrid.^alpha - k_new;
y_pfi = A.*kGrid.^alpha; % same as below since k_t is predetermined


%======= Actual (true) Policy Functions ============================

% Actual(true) capital policy function 
kprime_actual = alpha*beta*A*kGrid.^alpha;          % true k'

% Actual(true) consumption policy function
c_actual = A.*kGrid.^alpha - kprime_actual;         % true c

% Actual output
y_actual = A.*kGrid.^alpha  ;                  % true y (same as above since k_t is predetermined)

%===================================================================

% Euler Equation Errors
e = (k_new - kprime_actual)./kprime_actual; % percentage error between the true and the approximated policy functions
 
figure
subplot(2,2,1)
 hold on;
 plot(kGrid,k_new,'-c','Linewidth',2);
 plot(kGrid,kprime_actual,'b*-')
 hold off;
 xlabel('k');
 ylabel('g(k)');
 legend('PFI','True','Location','best');
 title('Capital Policy Funtion');
 grid on;

subplot(2,2,2)
 hold on;
 plot(kGrid,c,'-c','Linewidth',2);
 plot(kGrid,c_actual,'b*-');
 xlabel('k');
 ylabel('c');
 legend('PFI','True','Location','best');
 title('Consumption Policy Function');
 grid on;
 
subplot(2,2,3)
 hold on;
 plot(kGrid,y_pfi,'-c','Linewidth',2);
 plot(kGrid,y_actual,'b*-');
 xlabel('k');
 ylabel('y');
 legend('PFI','True','Location','best');
 title('Output');
 grid on;
 
subplot(2,2,4)
 hold on
 plot(kGrid, e, 'bs-', 'LineWidth', 2)
 hold off
 title('Euler Equation Errors')
 xlabel('k')
 ylabel('percent error')
 legend('PFI', 'Location','best')
 grid on

 
 function y = euler(kp) % this kp is k' found by fzero

 global beta alpha kGrid A k_prime k0 
 c = A*k0.^alpha - kp; % consumption function, k is k_t+1 = k'

 uprime = c.^(-1);
 k_doublep = interp1(kGrid,k_prime,kp,'linear'); % interpolate guess (k_prime) at true k' to get k''
 
 % Euler equation
 y= uprime - beta*(A*kp.^alpha - k_doublep).^(-1)*(A*alpha*(kp.^(alpha-1))); % iteratre upon the Euler equation
 
 end
