% Deterministic NGM with VFI + Grid search
% Panagiotis Veneris, U of Liverpool
% 18/11/2020

% NOTE 1: This is the deterministic NGM, i.e there is no aggregate (and/or)
% idiosyncratic uncertainty

% NOTE 2: Grid search restricts k' to be on the grid of k (kgrid), which is
% useful since we only know values of V(k') for k' \in kgrid

clear; 
close all; 
clc;

% Define parameters
alpha   = 0.3;
beta    = 0.7;
A       = 10;
tol_pol = 1e-6;
maxit   = 1000; %max # of iterations
dif     = 0.1;
it      = 0 ; % iteration counter

% Solve for the steady state 
k_ss = (alpha*beta*A)^(1/(1-alpha));
c_ss = ((alpha*beta*A)^(alpha/(1-alpha)))*A*(1-alpha*beta);

% Grid for state variable (k)
kGridnum = 30; 
kmin     = 0.25*k_ss;
kmax     = 1.75*k_ss;
kGrid    = linspace(kmin,kmax,kGridnum); % equally-spaced grid for k 30x1 matrix

% Output at each gridpoint for k
Y  = A.*kGrid.^alpha;                         % 30x1 matrix (1st column:Y(k1), 2nd column:Y(k2)....etc)
YY = ones(kGridnum,1)*Y;                     % 30x30 matrix

KK = ones(kGridnum,1)*kGrid;                 % 30x30 matrix (1st column:K_1, 2nd column K_2

% Consumption at each level(gridpoint) of k'
C = YY-KK';

% Calculate Utility
U = log(C);

% Initial guess for value function
V  = zeros(kGridnum,1);                     % V_0(k) 30x1 matrix of 0s (guess matrix for value function): one guess for each possible state (30x1)
VV = V*ones(1,kGridnum);                    % 30x30 guess matrix for value function (1st column:V(k1), 2nd column:V(k2).....etc)
W  = U + beta*VV;                           % return function 30x30 matrix
R  = W';                                    % 30x30 matrix
V  = max(W)';                               % V(k) true value function

tic
while dif>tol_pol && it<maxit
    it = it+1; %and keep iterating
    
    VV  = V*ones(1,kGridnum);
    W   = U + beta*VV;             % W=u(c)+βV(k') --> RHS of Bellman
    V1  = max(W)';                 % solve for V1 by searching for the k' \to K that max the objective function W -->finds the max of the LHS of Bellman (V)
    dif = max(abs(V1-V));          % if dif>tolerance then update the value function (line 80)
    V   = V1;                      % update the value function (the new (updated) guess is V1)
   
     %fprintf('Iteration: %d\tConv. tol_pol.: %g\n',it,dif); 
     fprintf('iteration = %d, dif = %e\n', it, dif);
end
toc

% Find k' that gives the max V for each k
% For each value of k find the k' that max the Bellman equation (V) given the
% initial guess of the value function

[val,ind]=max(W);                      % vector of incides where W takes the max value

% Use indices to find the corresponding value of k' that gives the max value
% of V(k)
k_prime = kGrid(ind);                  % k' : capital policy function

% Consumption and output policy functions
c = Y - k_prime;      % consumption at each gridpoint of k (c = F(k) - k')
y = A.*kGrid.^alpha;  % same as y_actual since k_t is predetermined; 


%========  Actual (true) Policy Functions =======================

% True capital policy function
kprime_actual =  alpha*beta*A*kGrid.^alpha; % true (actual) k'

% True consumption policy function
c_actual = Y - kprime_actual;

% Actual output
y_actual = A.*kGrid.^alpha;  % true y (same as above since k_t is predetermined)

% Actual(true) Value Function
P = log((kGrid.^alpha).*A.*(1-alpha.*beta));
Q = (alpha./(1-alpha.*beta)).*log(alpha.*beta.*A.*(kGrid.^alpha));
Z = (log(1-alpha.*beta)+log(A)+(((beta.*alpha)./(1-alpha.*beta)).*log(A.*alpha.*beta)))./(1-beta);

v_actual = P + beta.*(Q+Z);

%================================================================


% Euler Equation Errors
e = (k_prime - kprime_actual)./kprime_actual; % percentage error between the true and the approximated policy functions

% Plots
figure('Name','VFI with Grid Search')
subplot(2,2,1)
 hold on;
 plot(kGrid,V1,'-c','Linewidth',2);
 plot(kGrid,v_actual,'b*-');
 hold off;
 xlabel('k');
 ylabel('V(k)');
 legend('VFI','True','Location','best');
 title('Value Funtion');
 box on;
 grid on;

subplot(2,2,2)
 hold on;
 plot(kGrid,k_prime,'-c','Linewidth',2);
 plot(kGrid,kprime_actual,'b*-')
 hold off;
 xlabel('k');
 ylabel('g(k)');
 legend('VFI','True','Location','best');
 title('Capital Policy Funtion');
 box on;
 grid on;

subplot(2,2,3)
 hold on;
 plot(kGrid,c,'-c','Linewidth',2);
 plot(kGrid,c_actual,'b*-');
 xlabel('k');
 ylabel('c');
 legend('VFI','True','Location','best');
 title('Consumption Policy Function');
 grid on;
 
 subplot(2,2,4)
 hold on;
 plot(kGrid,y,'-c','Linewidth',2);
 plot(kGrid,y_actual,'b*-');
 xlabel('k');
 ylabel('y');
 legend('VFI','True','Location','best');
 title('Output');
 box on;
 grid on;
 
figure('Name','Euler Errors with Grid Search')
 hold on
 plot(kGrid, e, 'bs-', 'LineWidth', 2)
 hold off
 grid on
 box on
 title('Euler Equation Errors')
 xlabel('k')
 ylabel('percent error')
 legend('Grid Search VFI', 'Location','best')
