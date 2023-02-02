% Solve the Deterministic NGM with VFI and One dimensional min. algorithm + linear interp
% Panagiotis Veneris, U of Liverpool
% 18/11/2020 (a la Sims)

% NOTE 1: This is the deterministic NGM, i.e there is no aggregate (and/or)
% idiosyncratic uncertainty

% NOTE 2: Linear interpolation preserves monotonicity and concavity of value function, but value function not differentiable

clear; 
close all; 
clc;


%Denote global variables
global beta alpha kGrid k_init A V

alpha   = 0.3;
beta    = 0.7;
A       = 10;
tol_pol = 1e-6;
maxit   = 1000; %max # of iterations
dif     = 0.1;
it      = 0; % iteration counter

%Solve for the steady state 
k_ss = (alpha*beta*A)^(1/(1-alpha));     

%Grid for state variable (k)
kGridnum = 30; 
kmin     = 0.25*k_ss;
kmax     = 1.75*k_ss;
kGrid    = linspace(kmin,kmax,kGridnum)';    % equally-spaced grid for k (30x1 matrix)


%Initial guess for value function
V = zeros(kGridnum,1);                    % V_0(k) (30x1 matrix)



%Initialize auxiliary matrices
V1     = zeros(kGridnum,1);
K11    = zeros(kGridnum,1);
k_init = zeros(kGridnum,1);

 
tic
% Start loop that iterates on Bellman equation
while dif>tol_pol && it < maxit
    it = it+1;
    
 for i=1:kGridnum                         % loop over grid for capital
   k_init = kGrid(i,1);                   % for each i in the state space get the initial value of k_t (i=1 k1, i=2 k2 etc...)
   K1 = fminbnd(@valfunbnd,kmin,kmax);    % for each value of k_t (k1,k2 etc...) find k'(=K1) that max the Bellman equation
   V1(i,1) = -valfunbnd(K1);              % collect the optimized value (of the value function) into a new value function called V1(=V_1(k))
   K11(i,1) = K1;                         % get policy function for capital (optimal k'), i.e optimal savings policy
 end
 
 % Check for convergence
    dif = norm(V1-V);
    %fprintf('Iteration: %d\tConv. tol_pol.: %g\n',it,dif);  
    fprintf('iteration = %d, dif = %e\n', it, dif);
    
 % Update guess
 V = V1;
end
toc

% Consumption policy function
c = A.*kGrid.^alpha - K11;

% Find output policy function
y = A.*kGrid.^alpha;   % same as y_actual since k_t is predetermined; 


%--------- Actual(true/closed form) policy functions----------------------

% Actual(true) Value Function
P = log((kGrid.^alpha).*A.*(1-alpha.*beta));
Q = (alpha./(1-alpha.*beta)).*log(alpha.*beta.*A.*(kGrid.^alpha));
Z = (log(1-alpha.*beta)+log(A)+(((beta.*alpha)./(1-alpha.*beta)).*log(A.*alpha.*beta)))./(1-beta);
v_actual = P + beta.*(Q+Z);

% Actual(true) capital policy function (true k')
kprime_actual = alpha*beta*A*kGrid.^alpha;              % true k'

% Actual(true) Consumption Policy Function
c_actual = A.*kGrid.^alpha - kprime_actual;              % true c

% Actual output
y_actual = A.*kGrid.^alpha  ;                    % true y (same as above since k_t is predetermined)

%-------------------------------------------------------------------------

% Euler Equation Errors
e = (K11 - kprime_actual)./kprime_actual; % percentage error between the true and the approximated policy functions

% Plots
figure('Name','VFI with Linear Interpolation')
subplot(2,2,1)
 hold on;
 plot(kGrid,V1,'-c','Linewidth',2);
 plot(kGrid,v_actual,'b*-');
 hold off;
 xlabel('k');
 ylabel('V(k)');
 legend('VFI','True','Location','best');
 title('Value Funtion');
 grid on;

subplot(2,2,2)
 hold on;
 plot(kGrid,K11,'-c','Linewidth',2);
 plot(kGrid,kprime_actual,'b*-')
 hold off;
 xlabel('k');
 ylabel('g(k)');
 legend('VFI','True','Location','best');
 title('Capital Policy Funtion');
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
 grid on;

figure('Name','Euler Errors with Linear Interpolation')
 hold on
 plot(kGrid, e, 'bs-', 'LineWidth', 2)
 hold off
 grid on
 box on
 title('Euler Equation Errors')
 xlabel('k')
 ylabel('percent error')
 legend('Linear Interp VFI', 'Location','best')

 
function val=valfunbnd(kprime)                
 
% Set up the value function to feed matlab code
 
 global beta alpha kGrid k_init A V

 g = interp1(kGrid,V,kprime,'linear');        % interpolate(evaluate) function V at kprime (k'=K1 in the main body): mapping from k to V and evaluate V at k' to find V(k')
 c = A*k_init^alpha - kprime;                 % consumption function
 if c<0
 val = -888888888888888888-800*abs(c);        % keeps it from going negative
 else
 val= log(c) + beta*g;                        % V(k) = log(c) + β * V(k')
 end
 val = -val;                                  % make it negative since we're maximizing and code is to minimize.    
                                                                                      
end
