% Solve the Deterministic NGM with VFI and One dimensional min. algorithm + cubic spline
% Panagiotis Veneris, U of Liverpool
% 18/11/2020

% Notes:
  % -1. This is the deterministic NGM, i.e there is no aggregate (and/or)
  %     idiosyncratic uncertainty

  % -2. Cubic spline ensures that the value function is differentiable, but may not preserve monotonicity and concavity of value function
  
  % -3. While cubic splines generate an approximation with smooth first and second 
  %     derivatives, it may not preserve the concavity or monotonicity of the original function.
  
clear; 
close all;
clc;

%Denote global variables
global beta alpha kGrid k_init A V

alpha   = 0.3;
beta    = 0.7;
A       = 10;
tol_pol = 1e-6;
maxit   = 1000;    % max # of iterations
dif     = 0.1;
it      = 0;       % iteration counter


%Solve for the steady state 
k_ss = (alpha*beta*A)^(1/(1-alpha));

%Grid for state variable (k)
kGridnum = 30; 
kmin     = 0.25*k_ss;
kmax     = 1.75*k_ss;
kGrid    = linspace(kmin,kmax,kGridnum)'; % equally-spaced grid for k (30x1 matrix)


%Initial guess for value function
V = zeros(kGridnum,1);  % V_0(k) (30x1 matrix)

%Initialize auxiliary matrices
V1     = zeros(kGridnum,1);
K11    = zeros(kGridnum,1);
k_init = zeros(kGridnum,1);
 
tic
while dif>tol_pol && it < maxit
  it  = it+1;
    for i=1:kGridnum                            % loop over grid for capital
        k_init = kGrid(i,1);                    % for each i in the state space get the initial value of k_t (i=1 k1, i=2 k2 etc...)
        K1 = fminbnd(@valfunspline,kmin,kmax);  % for each value of k_t (k1,k2 etc...) find k' that max the Bellman equation --> we iterate on the Bellman Equation (i.e value function)
        V1(i,1) = -valfunspline(K1);            % collect the optimized value into a new value function called V1
        K11(i,1) = K1;                          % find policy function for capital (optimal k' = kprime = k_t+1 = Κ1)
    end
 
   dif = norm(V1-V);
   %dif = max(abs(V1-V));
   V   = V1;
   %fprintf('Iteration: %d\tConv. tol_pol.: %g\n',it,dif); 
    fprintf('iteration = %d, dif = %e\n', it, dif);
end
toc

% Find consumption and out[ut policy functions
c = A.*kGrid.^alpha - K11; 
y = A.*kGrid.^alpha;   % same as y_actual since k_t is predetermined; 


%========  Actual (true) Policy Functions =======================

% Actual(true) capital policy function (true k')
kprime_actual = alpha*beta*A*kGrid.^alpha;       % true k'
 
% Actual(true) Consumption Policy Function
c_actual = A.*kGrid.^alpha - kprime_actual;      % true c

% Actual(true) Value Function
P = log((kGrid.^alpha).*A.*(1-alpha.*beta));
Q = (alpha./(1-alpha.*beta)).*log(alpha.*beta.*A.*(kGrid.^alpha));
Z = (log(1-alpha.*beta)+log(A)+(((beta.*alpha)./(1-alpha.*beta)).*log(A.*alpha.*beta)))./(1-beta);

v_actual = P + beta.*(Q+Z);  % true V

% Actual output
y_actual = A.*kGrid.^alpha;  % true y (same as above since k_t is predetermined)

%================================================================

% Euler Equation Errors
e = (K11 - kprime_actual)./kprime_actual; % percentage error between the true and the approximated policy functions

% Plots
figure('Name','VFI with Spline')
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
 plot(kGrid,K11,'-c','Linewidth',2);
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
 hold off;
 xlabel('k');
 ylabel('c');
 legend('VFI','True','Location','best');
 title('Consumption Policy Function');
 box on;
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

 figure('Name','Euler Errors with Spline')
 hold on
 plot(kGrid, e, 'bs-', 'LineWidth', 2) % euler equation errors
 hold off
 grid on
 box on
 title('Euler Equation Errors')
 xlabel('k')
 ylabel('percent error')
 legend('Spline VFI', 'Location','best')


function val=valfunspline(kprime)              
 
 global beta alpha kGrid k_init A V

 g = interp1(kGrid,V,kprime,'spline');         % interpolate(evaluate) V at kprime (k') to get g(=V_t+1).
                                               % cubic spline interpolation to evaluate the value function at points not on the kGrid
 
% We essentially approximate the value function by a cubic spline defined on a set of grid points for
% capital, but do not restrict optimal choices for capital to lie on the
% grid --> g is our approximate value function

 c = A*k_init^alpha - kprime;                  % consumption function
 if c<0
 val = -888888888888888888-800*abs(c);         % value function; keeps it from going negative
 else
 val= log(c) + beta*g;                         % value function
 end
 val = -val;                                   % value function; make it negative since we're maximizing and code is to minimize.
end
