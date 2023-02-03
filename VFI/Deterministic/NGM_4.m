% Deterministic NGM with VFI + One dimensional min. algorithm + Schumaker spline
% Panagiotis Veneris, U of Liverpool
% 18/11/2020 

% Notes:
  % -1. This is the deterministic NGM, i.e there is no aggregate (and/or)
  % idiosyncratic uncertainty

  % -2. Schumaker spline ensures that the value function is differentiable, and preserves monotonicity and 
  % concavity of value function. In that sense, it is considered as a shape
  % preserving spline(shape=concavity/convexity+monotonicity).

  % -3. Since it calculates derivatives in each iteration, it comes at a cost of lower speed


clear; 
close all; 
clc;

%addpath('VFI_functions');

% Denote global variables
global beta alpha kGrid k_init A V 

alpha   = 0.3;
beta    = 0.7;
A       = 10;
tol_pol = 1e-6;
maxit   = 1000; %max # of iterations
dif     = 0.1;
it      = 0; % iteration counter

% Solve for the steady state 
k_ss = (alpha*beta*A)^(1/1-alpha);

% Grid for state variable (k)
kGridnum = 30; 
kmin     = 0.25*k_ss;
kmax     = 1.75*k_ss;
kGrid    = linspace(kmin,kmax,kGridnum)'; % equally-spaced grid for k (30x1 matrix)


% Initial guess for value function
V = zeros(kGridnum,1); % V_0(k) (30x1 matrix)

% Initialize auxiliary matrices
V1     = zeros(kGridnum,1);
K11    = zeros(kGridnum,1);
k_init = zeros(kGridnum,1);

tic
while dif>tol_pol && it < maxit 
   it = it+1; 
   
     for i=1:kGridnum
       k_init   = kGrid(i,1);
       K1       = fminbnd(@valfun,kmin,kmax); 
       V1(i,1)  = -valfun(K1);
       K11(i,1) = K1;
     end 
   dif = norm(V1-V);
   V=V1;

   %fprintf('Iteration: %d\tConv. tol_pol.: %g\n',it,dif); 
    fprintf('iteration = %d, dif = %e\n', it, dif);
end
toc

% Consumption and output policy function
c = A.*kGrid.^alpha - K11;
y = A.*kGrid.^alpha;   % same as y_actual since k_t is predetermined; 


%--------- Actual(true/closed form) policy functions----------------------

% Actual(true) Value Function
P = log((kGrid.^alpha).*A.*(1-alpha.*beta));
Q = (alpha./(1-alpha.*beta)).*log(alpha.*beta.*A.*(kGrid.^alpha));
Z = (log(1-alpha.*beta)+log(A)+(((beta.*alpha)./(1-alpha.*beta)).*log(A.*alpha.*beta)))./(1-beta);
v_actual = P + beta.*(Q+Z);

% Actual(true) capital policy function (true k')
kprime_actual = alpha*beta*A*kGrid.^alpha; % true k'
 

% Actual(true) Consumption Policy Function
c_actual = A.*kGrid.^alpha - kprime_actual; %(true c)

% Actual(true) Output
y_actual = A.*kGrid.^alpha  ;  % true y (same as above since k_t is predetermined)

%---------------------------------------------------------------------------

% Euler Equation Errors
e = (K11 - kprime_actual)./kprime_actual; % percentage error between the true and the approximated policy functions


figure('Name','VFI with Schumaker Spline')
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

%figure
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
 
%figure
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
 
 figure('Name','Euler Errors with Schumaker Spline')
 hold on
 plot(kGrid, e, 'bs-', 'LineWidth', 2)
 hold off
 grid on
 box on
 title('Euler Equation Errors')
 xlabel('k')
 ylabel('percent error')
 legend('Schumaker VFI', 'Location','best')
  
function val=valfun(kprime)               
 
% pp = schumaker(y,[],x] : calculates Schumaker's interpolation coefficients 
% x :  is the vector of gridpoints 
% y :  is the vector of values of value-function

% yq = ppval(pp,xq) 
% xq : vector of gridpoints at which i want to evaluate the values of value function
% yq : vector of interpolated values
% ppval : returns the values at the points xq of the piecewise polynomial contaiined in pp
 
 global beta alpha kGrid k_init A V 
 
 pp = schumaker(V, [] ,kGrid);          % this calculates Schumaker's interpolation coefficients, kGrid is the vector of gridpoints
                                        % V is the vector of value-function values
 g = ppval(pp, kprime);                 % kprime is the vector of grid points at which i want to evaluate/interpolate my value-function values
 c = A*k_init^alpha - kprime;           % consumption function
 if c<0
 val = -888888888888888888-800*abs(c);  % keeps it from going negative
 else
 val= log(c) + beta*g;
 end
 val = -val;                            % make it negative since we're maximizing and code is to minimize.
end
