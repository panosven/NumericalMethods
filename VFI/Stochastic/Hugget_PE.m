% Hugget Economy in Partial Equilibrium
% Solution Method: VFI with Schumaker Spline
% Panagiotis Veneris, U of Liverpool
% 28/1/2021 

clear all; clc;
close all; 

% Denote global variables

global beta n_y a_init y_init pr pp r sigma 


sigma   = 2;          % curvature of utility (inverse of EIS)
rho     = 0.06;       % discount rate
beta    = 1/(1+rho);  % discount factor
r       = 0.04;       % interest rate
sigma_y = 0.3;        % idiosyncratic earnings volatility
y1    = 1-sigma_y;    % low income
y2    = 1+sigma_y;    % high income
delta = 0.95;         % persistence of earning process
alpha = 3;            % scaling parameter for asset grid
tol   = 1e-6;         % tolerance level
maxit = 1000;         % max # of iterations
dif   = 0.1;
it    = 0;            % iteration counter


% Grid for income (y)
n_y = 2;                     % # gridpoints for income(exog.state variable)
endow = [y1,y2];             % grid for income, 1x2


% Earnings(endowment) follows a Markov process (2-state Markov)
% Initialize probability transition matrix P_ij (for stoh. income)
pr = zeros(n_y,n_y); 

% rows = future state, columns = current state
pr(1,1) = (1+delta)/2;         % transition from state 1(low income) to state 1(low income)
pr(2,1) = (1-delta)/2;         % transition from state 1(low income) to state 2(high income)
pr(1,2) = (1-delta)/2;         % transition from state 2(high income) to state 1(low income)
pr(2,2) = (1+delta)/2;         % transition form state 2(high income) to state 2(high income)

pr = [pr(1,1) pr(2,1); pr(1,2) pr(2,2)];  % prob. transition matrix, 2x2


% Grid for assets 
n_a = 200;         % # gridpoints for assets(endog.state variable)
amin = -1;         % lower bound on assets (min assets)
amax = 10;         % upper bound on assets (max assets)

agrid = (amin+(amax-amin)*(((1:n_a)/(n_a-1)).^(alpha)))';  % asset grid (denser for low values of a)


% Guess for Value Function
V = zeros(n_a,n_y);       % V_0(k,y) is a 200x2 matrix

% Initialize auxiliary matrices
V1 = zeros(n_a,n_y);
A11 = zeros(n_a,n_y);
a_init = zeros(n_a,n_y);
y_init = zeros(n_a,n_y);


pp = cell(n_y,1);

tic
while dif>tol && it<maxit
    
    it = it+1;
  
    for iy = 1:n_y
        pp{iy} = schumaker(V(:,iy), [], agrid);                 % compute schumaker interpolation coefficients for each y
    end
  
    for iy=1:n_y                                     % loop over income grid                                                  
        for ia=1:n_a                                 % loop over asset grid
            a_init = agrid(ia);                      % for each i in the state space get the initial value of a_t (ia=1 a1, ia=2 a2 etc) 
            y_init = endow(iy);                      % for each j in the state space get the initial value of y_t (iy=1 y1, iy=2 y2)
            A1 = fminbnd(@value_hugget,amin,amax);   % for each value of a_t find a' that max the Bellman
            V1(ia,iy) = -value_hugget(A1);           % collect optimized values of V in a matrix
            A11(ia,iy) = A1;                         % find asset policy function (a' policy function)
        end
    end
    
    dif = norm(V1-V);
    V=V1;
    fprintf('iteration = %d, dif = %e\n', it, dif);
  
end
toc

% Plots
figure
 hold on;
 plot(agrid,V1(:,1),'-.r','Linewidth',2);
 plot(agrid,V1(:,n_y),'b','Linewidth',2);
 hold off;
 xlabel('a');
 ylabel('V(a,y)');
 legend('V(a,y_{L})','V(a,y_{H})');
 title('Value Funtion');
 grid on;

figure
 hold on;
 plot(agrid,A11(:,1),'-.r','Linewidth',2);
 plot(agrid,A11(:,n_y),'b','Linewidth',2)
 hold off;
 xlabel('a');
 ylabel('g(a,y)');
 legend('g(a,y_{L})','g(a,y_{H})');
 title('Asset Policy Funtion');
 grid on;


function valhug = value_hugget(aprime)

% pp = schumaker(y,[],x) : calculates Schumaker's interpolation coefficients (pp)
% x :  is the vector of gridpoints 
% y :  is the vector of values of the value-function

% yq = ppval(pp,xq) 
% xq : vector of gridpoints at which i want to evaluate the values of value function
% yq : vector of interpolated values
% ppval : returns the values at the points xq of the piecewise polynomial contaiined in pp
 
 global beta n_y a_init y_init pr pp r sigma       
  
  for iy = 1:n_y                       % loop over income grid
    
    g(:,iy) = ppval(pp{iy},aprime);    % mind that inside the parentheses we need pp WITH curly brackets
                                       % returns the values of the polynomial pp evaluated at a'
  end
  
    c = y_init + (1+r)*a_init - aprime; % budget constraint(solved wrt c)
    aprime(aprime<-1)=0;
  
     if c<0
         valhug = -999999999999999999-900*abs(c); % keeps it from going negative
     else
         
        for iy = 1:n_y
           valhug=  ((c.^(1-sigma))/(1-sigma))+ beta*(pr(iy,:)*g');
        end
        
     end
         valhug = -valhug;  % make it negative since we're maximizing and code is to minimize.
    
  end

