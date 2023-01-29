% Ιncome-fluctuation model (McKay's website)
% Solution Method: Endogenous Gridpoint Method
% Panagiotis Veneris, U of Liverpool
% 1/11/2020 (Final Version 4/12/2020)

clear all;
clc;
close all;


% Define Parameters
P.gamma = 2;                % risk aversion
P.beta  = 0.98;             % discount factor
Pr.R    = 1/P.beta-1e-4;    % SS real interest rate (to have a well defined SS we need βR<1 in incomplete markets models)                 
tol     = 1e-15;            % tolerance level for policy function (if policy function>tol no convergence)

        
% Step 1: Grids for state variables (a',e)/we define our grid over assets
% tomorrow (a'=a_t+1=A) and income today (e = e_t)

%Grid for (end-of-period) assets (a'= a_t+1 = A)
P.gridsize = 200;        % number of gridpoints for a'
P.min      = 0;          % min value assets a' can take (ad-hoc borrowing limit a'>=0) 
P.max      = 20;         % max value assets a' can take
P.grid     = linspace(P.min.^0.25, P.max.^0.25, P.gridsize).^4; % equally-spaced grid for a'


% Grid for fluctuating income (e_t)
P.N     = 2;             % # of possible realizations of income
P.e1    = 0.8;           % low income  
P.e2    = 1.2;           % high income
P.endow = [P.e1,P.e2];

% Discretization of the stohastic process (fluct.income)
% Probability transition matrix P_ij (for stoh. income)
pr = zeros(P.N,P.N); 

% rows = future state, columns = current state
pr(1,1) = 0.5;     % transition from state 1 to state 1
pr(1,2) = 0.5;     % transition from state 2 to state 1
pr(2,1) = 0.5;     % transition from state 1 to state 2
pr(2,2) = 0.5;     % transition from state 2 to state 2


% Auxiliary matrices
P.tiledGrid  = repmat(P.grid,2,1); % 2x200
P.tiledEndow = repmat(P.endow',1,P.gridsize); % 2x200


% Initialize vector
G = 10+0.1*P.tiledGrid;   % guess for a_t (we give a guess for the object we are looking for (a_t).
                       
% Contraction Mapping Theorem assures that  whatever guess we give we will 
% converge to the true value, sooner of later

dif = 1;                  % criterion to check for convergence
it  = 0;                  % iteration counter

tic
% Start Loop that Iterates on Euler Equation until the fixed point is reached
while dif>tol && it<10000
    
    it = it+1;                  
    
    % Initialize c_t+1
    cp = zeros(P.N,P.gridsize); % c_t+1 2x200 matrix
    
    for i = 1:P.N
        
        % compute a_t+2
        app = interp1(G(i,:),P.grid,P.tiledGrid(i,:),'linear','extrap');  % a_t+2 = g(A,e'), 1x200
        % we have a mapping from G to P.grid (mapping from a_t to a_t+1)
        % and interpolate P.grid (a_t+1) to P.tiledGrid to we get app (a_t+2)
        % so interpolate (evaluate) P.grid at P.tiledGrid  
        
        % impose borrowing constraint
        app(app<0) = 0;
        
        % compute c_t+1
        cp(i,:) = (Pr.R*P.tiledGrid(i,:)+P.tiledEndow(i,:)-app); % c_t+1(a',e') = Ra'+e'-g(A,e'), 2x200
        
        %cp(cp<0) = 1e-6;       % to make code work for odd values of gamma                                                 
    end
   
   % compute next period's MU 
   upcp = cp.^(-P.gamma); % u'(c_t+1), 2x200
   Eupcp = [];

   for i_p = 1:P.N
      
      % compute next period's expected MU  
      Eupcp(i_p,:) = pr(:,i_p)' * upcp;   % E[u'(c_{t+1})], 2x200
                                          % (1,1): P(L|L)*c(A^L,e^L)+P(H|L)*c(A^L,e^H),exp. MU(c_t+1) at low inc. state (1st row is low inc.)
                                          % (2,1): P(L|H)*c(A^L,e^L)+P(H|H)*c(A^L,e^H),exp. MU(c_t+1) at high inc. state (2nd row is high inc.)                                       
   end
   
   % compute current period's MU
   upc = P.beta*Pr.R*Eupcp;   % u'(c_t)
   
   % compute current consumption 
   c = upc.^(-1/P.gamma);   % c_t
   
   % compute current assets  
   a = (P.tiledGrid+c-P.tiledEndow)/(Pr.R); % a_t 
   
   % Check for convergence
  if mod(it,50) == 0
    dif = max(max(abs(a-G)/(abs(a)+abs(G)+tol)));
    fprintf('it = %d, dif = %.3f\n',[it,dif]); 
  end
  
   % update guess for a_t
   G = a;  % set G = a (update your guess) and repeat the steps, otherwise stop/end loop
end
toc




% Get Policy functions (cosumption,savings)

% Consumption
figure;
subplot(1,2,1);
hold on;
plot(P.grid,c)
xlabel('Current Assets a');
ylabel('Consumption');
legend({'Low endowment','High endowment'});
legend('High/Low endowment');
grid on;

% Savings 
subplot(1,2,2);
hold on;
plot(P.grid,interp1(G(1,:),P.grid,P.grid,'linear','extrap')-P.grid) % y-axis a_t+1, x-axis a_t
plot(P.grid,interp1(G(2,:),P.grid,P.grid,'linear','extrap')-P.grid) % we can omit the ''-P.grid'' to get the usual plot (and above the same)
xlabel('Current Assets a')
ylabel('Savings a΄')
legend('Low endowment','High endowment');
grid on;


%================================================
% Panos version

% figure;
% subplot(1,2,1);
% hold on;
% plot(a(1,:),c(1,:));
% plot(a(2,:),c(2,:));
% xlabel('Current Assets a');
% ylabel('Consumption');
% legend({'Low endowment','High endowment'});
% grid on;
% 
% subplot(1,2,2);
% hold on;
% plot(a(1,:),P.grid,'r');
% plot(a(2,:),P.grid,'b');
% xlabel('Current Assets a')
% ylabel('Savings a΄')
% legend('Low endowment','High endowment');
% grid on;
%================================================

