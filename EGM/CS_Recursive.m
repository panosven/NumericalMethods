% Ιncome-fluctuation problem with CRRA utility in recursive formulation
% Solution Method: EGM
% Panagiotis Veneris, U of Liverpool
% 20/2/2021 

clear all;
close all;
clc;

% Define Parameters
P.gamma = 2;             % risk aversion
P.beta  = 0.98;          % discount factor
Pr.R    = 1/P.beta-1e-4; % SS real interest rate    
tol     = 1e-15;             
maxit   = 10000;
delta   = 0.95;          % persistence of earning process                 


% Step 1: Grids for state variables (a',e)/we define our grid over assets
% tomorrow (a'=a_t+1=A)

%Grid for (end-of-period) Assets (a'= a_t+1 = A)
P.gridsize = 200;        % number of gridpoints for a'
P.min      = 0;          % min value assets a' can take (ad-hoc borrowing limit a'>=0)
P.max      = 20;         % max value assets a' can take
P.grid     = linspace(P.min,P.max,P.gridsize); % grid for a' (a_t+1)


% Grid for fluctuating income (e)
P.N     = 2;                % # of possible realizations of income
P.e1    = 0.8;              % low income  
P.e2    = 1.2;              % high income
P.endow = [P.e1,P.e2];

pr = zeros(P.N,P.N);    % Probability transition matrix P_ij (for stoh. income)

% rows = future state, columns = current state

pr(1,1) = 0.5;          % transition from state 1 to state 1
pr(2,1) = 0.5;          % transition from state 1 to state 2 
pr(1,2) = 0.5;          % transition from state 2 to state 1
pr(2,2) = 0.5;          % transition from state 2 to state 2

% Alternative way for transition matrix

%pr(1,1) = (1+delta)/2;            % transition from state 1 (low) to state 1 (low)
%pr(2,1) = (1-delta)/2;            % transition from state 1 (low) to state 2 (high)
%pr(1,2) = (1-delta)/2;            % transition from state 2 (high) to state 1 (low)
%pr(2,2) = (1+delta)/2;            % transition from state 2 (high) to state 2 (high)

% Auxiliary matrices
P.tiledGrid = repmat(P.grid,2,1); % 2x200
P.tiledEndow = repmat(P.endow',1,P.gridsize); % 2x200


% Guess for a
G   = 10+0.1*P.tiledGrid;   % guess for a_t (we give a guess for the object we are looking for (a_t, V_t in VFI etc).                                                


dif = 1;                    % criterion to check for convergence
it  = 0;                    % iteration counter
disp('solving for policy rules')

% Start Loop that Iterates on Euler Equation until the fixed point is reached
while dif>tol && it<maxit
    
    it = it+1;                  
    
    % Create c_t+1
    cp = zeros(P.N,P.gridsize);   % c_t+1, 2x200 matrix (cp=cprime=c'=c_t+1)
    
    for i=1:P.N
        app = interp1(G(i,:),P.grid,P.tiledGrid(i,:),'linear','extrap'); % a_t+2=g(A,e'), 1x200     
        app(app<0) = 0;           
        cp(i,:) = Pr.R*P.tiledGrid(i,:)+P.tiledEndow(i,:)-app; % c_t+1(a',e') = Ra'+e'-g(A,e'), 2x200
        %cp(cp<0) = 1e-6;       % to make code work for odd values of P.gamma  
    end
    
   upcp  = cp.^(-P.gamma); % u'(c_t+1), 2x200
   Eupcp = [];
   
   for i_p = 1:P.N 
      Eupcp(i_p,:) = pr(:,i_p)' * upcp;  
   end

   upc = P.beta*Pr.R*Eupcp; % u'(c_t)
   c   = upc.^(-1/P.gamma);   % c_t
   a   = (P.tiledGrid+c-P.tiledEndow)/Pr.R; % a_t 
   
   % Check for convergence
   
  if mod(it,50) == 0
    dif = max(max(abs(a-G)/(abs(a)+abs(G)+tol)));
    fprintf('it = %d, dif = %.3f\n',[it,dif]); 
  end
  
   G = a;  % set G=a (update your guess) and repeat the steps, otherwise stop/end loop

end


% Find Value Function (since it is auxilliary we can find it separately
% from the rest of the model)

ap(1,:) = interp1(G(1,:),P.grid,P.grid,'linear','extrap'); % a_t+1, e^L
ap(2,:) = interp1(G(2,:),P.grid,P.grid,'linear','extrap'); % a_t+1, e^H
ap(ap<P.min) = P.min; 

c_pr(1,:) = interp1(G(1,:),c(1,:),P.grid,'linear','extrap'); % c_t+1, e^L,(c_pr=cprime=c')
c_pr(2,:) = interp1(G(2,:),c(2,:),P.grid,'linear','extrap'); % c_t+1, e^H
u = (c_pr.^(1-P.gamma)-1)./(1-P.gamma); % utility
%u = (c.^(1-P.gamma)-1)./(1-P.gamma); % utility (same results as above line)

% Guess for V
V = 0*P.tiledGrid;    % Guess for Value Function (V_t)

dif = 1;
it  = 0;

disp('solving for value function')
while dif>tol && it<maxit
    
   it = it+1;
    
    for i = 1:P.N
        Vp(i,:) = interp1(P.grid,V(i,:),ap(i,:),'linear','extrap'); % V_t+1
    end
    
    E_Vp = []; % Initialize E_Vp
    
    for i = 1:P.N      
        E_Vp(i,:) = pr(:,i)'*Vp;     % E_t V_t+1
    end
    
    V_new = u + P.beta*E_Vp;         % V_t
     
     if mod(it,100) == 0
         dif = max(max(abs(V-V_new)/(abs(V)+abs(V_new)+tol)));
         fprintf('iter = %d, dif = %.12f\n',[it,dif])
     end
     
    V = V_new;
    
end


% Get Policy functions (cosumption,savings)

% Consumption
figure;
subplot(1,2,1);
hold on;
plot(P.grid,c_pr)
xlabel('Current Assets a');
ylabel('Consumption');
legend({'Low endowment','High endowment'});
%legend('High/Low endowment');
grid on;

% Savings 
subplot(1,2,2);
hold on;
plot(P.grid,ap-P.grid) % y-axis a_t+1, x-axis a_t
xlabel('Current Assets a')
ylabel('Savings a΄')
legend('Low endowment','High endowment');
grid on;
       
% Value Function
figure;
plot(P.grid,V_new);
xlabel('Current Assets a');
ylabel('V_{t}');
legend('Low endowment','High endowment');
title('Value Function');
grid on;

save CRRApar.mat V_new  % save V_new to use it as a guess in EZ_CS case


%==============================================================================
% An increase in the persistence (alternative trans.function) of the shock 
%(delta = 0.95 vs delta=0.90) makes high income HHs save less consume more, 
% and low income HHs save more and consume less

% An increase in RA (RA=6) makes both high and low income HHs save more
% (more precautionary savings) and consume less
%===============================================================================
