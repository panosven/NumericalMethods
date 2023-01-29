% Income-Fluctuation with Epstein-Zin preferences
% Solution Method: Endogenous Gridpoint Method (Carroll 2005) 
% Panagiotis Veneris, U of Liverpool
% 15/1/2021

clear all;
clc;
close all;

addpath('EGM_functions')

%====================
% Parameters
%====================

P.gamma = 2;             
P.beta  = 0.98;                  % discount factor
Pr.R    = 1/P.beta-1e-4;         % SS real interest rate   
P.N     = 2;                     % number of possible realizations of income
tol     = 1e-15;                 % tolerance level for policy function 
P.RA    = 5;                     % risk aversion parameter
P.alpha = 1-P.RA/P.gamma;        % a=0 (P.RA=P.gamma) reduces to CRRA case
maxit   = 50000;                 % maximum number of iterations


% Grid for (end-of-period) Assets (a'= a_t+1 = A)
P.gridsize = 200;                                                      % # of gridpoints for a'
P.min      = 0;                                                        % min value assets a' can take (ad-hoc borrowing limit a'>=0)
P.max      = 20;                                                       % max value assets a' can take
%P.grid    = linspace(P.min.^0.25, P.max.^0.25, P.gridsize).^4;        % McKay's way (doesn't work in this code)
P.grid     = linspace(P.min,P.max,P.gridsize);                         % equally-spaced grid for a'
%P.grid    = P.min + ((1:P.gridsize)/P.gridsize).^2 * (P.max-P.min);   % Lorenzoni's way (assign more gridpoints where the curvature is higher), works better for 1e-14

% Grid for fluctuating income (e)
P.N     = 2;                % # gridpoints for income
P.e1    = 0.8;              % low income  
P.e2    = 1.2;              % high income
P.endow = [P.e1,P.e2];


pr = zeros(P.N,P.N);      % Probability transition matrix P_ij (for stoh. income)

% rows = future state, columns = current state

pr(1,1) = 0.5;            % transition from state 1 to state 1
pr(2,1) = 0.5;            % transition from state 1 to state 2
pr(1,2) = 0.5;            % transition from state 2 to state 1
pr(2,2) = 0.5;            % transition from state 2 to state 2


%========================
% Initialize matrices
%========================

% Auxiliary matrices
P.tiledGrid  = repmat(P.grid,2,1);                  % 2x200
P.tiledEndow = repmat(P.endow',1,P.gridsize);       % 2x200


% Guess for a (a = a_t)             
G = 10+0.1*P.tiledGrid;  % guess for a_t: we give a guess for the object we are looking for - a_t -  
                         % guess is that a_t = 10 + 0.1 * a_t+1
                                         
% Guess for V            % guess for V_t: needed in order to find V_t+1                         
load CRRApar             % loads V_new drawn from CRRA case ( ''CS_CRRA'' code)
V1 = V_new;               

%load CRRApanosversion
%V1 = PF.V;
dif = 1;                 % criterion to check convergence
it  = 0;                 % iteration counter


disp('solving for policy rules')
tic
% Start Loop that Iterates on Euler Equation until the fixed point is reached
while dif>tol  && it<maxit
    
    it = it+1;     
    
    % Initialize matrices 
    cp  = zeros(P.N,P.gridsize);  % c_t+1 2x200 matrix, cp=c'=c_t+1
    v_v = zeros(P.N,P.gridsize);  % V_t+1 2x200 matrix
    
    for i = 1:P.N
       
        app = interp1(G(i,:),P.grid,P.tiledGrid(i,:),'linear','extrap');     % a_t+2=g(A,e'), 1x200    
        Vp  = interp1(G(i,:),V1(i,:),P.tiledGrid(i,:),'linear','extrap');    % find V_t+1(a_t+1,e_t+1), 1x200
        
        Vp(Vp>0)   = -1e-6;                                                  % restric V_t+1 to be < 0 everywhere (and then put (-) in front)
        app(app<0) = 0;                                                      % a_t+2 >= 0 from borrowing constraint
       
        cp(i,:)  = (Pr.R*P.tiledGrid(i,:)+P.tiledEndow(i,:)-app);            % c_t+1(a',e')= Ra'+e'-g(A,e'), 2x200
        cp(cp<0) = 1e-6;                                                     % c_t+1 >= 0
        v_v(i,:) = Vp;                                                       % V_t+1, 2x200    <0 everywhere                                                   
                                                         
    end
    
   % Find u'(c_t+1) 
   upcp = cp.^(-P.gamma);  % u'(c_t+1) = (c_t+1)^(-γ), 2x200 
   
   % Initialize matrices
   E_v1 = zeros(P.N,P.gridsize);        % 2x200
   E_v2 = zeros(P.N,P.gridsize);        % 2x200
     
   for i_p = 1:P.N
                                     
      E_v1(i_p,:) = pr(:,i_p)'*(((-v_v).^(-P.alpha)).*upcp);   % 2x200,      E [ (-V_t+1)^(-a) * u'(c_t+1) ], since u<=0 everywhere 
      E_v2(i_p,:) = pr(:,i_p)'*((-v_v).^(1-P.alpha));          % 2x200,      E [ (-V_t+1)^(1-a) ], >0 everywhere
      
   end
       
   % Back out  true c_t from Euler and true a_t from budget constraint
   upc = P.beta*Pr.R*((E_v2).^(P.alpha/(1-P.alpha))).*(E_v1);       % u'(c_t)= beta*R*[(E_t*(-V_t+1)^(1-alpha)]^(alpha/1-alpha) * Ε_t [(-V_t+1)^(-alpha) * u'(c_t+1)]
   c   = upc.^(-1/P.gamma);                                         % true c_t (using Euler for consumption)
   a   = (P.tiledGrid+c-P.tiledEndow)/(Pr.R);                       % true a_t (using BC, solved for a_t)
    
   % Get true V_t 
   V_update = c.^(1-P.gamma).*((1/(1-P.gamma))) - P.beta * ((E_v2).^(1/(1-P.alpha))); % true V_t 
        
  % Check for convergence
  if mod(it,50) == 0
      dif1 = max(max(abs(a-G)/(abs(a)+abs(G)+tol)));
      dif2 = max(max(abs(V_update-V1)./(abs(V_update)+abs(V1)+tol)));
      dif  = max([dif1,dif2]);
      fprintf('it = %d, dif = %.15f\n',[it,dif]);      
  end
  
  % update
   G  = a;         % set G = a (update your guess for a_t) and repeat the steps, otherwise stop/end loop
   V1 = V_update;  % update the Value function guess
  
end
toc

% Get output
y = Pr.R*a + P.tiledEndow; %% ADDITION 19/7/2021

%kk = interp1(G(1,:),P.grid,P.grid,'linear','extrap');
%===============
% Plots
%===============

% Get Policy functions (cosumption,savings)

% Consumption
figure;
subplot(2,2,1);
hold on;
plot(P.grid,c)
xlabel('Current Assets (a)');
ylabel('Consumption');
legend({'Low endowment','High endowment'},'Location','best');
title('Consumption Policy Function');
grid on;

% Savings 
subplot(2,2,2);
hold on;
plot(P.grid,interp1(G(1,:),P.grid,P.grid,'linear','extrap')-P.grid) % y-axis a_t+1, x-axis a_t
plot(P.grid,interp1(G(2,:),P.grid,P.grid,'linear','extrap')-P.grid) % we can omit the ''-P.grid'' to get the usual plot (and above the same)
xlabel('Current Assets (a)')
ylabel('Savings a(t+1)')
legend({'Low endowment','High endowment'},'Location','best');
title('Asset Policy Function')
grid on;
   
% Value Function
subplot(2,2,3);
plot(P.grid,V_update);
xlabel('Current Assets (a)');
ylabel('V(t)');
legend({'Low endowment','High endowment'},'Location','best');
title('Value Function');
grid on;

% Output
subplot(2,2,4);
plot(P.grid,y);
xlabel('Current Assets (a)');
ylabel('y');
legend({'Low endowment','High endowment'},'Location','best');
title('Output');
grid on;

%return

%==========================================================================
% Monte Carlo simulation % in MC simulation we go forward --> from b_t we
% go to b_t+1, from c_t to c_t+1 etc
%==========================================================================

% 1) Simulation Parameters
nTime = 3000;                            % length of simulation (simulation periods)
nCut  = 30;                              % burn first 300 periods to make sure our results do not depend on initial conditions
nInitialBond = P.grid((P.gridsize)/2);   % set initial state: initial value for (endogenous) state variable
                                       

% 2) Simulating a markov chain 
% Simulate to find next periods' endowment
rng('default');                                      % random number generator
[endowment,~] = fMarkov(pr,nTime+1,1,1:length(pr));  % simulation of a markov chain (endowment follows a Markov process) to find next periods' endowment
%[endowment,state] = fMarkov(pr,nTime+1,1,1:length(pr));  % low income if endowment=1, high income if endowment=2                                                   
                                                     
% 3) Initialize matrices for simulated assets, endowment and output
simBPrime = zeros(nTime,1);
simC      = zeros(nTime,1);
simY      = zeros(nTime,1);

a_today = a';    % period-t assets
cpol    = cp;    % c_t+1
                               

% 4) Run the simulation
for iTime = 1:nTime
 
   % simProd(iTime,1) = P.endow(endowment(iTime));
    
    if iTime == 1 
         %simBPrime(iTime,1) = interp1(vBGrid,mBPrime(:,vStateIndex(iTime)),nInitialBond,'linear','extrap');
         simBPrime(iTime,1) = interp1(P.grid,a_today(:,endowment(iTime)),nInitialBond,'linear','extrap');  % (simulated) path of assets (period t=1 assets)
         %simC(iTime,1)      = interp1(P.grid,c(endowment(iTime),:),nInitialBond,'linear','extrap');       % (simulated) path of consumption                       
         simC(iTime,1)      = interp1(P.grid,cpol(endowment(iTime),:),nInitialBond,'linear','extrap');
         simY(iTime,1)      = Pr.R*simBPrime(iTime)+ P.endow(:,endowment(iTime)); %% ADDITION 19/7/2021: (simulated) path of output
    else

        %simBPrime(iTime,1) = interp1(vBGrid,mBPrime(:,vStateIndex(iTime)),simBPrime(iTime-1),'linear','extrap');
        simBPrime(iTime,1) = interp1(P.grid,a_today(:,endowment(iTime)),simBPrime(iTime-1),'linear','extrap'); % period t=2 till t=3000 assets
        %simC(iTime,1)      = interp1(P.grid,c(endowment(iTime),:),simBPrime(iTime-1),'linear','extrap');
        simC(iTime,1)      = interp1(P.grid,cpol(endowment(iTime),:),simBPrime(iTime-1),'linear','extrap');
        simY(iTime,1)      = Pr.R*simBPrime(iTime-1) + P.endow(:,endowment(iTime-1)) ;  %% ADDITION 19/7/2021
    end
   
end

% Plot bond distribution
figure;
[f,xi]=ksdensity(simBPrime(:,1),'function','pdf','support',[min(simBPrime(:,1))-0.001 max(simBPrime(:,1))+0.1] );
plot(xi,f,'linewidth',2);
xlabel('Assets','FontSize',12,'fontname','times');
ylabel('\Phi(a,y)','FontSize',12,'fontname','times')
title('Bond distribution');

% xi = [min(simBPrime(:,1))-0.001 max(simBPrime(:,1))+0.1 : points (100) at which
%       we evaluate the function f

% f: probability denstity estimate of the (unobservable) density function 
%    of simBprime(bond holdings)

fprintf('\n');
fprintf('Averages \n');
%fprintf('Average b/y = %.3f \n',mean(simBPrime(nCut:end-1)));
fprintf('Average c/y = %.3f \n',mean(simC(nCut:end-1)));

figure;
subplot(3,1,1);
plot(simBPrime(:,1),'linewidth',1.5,'color','r');
set(gca,'FontSize',12,'fontname','times')
xlabel('Period','FontSize',12,'fontname','times');
ylabel('a','FontSize',12,'fontname','times')
title('assets');
grid on

subplot(3,1,2);
plot(simC(:,1),'linewidth',1.5,'color','g');
set(gca,'FontSize',12,'fontname','times')
xlabel('Period','FontSize',12,'fontname','times');
ylabel('c','FontSize',12,'fontname','times');
title('Consumption');
grid on

subplot(3,1,3);
plot(simY(:,1),'linewidth',1.5,'color','k');
set(gca,'FontSize',12,'fontname','times')
xlabel('Period','FontSize',12,'fontname','times');
ylabel('y','FontSize',12,'fontname','times');
title('Output');
grid on


% Compute long run moments
ev1     = sprintf('\n Over the simulation the means and variances were  \n');
ev2     = sprintf('Variable  \t Mean       \t Std.dev \t Rel. to GDP \t Corr. with GDP \t Min    \t Max  \t Autocorr \n');

%%% consumption 
X       = simC(nCut:nTime-1);
Xl      = simC(nCut-1:nTime-2);
Y       = simY(nCut:nTime-1);
ev3     = sprintf('C        \t %f   \t %f \t %f \t %f \t %f \t %f  \t %f \n',mean(X),std(X),std(X)/std(Y), corr(X,Y), min(X),max(X),corr(X,Xl));

% std(X)/std(Y) is the volatility of X relative to Y

%%% assets
X       = simBPrime(nCut:nTime-1);
X1      = simBPrime(nCut-1:nTime-2);
Y       = simY(nCut:nTime-1);
ev4     = sprintf('B        \t %f   \t %f \t %f \t %f \t %f \t %f  \t %f \n',mean(X),std(X),std(X)/std(Y), corr(X,Y), min(X),max(X),corr(X,Xl));

disp([ev1 ev2 ev3 ev4]);


% Plotting Histograms (optional)
figure;
subplot(2,1,1)
histogram(simBPrime(:,1),'Normalization','probability');
title('assets-to-GDP ratio');

subplot(2,1,2)
histogram(simC(:,1),'Normalization','probability');
title('Consumption');
