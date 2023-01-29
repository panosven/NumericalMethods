%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stochastic NGM + Chebyshev Approximation (collocation) of V(k,z)

% Solution method:
% - Chebyshev polynomial approximation (collocation) of consumption function [c(k,z)]
% - Gauss-Hermite quadrature approximation of integral (expectation term)

% We use GH quadrature for the expectation since the shock is drawn from 
% a Gaussian distribution (check Collard's notes on numerical integration)

% Panagiotis Veneris, U of Liverpool
% 24/11/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------- Model Primitives----------------------------------
% V(k,z) = max_{c,k'} [ c^(1-sigma)-1/(1-sigma) + E V(k',z')] 
%
% with
% c+k' = e^(z) * F(k) + (1-delta) * k  % budget constraint
% c^(-sigma) = beta*E[ (c')^(-sigma) * (a * e^(z') * (k')^(a-1) + (1-delta))] 
% log(Z') = rho_z * log(Z) + eps_z';
%-----------------------------------------------------------------


% Housekeeping
clear;
close all;
clc;
% addpath('C:\Users\pvene\Dropbox\Codes_Library\Numerical Methods\Chebyshev Collocation\Chebyshev_functions')
addpath('Chebyshev_functions')

global param 
global kpolynomial_nodesnum kpolynomial_order k_i T_x_k
global zpolynomial_nodesnum zpolynomial_order z_i T_x_z
global gh_nodesnum gh_weights gh_nodes


% Parameters
p.alpha = 0.33;  % capital share
p.beta  = 0.95;  % discount factor
p.delta = 0.1;   % capital depreciation rate 
p.sigma = 1.5;   % inverse EIS (log utility)


% Steady state (analytic)
Zss   = 1;
Kss   = (((1/p.beta)+p.delta-1)/p.alpha)^(1/(p.alpha-1)); 
Yss   = Zss*Kss^p.alpha;
Css   = Yss - p.delta*Kss;
Invss = p.delta*Kss;


% Grid for capital, and Chebyshev polynomial for capital
KGridwidth = 0.5;
p.kmin     = log((1-KGridwidth)*Kss);   % lower bound in logs
p.kmax     = log((1+KGridwidth)*Kss);   % upper bound in logs
p.Ngrid    = 30;

kpolynomial_order      = 4;                                         % order of Chebyshev polynomial for k
kpolynomial_nodesnum   = kpolynomial_order+1;                       % # of Chebyshev nodes for k
kpolynomial_chebynodes = chebyshev_nodes(kpolynomial_order);        % Chebyshev nodes (roots) for k in [-1,1]
T_x_k = compute_chebypol(kpolynomial_chebynodes,kpolynomial_order); % Chebyshev polynomials for k 
k_i = exp(transformMethod1_x_to_k(kpolynomial_chebynodes,p));    


% T_x_k in [-1,1] are the basis functions for k constructed at the Chebyshev
% nodes for k (kpolynomial_chebynodes)

% transform_x_to_k maps nodes {kpolynomial_chebynodes}_i=1...,kpolynomial_nodesnum in [-1,1] 
% to {k_i}_i=1,...kpolynomial_nodesnum in [p.kmin,p.kmax]


% Technology shock process(z):

% Under the assumption that the productivity shock is normally distributed (i.e stochastic process
% with Gaussian distribution) we can use Gauss-Hermite quadrature techniques to compute the 
% expectation term E(V(k',z')). This is called numerical integration (remember,the expectation
% is an integral).

% An easier way to compute the expectation is to discretize the process for z using Tauchen's
% method. In this way we obtain a grid for z and the associated probability
% transition matrix, and the E(V(k',z')) can be expressed as the sum of
% the probability trans. matrix * V(k',z').

sigma_e     = 0.02;   % std. deviation of stohastic process
rho_z       = 0.8;    % persistence of stohastic process
mu          = 0;      % mean of stohastic process
gh_nodesnum = 12;     % # of nodes, GH quadrature; order of GH polynomial

[gh_nodes,gh_weights] = GHernodes(gh_nodesnum); 
gh_nodes = gh_nodes * sqrt(2)* sigma_e; % GH nodes, with gh_weights as associated weights


% Productivity grid, and polynomials
p.zmin = (mu - gh_nodes(1)); % lower bound for productivity grid
p.zmax = (mu + gh_nodes(1)); % upper bound for productivity grid

zpolynomial_order      = 2;                   % order of chebyshev polynomial for productivity
zpolynomial_nodesnum   = zpolynomial_order+1; % # chebyshev nodes for productivity (technology shock)
zpolynomial_chebynodes = chebyshev_nodes(zpolynomial_order); % Chebyshev nodes (roots) for z in [-1,1]
T_x_z = compute_chebypol(zpolynomial_chebynodes,zpolynomial_order); % Cheb polynomials for z 
z_i = transformMethod1_x_to_z(zpolynomial_chebynodes,p);


% Τ_x_z in [-1,1] are the basis functions for z constructed at the
% Chebyshev nodes for z (zpolynomial_chebynodes)

% transform_x_to_z map nodes {zpolynomial_chebynodes}_i=1...,zpolynomial_nodesnum in [-1,1]
% to {z_i}_i=1,...zpolynomial_nodesnum in [p.zmin,p.zmax]

% Store parameters in matrix (for use inside functions)
param = [p.alpha p.beta p.delta p.kmin p.kmax p.zmin p.zmax rho_z sigma_e mu p.sigma]';

% Initial guess of Chebyshev weights (comes from Freund)
% theta0 = 2*ones(12,1);
theta0 = [ 
    0.7212
    0.0493
    0.0017
    0.6527
   -0.0209
   -0.0005
    0.0138
    0.0016
    0.0000
   -0.0000
   -0.0000
    0.0000
];
theta0 = theta0(:);

% Mind that the stability of the algorithm is highly dependent on the
% initial guess. For poor initial guesses, fsolve can't spot roots.

% Main loop
theta = fsolve('residuals_SNGM', theta0, [] , param);



%------ Computing Policy Functions using finer grid for accuracy---------

% Capital
kgrid_num = 100;
kgrid_fine = linspace(p.kmin,p.kmax,kgrid_num)';   % finer grid in [p.kmin,p.kmax]
vkgrid_fine = exp(kgrid_fine);
kgrid_fine_trans =  (2*(kgrid_fine-p.kmin))/(p.kmax-p.kmin) - 1; % map into [-1,1]

T_x_k = compute_chebypol(kgrid_fine_trans,kpolynomial_order); % Chebyshev polynomials for k 

% Productivity
zgrid_num = 10;
zgrid_fine = linspace(p.zmin,p.zmax,zgrid_num)';              % finer grid in [p.zmin,p.zmax]
vzgrid_fine = exp(zgrid_fine);
zgrid_fine_trans = (2*(zgrid_fine-p.zmin))/(p.zmax-p.zmin) - 1; % map into [-1,1]

T_x_z = compute_chebypol(zgrid_fine_trans,zpolynomial_order); % Cheb polynomials for z 


% Initialize
C      = [];           
Y      = [];
KPrime = [];
Inv    = [];


for iz = 1:zgrid_num
    Phi    = polynomial_comb(T_x_z(iz,:),T_x_k);              % approximated object
    C      = [C exp(Phi*theta)];                                   % consumption
    Y      = [Y exp(zgrid_fine(iz))*vkgrid_fine.^(p.alpha)];   % output
    KPrime = [KPrime Y(:,iz)+(1-p.delta)*vkgrid_fine-C(:,iz)]; % capital tomorrow
    Inv    = [Inv exp(zgrid_fine(iz))*vkgrid_fine.^(p.alpha)-C(:,iz)]; % investment
end


% Plot Capital Policy function 
figure('Name','Policy Functions: Capital')
p1 = plot(vkgrid_fine,KPrime(:,1),'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vkgrid_fine,KPrime(:,round(zpolynomial_nodesnum/2+0.1)),'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(vkgrid_fine,KPrime(:,end),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold on
p4 = plot(vkgrid_fine,vkgrid_fine,'--','color','k');
hold off
xlabel('Capital today,k','FontSize',12,'fontname','times')
ylabel('Capital tomorrow, k''','FontSize',12,'fontname','times')
grid on
legend1 = legend('Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)','$45^{\circ}$ line');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('capital policy');

% Plot Consumption Policy function 
figure('Name','Policy Functions: Consumption')
p1 = plot(vkgrid_fine,C(:,1),'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vkgrid_fine,C(:,round(zpolynomial_nodesnum/2+0.1)),'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(vkgrid_fine,C(:,end),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold off
xlabel('Capital today,K','FontSize',12,'fontname','times')
ylabel('Consumption today, C''','FontSize',12,'fontname','times')
grid on
legend1 = legend('Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('consumption policy');

% Plot investment policy function 
figure('Name','Policy Function: Investment')
p1 = plot(vkgrid_fine,Inv(:,1),'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vkgrid_fine,Inv(:,round(zpolynomial_nodesnum/2+0.1)),'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(vkgrid_fine,Inv(:,end),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold off
xlabel('Capital today, K','FontSize',12,'fontname','times')
ylabel('Investment, I''','FontSize',12,'fontname','times')
grid on
legend1 = legend('Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('investment policy');


     
     
     
