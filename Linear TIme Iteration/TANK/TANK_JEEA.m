% Basic TANK (Gali et al, 2007, JEEA)
% Main differences: 
%    1. investmet adjustment costs 
%    2. Rotemberg pricing
%    3. Taylor rule 
% Solved with LTI (Rendahl, 2017)
% Panagiotis Veneris, U of Liverpool
% 11/04/2023

clear variables
% close all
clc

% Specify shock
% shock = 'MP';
shock = 'Government Spending';


% Parameters (#19,quarterly calibration)
p.beta    = 0.99;     % hh discount factor
p.alpha   = 1/3;      % capital share
p.epsilon = 6;        % elasticity of subs. among differentiated goods
p.delta   = 0.025;    % depreciation rate
p.sigma   = 1;        % risk aversion
p.phi     = 0.2;      % inverse Frisch
p.eta     = 0.5;      % fraction of rule-of-thumb consumers
p.phipi   = 1.5;      % mp response to devs of inflation from target
p.phiy    = 0.0;      % mp response to devs of output from target
p.kappaI  = 2.48;     % investment adjustment cost parameter (as in CEE)
p.rhoa    = 0.9;      % persistence of productivity shock
p.rhog    = 0.9;      % persistence of gov. spending shock
p.rhonu   = 0.0;      % persistence of mp shock (no persistence = inertia)
p.phig    = 0.1;      % tax reaction to devs of gov.spending from target
p.phid    = 0.33;     % tax reaction to devs of debt from target
calvo     = 0.75;     % calvo parameter
p.sig_nu  = 0.0025;   % size of monetary policy shock (% deviation)
p.sig_a   = 0.01;     % size of technology shock (% deviation)
p.sig_g   = 0.01;     % size of government spending  

% Calibration targets
p.g       = 0.2;      % public spending over gdp ratio:(G_t-G_ss)/Y=0.2, similar to Gali 2007
p.D       = 1;        % annual public debt to gdp ratio (in general, d is equivalent to b in Gali's paper)

% Deterministic S.S values of model variables(#21 variables)
s.g_ss       = p.g;      % s.s gov. spending
s.Pi_ss      = 1;        % s.s gross inflation rate
s.y_ss       = 1;        % s.s aggregate output (gdp); productivity (a) calibrated to yield y_ss = 1
s.h_ss       = 1/3;      % s.s aggregate hours worked
s.Rreal_ss   = 1/p.beta; % s.s real interest rate
s.Rnom_ss    = s.Pi_ss/p.beta; % s.s real interest rate (as net inflation is zero in s.s)
s.q_ss       = 1;              % marginal value of investment
s.rk_ss      = 1/p.beta - (1-p.delta);  % rental rate of capital
s.MC_real_ss = (p.epsilon-1)/p.epsilon; % s.s real marginal cost
s.k_ss       = p.alpha*s.MC_real_ss*s.y_ss/s.rk_ss; % s.s capital stock
s.w_real_ss  = (1-p.alpha)*s.MC_real_ss*s.y_ss / s.h_ss; % s.s real wage
s.invest_ss  = p.delta*s.k_ss;             % s.s aggregate investment
s.c_ss       = s.y_ss - s.invest_ss - s.g_ss; % s.s aggregate consumption
s.t_r_ss     = s.w_real_ss * s.h_ss - s.c_ss; % s.s taxes on rule-of-thumb consumers
s.d_ss       = 4*p.D;                         % s.s public debt to gdp ratio
s.t_o_ss     = (s.g_ss + (1/p.beta - 1) * s.d_ss - p.eta*s.t_r_ss)/(1-p.eta); % s.s taxes on optimizing consumers
s.a_ss       = s.y_ss / (s.k_ss^(p.alpha)*s.h_ss^(1-p.alpha)); % s.s productivity (calibrated to yield y_ss=1)
s.c_o_ss     = s.c_ss; % s.s consumption of optimizing hhs (see Gali et al,2007,JEEA, pg243)
s.c_r_ss     = s.c_ss; % s.s consumption of rule-of-thumb hhs (see Gali et al,2007,JEEA, pg243)
s.h_o_ss     = s.h_ss; % s.s hours worked of optimizing hhs (see Gali et al,2007,JEEA, pg243)
s.h_r_ss     = s.h_ss; % s.s hours worked of rule-of-thumb hhs (see Gali et al,2007,JEEA, pg243)

% Parameters(continued)
p.kappaL     = s.w_real_ss / (s.h_ss^(p.phi)*s.c_ss^(p.sigma)); % disutility from labor parameter
p.kappaP     = (p.epsilon-1)*calvo/(s.Pi_ss^2*(1-calvo)*(1-p.beta*calvo)); % adjustment cost parameter

% Set up non-linear stochastic system symbolically and separate forward,lagged and contemporaneous variables (#21 variables)
syms c Rreal rk w_real h y k q invest Rnom MC_real Pi g a c_o c_r h_o h_r t_o t_r d % contemporaneous                  
syms lc lRreal lrk lw_real lh ly lk lq linvest lRnom lMC_real lPi lg la lc_o lc_r lh_o lh_r lt_o lt_r ld  % lagged    
syms fc fRreal frk fw_real fh fy fk fq finvest fRnom fMC_real fPi fg fa fc_o fc_r fh_o fh_r ft_o ft_r fd  % forward   

% Write non-linear system equations in symbolic format (#21 equations in #21 unknown variables)
%---households---
e.eq1  = c_o.^(-p.sigma) - p.beta.* ((fc_o).^(-p.sigma) .* Rnom./fPi);       % Euler_b for optimizing hhs
e.eq2  = 1 - ((p.beta./q).* (fc_o./c_o).^(-p.sigma).*(frk+(1-p.delta).*fq)); % Euler_k for optimizing hhs
e.eq3  = p.kappaL.*(h_o).^(p.phi) - w_real.*(c_o).^(-p.sigma);               % Euler_h for optimizing hhs
e.eq4  = 1 - q.*(1-p.kappaI/2.*((invest./linvest)-1).^2 - p.kappaI.*(invest./linvest).*((invest./linvest)-1))- (p.kappaI*p.beta).*fq.*(fc_o./c_o).^(-p.sigma).*((finvest./invest).^2).*((finvest./invest)-1); % FOCi for optimizing hhs since those HHs are the owners of firms  
e.eq5  = p.kappaL.*h_r.^(p.phi) - w_real.*c_r.^(-p.sigma);                   % Euler_h for rule of thumb hhs
e.eq6  = k - (1-p.delta).*lk - (1-(p.kappaI/2).*((invest./linvest)-1).^2).*invest;    % LOM of capital
%---firms---
e.eq7  = y - a.*lk.^(p.alpha).*h.^(1-p.alpha);    % aggregate output
e.eq8  = w_real.*h - MC_real.*(1-p.alpha).*y;     % MPL
e.eq9  = rk.*lk - MC_real.*p.alpha.*y;            % MPK
e.eq10 = (Pi-s.Pi_ss).*Pi - p.beta.* (((fc_o./c_o).^(-p.sigma)).*fPi.*(fPi-s.Pi_ss).*(fy./y)) - (p.epsilon/p.kappaP).*(MC_real - (p.epsilon-1)/p.epsilon); % NKPC
%---shocks---
e.eq11 = -log(g) + (p.rhog.*log(lg)) + (1-p.rhog).*log(s.g_ss); % AR(1) in logs gov.spending shock process 
e.eq12 = -log(a) + (p.rhoa.*log(la)) + (1-p.rhoa).*log(s.a_ss); % AR(1) in logs technology shock process
%---market clearing (market for capital clears residually from Walras law)---
e.eq13 = y-c-invest-g-(p.kappaP/2).*((Pi-s.Pi_ss).^2).*y; % final good's market clearing condition (aggregate resource constraint)
e.eq14 = c - (1-p.eta).*c_o - p.eta.*c_r;  % aggregate consumption
e.eq15 = h - (1-p.eta).*h_o - p.eta.*h_r;  % labor market clearing condition (aggregate labor supply or hours worked)
%---auxiliary---
e.eq16 = c_r - w_real.*h_r + t_r;  % budget constraint of rule of thumb hhs
e.eq17 = Rnom - (Rreal.*fPi);      % fisher equation (nominal interest rate)
%---policy (tax rules in linear form)---
e.eq18 = -Rnom./s.Rnom_ss + ((Pi./s.Pi_ss).^(p.phipi).*(y./s.y_ss).^(p.phiy)).^(1-p.rhonu).*(lRnom./s.Rnom_ss).^(p.rhonu); % Taylor rule (# Gali 2007,eq18)   
e.eq19 = t_r - s.t_r_ss - p.phig.*(g-s.g_ss) - p.phid.*(ld-s.d_ss);  % taxation of rule of thumb hhs (Gali 2007, eq20; calibrated to yield c_o=c_r=c)
e.eq20 = t_o - s.t_o_ss - p.phig.*(g-s.g_ss) - p.phid.*(ld-s.d_ss);  % taxation of optimizing hhs (Gali 2007, eq20;calibrated to yield c_o=c_r=c); 
e.eq21 = g + (lRnom./Pi).*ld - (1-p.eta).*t_o - p.eta.*t_r - d;      % gov. BC    

% Initialize vectors and collect equations and variables (all variables should appear in the same order across matrices)
m.Vxss = [s.c_ss,s.Rreal_ss,s.rk_ss,s.w_real_ss,s.h_ss,s.y_ss,s.k_ss,s.q_ss,s.invest_ss,s.Rnom_ss,s.MC_real_ss,s.Pi_ss,s.g_ss...     % 1x21 column vector of S.S values of variables                
          s.a_ss,s.c_o_ss,s.c_r_ss,s.h_o_ss,s.h_r_ss,s.t_o_ss,s.t_r_ss,s.d_ss]; 
m.Vxf = [fc fRreal frk fw_real fh fy fk fq finvest fRnom fMC_real fPi fg fa fc_o fc_r fh_o fh_r ft_o ft_r fd]; % 1x21 vector of forward symbolic variables   
m.Vxl = [lc lRreal lrk lw_real lh ly lk lq linvest lRnom lMC_real lPi lg la lc_o lc_r lh_o lh_r lt_o lt_r ld ]; % 1x21 vector of lagged symbolic variables    
m.Vxc = [c Rreal rk w_real h y k q invest Rnom MC_real Pi g a c_o c_r h_o h_r t_o t_r d ]; % 1x21 vector of contemporaneous symbolic variables                 

m.system     = [e.eq1;e.eq2;e.eq3;e.eq4;e.eq5;e.eq6;e.eq7;e.eq8;e.eq9;e.eq10;e.eq11;e.eq12;e.eq13;e.eq14;... % system of equations in row vector form                
                e.eq15;e.eq16;e.eq17;e.eq18;e.eq19;e.eq20;e.eq21];
m.Vxss_total = repmat(m.Vxss,3,1); % 3x21 matrix
m.Vxss_total = m.Vxss_total(:)';   % 1x63 column vector S.S values for all variables (22 contemp.,22 lagged,22 forward), 
varnum       = size(m.Vxss,2);     % number of variables

% Auxiliary matrix
m.V_auxil = reshape([m.Vxl;m.Vxc;m.Vxf],size(m.Vxc,1),[]);

% Linearize system using Jacobian linearization 
m.C = jacobian(m.system,m.Vxl);
m.C = double(subs(m.C,m.V_auxil,m.Vxss_total));
m.B = jacobian(m.system,m.Vxc);
m.B = double(subs(m.B,m.V_auxil,m.Vxss_total));
m.A = jacobian(m.system,m.Vxf);
m.A = double(subs(m.A,m.V_auxil,m.Vxss_total));
% the state-space system is expressed as A*X(t+1)+B*X(t)+C*X(t-1)

% Convert the linearized system into log-linear (variables in p.p deviations from their S.S values)
m.L = ones(21,1)*m.Vxss;                                                                            
m.A = m.A.*m.L;
m.B = m.B.*m.L;
m.C = m.C.*m.L;

% LTI solver
p.Tol   = 1e-13;    
metric  = 0.1;         
it      = 0;       

P = 0;  % initial guess for P
tic
while metric>p.Tol
    it     = it+1;   
    P      = -(m.A*P+m.B)\m.C;                              
    metric = max(max(abs(m.A*P*P+m.B*P+m.C)));                                                                                                
    fprintf('iteration = %d, metric = %e\n', it, metric);                
end
toc
% Check for stable equilibrium (all eigenvalues of solvent P must be < 1)
eig_P = max(abs(eig(P)));
if eig_P < 1
    disp('Stable solution exists')  % all eigenvalues must be < 1 in absolute value
end
Q  = -inv(m.A*P+m.B);

% Strings for shock specification
if strcmp(shock,'MP')
    size   = p.sig_nu; % size of monetary policy shock 
    eqnumb = 18;     % equation that we shock            

else strcmp(shock,'Government Spending')
    size   = p.sig_g;  % size of government spending shock
    eqnumb = 11;     % equation that we shock          
end
    
% IRFs
p.T = 15;             % horizon of responses
u = zeros(varnum,1);  % initialize vector of disturbances
u(eqnumb) = size; % size of shock


x(:,1) = Q*u;

for t = 1:p.T-1
    x(:,t+1) = P*x(:,t);
end

% Plots
figure('name',shock)
subplot(3,2,1)
plot(x(10,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % nominal interest rate   
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Nominal Interest Rate')
box on

subplot(3,2,2)
plot(x(12,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % inflation               
xlabel('Period','FontSize',12,'fontname','times'); 
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Inflation')
box on

subplot(3,2,3)
hold on
plot(x(2,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % real interest rate     
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Real Interest Rate')
box on

subplot(3,2,4)
plot(x(1,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % consumption            
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Consumption')
box on

subplot(3,2,5)
hold on
plot(x(5,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % labor supply           
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Labor Supply')
box on

subplot(3,2,6)
hold on
plot(x(6,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % output                
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Output')
box on

figure('name',shock)
subplot(3,2,1)
plot(x(7,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % aggregate capital     
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Capital')
box on

% subplot(3,2,2)
% plot(x(9,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % tobin's q
% xlabel('Period','FontSize',12,'fontname','times');
% ylabel('% dev. from SS','FontSize',12,'fontname','times')
% title('Tobin Q')
% box on

subplot(3,2,2)
hold on
plot(x(9,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % investment          
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Investment')
box on

subplot(3,2,3)
plot(x(15,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % consumption optimizers  
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Consumption Optimizers')
box on

subplot(3,2,4)
hold on
plot(x(16,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % consumption HtM (rule of thumb)         
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Consumption HtM')
box on

subplot(3,2,5)
hold on
plot(x(17,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % labor supply optimizers 
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Labor Supply Optimizers')
box on

subplot(3,2,6)
hold on
plot(x(18,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % labor supply HtM      
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Labor Supply HtM')
box on
