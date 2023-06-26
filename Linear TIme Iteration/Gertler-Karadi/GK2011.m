% Replication of GK(2011,JME)
% Main differences:
%    1. No retail firms
%    2. No variable capital utilization
%    3. No cost of government intermediation
%    4. No price indexation
%    5. No habit persistence
% Panagiotis Veneris, U of Liverpool
% 15/05/2023

clear variables
% close all
clc

% Specify shock
% shock = 'MP';
% shock = 'Government Spending';
% shock = 'Technology';
% shock = 'capital quality';
shock = 'QE';

% Structural Parameters (#21)
p.alpha   = 0.33;           % capital share
p.beta    = 0.99;           % discount factor
p.delta   = 0.025;          % depreciation rate
p.epsilon = 6;              % elasticity of substitution btw differentiated goods
p.eta     = 2.48;           % investment adjustment cost 
p.g       = 0.2;            % share of public spending in ss (as in GK 2011, g/y =g as y=1 in SS)
p.phi     = 1;              % inverse of Frisch elasticity
p.phipi   = 1.5;            % mp response to inflation
p.phiy    = 0.125;          % mp response to output
p.phiqe   = 0.0;            % qe response to spread                            
p.QE      = 0.1;            % share of government assets over total bonds
p.rhoa    = 0.9;            % persistence of productivity shock
p.rhog    = 0.9;            % persistence of gov.spending shock
p.rhonu   = 0.8;            % peristence of mp shock
p.rhoqe   = 0.9;            % qe persistence
p.rhoqe   = 0.9;            % qe persistence
p.rhoz    = 0.66;           % capital shock persistence
p.sig     = 2;              % inverse EIS
p.theta   = 0.972;          % survivial rate of bankers (as in GK 2011
calvo     = 0.66;           % price rigidity in calvo framework
p.sig_nu  = 0.0025;         % size of monetary policy shock (% deviation)
p.sig_a   = 0.01;           % size of technology shock (% deviation)
p.sig_g   = 0.01;           % size of government spending  
p.sig_qe  = 0.1;            % size of qe shock
p.sig_z   = 0.05;           % size of capital quality shock


% Deterministic S.S values of model variables(#23 variables: discard (varrho) relative to dynare due to Eulers instead of FOCs)
s.g_ss       = p.g;                              % s.s gov. spending
s.Pi_ss      = 1;                                % s.s inflation targeting (quarterly calibration)
s.y_ss       = 1;                                % s.s output
s.l_ss       = 1/3;                              % s.s labor supply (hours of work)
s.sp_ss      = 0.0025;                           % s.s credit spread (spread between loan and deposit rates)
s.lev_ss     = 4;                                % s.s leverage multiple
s.r_ss       = 1/p.beta;                         % s.s real interest rate
s.ir_ss      = s.Pi_ss/p.beta;                   % s.s nominal interest rate
s.q_ss       = 1;                                % s.s price of capital  
s.rk_ss      = s.r_ss+s.sp_ss;                   % s.s loan rate (gross real rate of return on capital) or s.s bank int.rate
s.rR_ss      = s.rk_ss-(1-p.delta);              % s. rental rate of capital
s.MC_real_ss = (p.epsilon-1)/p.epsilon;          % s.s real marginal costs
s.k_ss       = p.alpha*s.MC_real_ss*s.y_ss/s.rR_ss;          % s.s capital stock
s.w_real_ss  = (1-p.alpha)*s.MC_real_ss*s.y_ss/s.l_ss;       % s.s real wage
s.invest_ss  = p.delta*s.k_ss;                               % s.s investment
s.c_ss       = s.y_ss-s.invest_ss-s.g_ss;                    % s.s consumption
s.a_ss       = s.y_ss/(s.k_ss^(p.alpha)*s.l_ss^(1-p.alpha)); % s.s productivity (TFP)
s.sg_ss      = p.QE*s.k_ss;                                  % s.s government's assets
s.sf_ss      = s.k_ss-s.sg_ss;                               % s.s bank assets
s.n_ss       = s.sf_ss/s.lev_ss;                             % s.s net worth
s.bF_ss      = s.sf_ss-s.n_ss;                               % s.s HH deposits to financial institutions(banks)
s.ksi_ss     = (1-p.theta)/(1-p.theta*p.beta*(s.sp_ss*s.lev_ss+s.r_ss)); % s.s bank discount factor
s.zeta_ss    = 1;                                            % s.s capital quality

% Parameters (continued)
p.chi     = s.w_real_ss/(s.l_ss^(p.phi)*s.c_ss^(p.sig));  % labor disutility parameter without habit persistence
p.omega   = s.n_ss*(1-p.theta*(s.sp_ss*s.lev_ss+s.r_ss))/s.sf_ss;       % s.s tranfer to new bankers (remember q=1 in s.s)
p.lambda  = s.ksi_ss/s.lev_ss+p.beta*s.ksi_ss*s.sp_ss;                  % s.s fraction of divertable assets (in SS r = 1/beta)
p.gammaP  = (p.epsilon-1)*calvo/(s.Pi_ss^2*(1-calvo)*(1-p.beta*calvo)); % adjustment cost parameter

% Set up non-linear stochastic system symbolically and separate forward,lagged and contemporaneous variables (#23 variables)
syms c  rR  w_real  l  y  k  q  invest  r  ir  Pi  MC_real  n  bF  ksi  lev  sp  rk  sf  sg  g  a  zeta  % #23 contemporaneous variables
syms lc lrR lw_real ll ly lk lq linvest lr lir lPi lMC_real ln lbF lksi llev lsp lrk lsf lsg lg la lzeta % #23 lagged variables
syms fc frR fw_real fl fy fk fq finvest fr fir fPi fMC_real fn fbF fksi flev fsp frk fsf fsg fg fa fzeta % #23 forward variables

% Write non-linear system equations in symbolic format (#23 equations in #23 unknown variables)
%---households---
e.eq1 = c.^(-p.sig) - p.beta.*((fc).^(-p.sig).*ir./fPi); % Euler_bonds
e.eq2 = w_real.* c.^(-p.sig) - p.chi.*(l.^(p.phi));      % Euler_labor
%---capital firms---
e.eq3 = 1 - q.*(1-(p.eta/2).*(((invest./linvest)-1).^2) - p.eta.*(invest./linvest).*((invest./linvest)-1))- (p.eta*p.beta).*fq.*((fc./c).^(-p.sig)).*((finvest./invest).^2).*(finvest./invest - 1); % FOC for investment
e.eq4 = k - (1-p.delta).*zeta.*lk - invest.*(1- (p.eta/2).*((invest./linvest - 1).^2)); % LOM for capital (one lag for k since K_t+1 is due to a decision about K made at period t)
%---intermediate firms---
e.eq5 = y - a.*(zeta.*lk).^(p.alpha).*l.^(1-p.alpha);  % Production function (eq.20)
e.eq6 = w_real.*l - (1-p.alpha).*MC_real.*y;           % FOC for labor demand (eq.26)
e.eq7 = rR.*lk - p.alpha.*MC_real.*y;                  % FOC for capital demand (eq.25)
e.eq8 = rR - rk.*lq + (1-p.delta).*zeta.*q;            % Rental rate of capital (eq.24)
e.eq9 = (Pi-s.Pi_ss).*Pi - p.beta.* (((fc./c).^(-p.sig)).*fPi.*(fPi-s.Pi_ss).*(fy./y)) - (p.epsilon/p.gammaP).*(MC_real - (p.epsilon-1)/p.epsilon); % NKPC
%---banking sector---
e.eq10 = lev - ((p.beta.*((fc./c).^(-p.sig)).*fksi.*(ir./fPi))./ (p.lambda - (p.beta.* ((fc./c).^(-p.sig)).*fksi.*(frk-ir./fPi)))); % Leverage multiple
e.eq11 = lev - (q.*sf)./n;   % Bank's demand for assets (demand for capital claims)
e.eq12 = q.*sf - bF - n;     % Bank's balance sheet (bF enters with one period lag since the deposits in t+1 have been decided at period t)    
e.eq13 = n - p.theta.*((rk-lir./Pi).*llev+(lir./Pi)).*ln - p.omega.*q.*lsf; % LOM of bank's net worth                                        
e.eq14 = ksi - (1-p.theta) - (p.theta*p.beta).* (((fc./c).^(-p.sig)).*fksi.*((frk - ir./fPi).*lev + ir./fPi)); % banks' discount factor
%--shock processes---
e.eq15 = -log(a) + (p.rhoa.*log(la)) + (1-p.rhoa).*log(s.a_ss); % AR(1) in logs technology shock process
e.eq16 = -log(g) + (p.rhog.*log(lg)) + (1-p.rhog).*log(s.g_ss); % AR(1) in logs gov.spending shock process 
e.eq17 = -log(zeta) + (p.rhoz.*log(lzeta)); % AR(1) in logs capital quality shock process 
%---market clearing conditions---
e.eq18 = y-c-invest-g-(p.gammaP/2).* (((Pi-s.Pi_ss).^2).*y);  % Goods-market clearing (aggregate resource constraint)
e.eq19 = q.*k - q.*sf - q.*sg;                             % Bonds(deposits)-market clearing; k is with one period lag                    
%---auxiliary---
e.eq20 = sp - (frk - (ir./fPi));  % Spread method 1
e.eq21 = r - (ir./fPi);           % Fisher equation (nominal interest rate) method 1
%---policy rules---
e.eq22 =  -(ir./s.ir_ss) + (lir./s.ir_ss).^(p.rhonu).*((Pi./s.Pi_ss).^(p.phipi).*(y./s.y_ss).^(p.phiy)).^(1-p.rhonu); % Taylor rule 
e.eq23 =  -(sg./s.sg_ss) + (lsg./s.sg_ss).^(p.rhoqe).*((sp./s.sp_ss).^(p.phiqe)).^(1-p.rhoqe);  % QE rule

% Initialize vectors and collect equations and variables (all variables should appear in the same order)
m.Vxss = [s.c_ss,s.rR_ss,s.w_real_ss,s.l_ss,s.y_ss,s.k_ss,s.q_ss,s.invest_ss,s.r_ss,s.ir_ss,s.Pi_ss,s.MC_real_ss,s.n_ss,... % 1x23 column vector of S.S values of symbolic variables 
          s.bF_ss,s.ksi_ss,s.lev_ss,s.sp_ss,s.rk_ss,s.sf_ss,s.sg_ss,s.g_ss,s.a_ss,s.zeta_ss];      
m.Vxf  = [fc frR fw_real fl fy fk fq finvest fr fir fPi fMC_real fn fbF fksi flev fsp frk fsf fsg fg fa fzeta]; % 1x23 column vector of forward symbolic variables 
m.Vxl  = [lc lrR lw_real ll ly lk lq linvest lr lir lPi lMC_real ln lbF lksi llev lsp lrk lsf lsg lg la lzeta]; % 1x23 column vector of lagged symbolic variables 
m.Vxc  = [c  rR  w_real  l  y  k  q  invest  r  ir  Pi  MC_real  n  bF  ksi  lev  sp  rk  sf  sg  g  a  zeta ]; % 1x23 column vector of contemporaneous symbolic variables                 

m.system     = [e.eq1;e.eq2;e.eq3;e.eq4;e.eq5;e.eq6;e.eq7;e.eq8;e.eq9;e.eq10;e.eq11;e.eq12;e.eq13;e.eq14;e.eq15;...
                e.eq16;e.eq17;e.eq18;e.eq19;e.eq20;e.eq21;e.eq22;e.eq23];
m.Vxss_total = repmat(m.Vxss,3,1); % 3x23 matrix
m.Vxss_total = repmat(m.Vxss,3,1); % 3x23 matrix
m.Vxss_total = m.Vxss_total(:)';   % 1x69 column vector S.S values for all variables (23 contemp.,23 lagged,23 forward), 
varnum       = size(m.Vxss,2);     % number of variables

% Auxiliary matrix
m.V_auxil = reshape([m.Vxl;m.Vxc;m.Vxf],size(m.Vxc,1),[]);

% Linearize system using Jacobian linearization 
m.A = jacobian(m.system,m.Vxl);                 
m.A = double(subs(m.A,m.V_auxil,m.Vxss_total));
m.B = jacobian(m.system,m.Vxc);
m.B = double(subs(m.B,m.V_auxil,m.Vxss_total));
m.C = jacobian(m.system,m.Vxf);
m.C = double(subs(m.C,m.V_auxil,m.Vxss_total));
% the state-space system is expressed as A*X(t-1)+B*X(t)+C*X(t+1)

% Convert the linearized system into log-linear (variables in p.p deviations from their S.S values)
m.L = ones(23,1)*m.Vxss;                                                                            
m.A = m.A.*m.L;
m.B = m.B.*m.L;
m.C = m.C.*m.L;

% LTI solver
metric  = 0.1;    % initial value      
it      = 0; 
maxit   = 500;
p.Tol   = 1e-11; 

P       = 0;  % initial guess for dominant solvent
S       = 0;  % initial guess for minimal solvent

eye = zeros(varnum,varnum);
for i = 1:varnum
    eye(i,i) = 1; % identity matrix
end
I = eye*p.Tol;
m.Ch = m.C;                 % forward
m.Bh = (m.B+m.C*2*I);       % contemp
m.Ah = (m.C*I^2+m.B*I+m.A); % lagged 

if rcond(m.A)<1e-16
    disp('Matrix m.A is singular')
end

tic
while metric>p.Tol 
    it     = it+1;
    P = -(m.Bh+m.Ch*P)\m.Ah;
    S = -(m.Bh+m.Ah*S)\m.Ch;
    metric1 = max(max(abs(m.Ah+m.Bh*P+m.Ch*P*P)));
    metric2 = max(max(abs(m.Ah*S*S+m.Bh*S+m.Ch)));
    metric = max(metric1, metric2);  
    fprintf('iteration = %d, metric = %e\n', it, metric);
end
% Check for stable equilibrium (all eigenvalues of solvent P must be < 1)
eig_F = max(abs(eig(P))); 
eig_S = max(abs(eig(S)));

if eig_F>1 || eig_S>1 || p.Tol>1-eig_S
    disp('Conditions of Proposition 3 violated')
end
P = P+I;
Q = -inv(m.B+m.C*P);

% Strings for shock specification
if strcmp(shock,'MP')
    size = p.sig_nu; % size of monetary policy shock 
    eqnumb = 22;    % equation that we shock            
elseif strcmp(shock,'Technology')
    size = p.sig_a;  % size of technology shock
    eqnumb = 15;     % equation that we shock
elseif strcmp(shock,'QE')
    size = p.sig_qe;
    eqnumb = 23;
elseif strcmp(shock,'capital quality')
    size = p.sig_z;
    eqnumb = 17;
else strcmp(shock,'Government Spending')
    size = p.sig_g;  % size of government spending shock
    eqnumb = 16;     % equation that we shock          
end
    
% IRFs
p.T = 20;             % horizon of responses
u = zeros(varnum,1);  % initialize vector of disturbances
u(eqnumb) = size;     % size of shock

x(:,1) = Q*u;

for t = 1:p.T-1
    x(:,t+1) = P*x(:,t);
end

% Plots
figure('name',shock)
subplot(3,3,1)
hold on
plot(x(5,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % output                
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Output')
box on; grid on;

subplot(3,3,2)
plot(x(1,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % consumption            
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Consumption')
box on; grid on;

subplot(3,3,3)
plot(x(8,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % investment 
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Net Investment')
box on; grid on;

subplot(3,3,4)
plot(x(6,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % capital               
xlabel('Period','FontSize',12,'fontname','times'); 
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Capital')
box on; grid on;

subplot(3,3,5)
hold on
plot(x(7,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % Tobin's q    
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Tobin Q')
box on; grid on;

subplot(3,3,6)
hold on
plot(x(11,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % inflation           
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Inflation')
box on; grid on;

subplot(3,3,7)
hold on
plot(x(17,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % spread          
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Spread')
box on; grid on;

subplot(3,3,8)
hold on
plot(x(16,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % leverage          
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Leverage')
box on; grid on;

subplot(3,3,9)
hold on
plot(x(13,:),'-o','linewidth',1.5,'color',[0 0.4470 0.7410]); % net worth         
xlabel('Period','FontSize',12,'fontname','times');
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title('Net Worth')
box on; grid on;
