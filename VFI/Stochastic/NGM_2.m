% Problem Set 6: Advanced Macroeconomics
% Stohastic NGM: VFI + Schumaker spline
% Panagiotis Veneris, U of Liverpool
% 11/1/2021 

% NOTE 1: This is the stohastic NGM, that is, there is aggregate
% uncertainty through (aggregate) shocks in the productivity of capital (z)

clear;
close all; 
clc;

% Denote global variables
global beta alpha kGrid k_init V n delta z_init Pi pp z 

alpha   = 0.3;
beta    = 0.96;     % discount factor
tol     = 1e-6;     % tolerance level
maxit   = 500;      % max # of iterations
dif     = 0.1;
it      = 0;        % iteration counter
n       = 9;        % # gridpoints for productivity(exog.state variable)
sigma_e = 0.007;    % std. deviation of stohastic process
rho     = 0.95;     % persistence of stohastic process
mu      = 0;        % mean of stohastic process
w       = 3 ;       % how many std. dev. above/below mean
delta   = 0.1;      % depreciation rate


% Grid for state variable (k)
kGridnum = 30;                            % # of gridpoints for capital k
kmin     = 1;
kmax     = 5;
kGrid    = linspace(kmin,kmax,kGridnum)'; % equally-spaced grid for k (30x1 matrix)


% Grid for state variable (z)
[Pi,z] = tauchen(n,mu,rho,sigma_e,w);   % z : grid for productivity
                                        % Pi: prob. trans. matrix

% Initial guess for value function
V = zeros(kGridnum,n);             % V_0(k,z) (30x1 matrix)


% Initialize auxiliary matrices
V1     = zeros(kGridnum,n);
K11    = zeros(kGridnum,n);
k_init = zeros(kGridnum,n);
z_init = zeros(kGridnum,n);
pp     = cell(n,1);

tic
while dif>tol && it < maxit 
   it = it+1;
   
    for i_p = 1:n
        pp{i_p} = schumaker(V(:,i_p), [] ,kGrid);     % compute schumaker interpolation coefficients (coeff of polynomial) for each z
    end
   
    for i_p=1:n                                        % loop over grid for z (productivity)
       for i=1:kGridnum                                % loop over grid for assets           
           k_init = kGrid(i,1);                        % for each i in the state space get the initial value of k_t (i=1 k1, i=2 k2 etc...)
           z_init = z(i_p,1);                          % for each j in the state space get the initial value of z_t (j=1 z1, j=2 z2 etc...)
           K1 = fminbnd(@val_stoch,kmin,kmax);         % for each value of ((k_t,z_t) (k1,z1),(k2,z2) etc...) find k' that max the Bellman equation
           V1(i,i_p) = -val_stoch(K1);                 % collect the optimized values into a new value function called V1
           K11(i,i_p) = K1;                            % find policy function for capital (optimal k')-->function of 2 state variables (k,z)       
       end
    end
     
   dif = norm(V1-V);
   V=V1;

   %fprintf('Iteration: %d\tConv. tol_pol.: %g\n',[it,dif]);
    fprintf('iteration = %d, dif = %e\n', it, dif);  
end
toc
% with two state variables fminbnd computes the optimal k' for n_z X n_k different combinations
% which makes the code much slower --> Curse of Dimensionality

for i_p = 1:n
    kapital(:,i_p) = K11(:,i_p);
    c(:,i_p) = exp(z(i_p,1)).*(kGrid.^alpha)+((1-delta).*kGrid)-kapital(:,i_p); %need to check policy function
    y_output(:,i_p) = exp(z(i_p,1)).*(kGrid.^alpha);
end

%********** Alternative way to get output ********************************
%y_output1 = exp(z)*(kGrid.^alpha)';
%*************************************************************************

figure;
 subplot(2,2,1);
 hold on;
 plot(kGrid,V1(:,1),'c-*','Linewidth',2);
 plot(kGrid,V1(:,n),'b*-','Linewidth',2);
 hold off;
 xlabel('k');
 ylabel('V(k)');
 legend('n=1','n=9','Location','best');
 title('Value Funtion');
 grid on;
 
 %figure
 subplot(2,2,2);
 hold on;
 plot(kGrid,y_output(:,1),'c-*','Linewidth',2);
 plot(kGrid,y_output(:,n),'b-*','Linewidth',2);
 hold off;
 xlabel('k');
 ylabel('y');
 legend('n=1','n=9','Location','best');
 title('Output');
 grid on;
 
 %*************** Plot of output for alternative way **********************
%  figure
%  hold on;
%  plot(kGrid,y_output1(1,:),'c-*','Linewidth',2);
%  plot(kGrid,y_output1(n,:),'b-*','Linewidth',2);
%  hold off;
%  xlabel('k');
%  ylabel('y');
%  legend('n=1','n=9','Location','best');
%  title('Output');
%  grid on;
 %*************************************************************************

%figure;
 subplot(2,2,3);
 hold on;
 plot(kGrid,K11(:,1),'c-*','Linewidth',2);
 plot(kGrid,K11(:,n),'b*-','Linewidth',2)
 hold off;
 xlabel('k');
 ylabel('g(k)');
 legend('n=1','n=9','Location','best');
 title('Capital Policy Funtion');
 grid on;
 
 %figure
 subplot(2,2,4)
 hold on;
 plot(kGrid,c(:,1),'c-*','Linewidth',2);
 plot(kGrid,c(:,n),'b*-','Linewidth',2);
 xlabel('k');
 ylabel('c');
 legend('n=1','n=9','Location','best');
 title('Consumption Policy Function');
 grid on;


function val=val_stoch(kprime)   % kprime is k'
 
% pp = schumaker(y,[],x) : calculates Schumaker's interpolation coefficients (pp)
% x :  is the vector of gridpoints 
% y :  is the vector of values of the value-function

% yq = ppval(pp,xq) 
% xq : vector of gridpoints at which i want to evaluate the values of value function
% yq : vector of interpolated values
% ppval : returns the values at the points xq of the piecewise polynomial contaiined in pp
 
 global beta alpha k_init n z_init delta Pi pp  
          
  for i_p = 1:n                                                 % loop over productivity grid
    g(:,i_p) = ppval(pp{i_p},kprime);                           % mind that inside the parentheses we need pp WITH curly brackets                                                               % returns the values of the polynomial pp evaluated at k=k'
  end
    c = exp(z_init)*(k_init^alpha)+(1-delta)*k_init-kprime;     % budget constraint(solved wrt c)
  
   if c<0
       val = -888888888888888888-800*abs(c);                    % keeps it from going negative
   else         
       for i_p=1:n
          val= log(c) + beta*(Pi(i_p,:)*g');
       end  
  end
         val = -val;                                            % make it negative since we're maximizing and code is to minimize.  
end


function [Pi,z] = tauchen(n,mu,rho,sigma_e,w)

   z    = zeros(n,1);      % nodes for stohastic process (possible different states of productivity)
   Pi   = zeros(n,n);      % prob. transition matrix
   a    = (1-rho)*mu;
   z(n) = w*(sigma_e/sqrt(1-rho^2));
   z(1) = -z(n);
   step = (z(n)-z(1))/(n-1); 

   for i=2:(n-1)
       z(i) = z(1) + step * (i - 1);
   end 
   z = z + a / (1-rho);

   for j = 1:n
       for k = 1:n
           if k == 1
               Pi(j,k) = normcdf((z(1) - a - rho * z(j) + step / 2) / sigma_e);
           elseif k == n
               Pi(j,k) = 1 - normcdf((z(n) - a - rho * z(j) - step / 2) / sigma_e);
           else
               Pi(j,k) = normcdf((z(k) - a - rho * z(j) + step / 2) / sigma_e) - ...
                         normcdf((z(k) - a - rho * z(j) - step / 2) / sigma_e);
           end
       end
   end

end

% Schumaker is a shape preserving spline,in the sense that it preserves both
% concavity and monotonicity of the value function
% (shape=concavity/convexity+monotonicity)
% Since it calculates derivatives in each iteration, it comes at a cost of 
% lower speed
% However, here we have modified the code so as to calculate the
% coefficients pp only once at the beginning (lines 47-55) and not in every
% iteration, thus making the code much faster
