% Stohastic NGM: VFI + Cubic spline
% Panagiotis Veneris, U of Liverpool
% 23/1/2021 

clear;
close all;
clc;

%Denote global variables
global beta alpha kGrid k_init V n delta z_init Pi z

alpha   = 0.3;
beta    = 0.96;     % discount factor
tol_pol = 1e-6;     % tolerance level
maxit   = 500;      % max # of iterations
dif     = 0.1;
it      = 0;        % iteration counter
n       = 9;        % number of gridpoints for productivity
sigma_e = 0.007;    % std. deviation of stohastic process
rho     = 0.95;     % persistence of stohastic process
mu      = 0;        % mean of stohastic process
w       = 3 ;       % how many std. dev. above/below mean
delta   = 0.1;      % depreciation rate

% Grid for state variable (k)
kGridnum = 30;        % number of gridpoints for capital
kmin     = 1;
kmax     = 5;
kGrid    = linspace(kmin,kmax,kGridnum)';  % equally-spaced grid for k (30x1 matrix)

% Grid for state variable (e)
[Pi,z] = tauchen(n,mu,rho,sigma_e,w);

% Initial guess for value function
V = zeros(kGridnum,n);   % V_0(k) (30x1 matrix)

% Initialize auxiliary matrices
V1 = zeros(kGridnum,n);
K11 = zeros(kGridnum,n);
k_init = zeros(kGridnum,n);
z_init = zeros(kGridnum,n);

tic
while dif>tol_pol && it<maxit   
    it = it+1;
    
    for j=1:n
        for i=1:kGridnum
            k_init = kGrid(i,1);
            z_init = z(j,1);
            K1 = fminbnd(@stoh_spline,kmin,kmax);
            V1(i,j) = -stoh_spline(K1);
            K11(i,j) = K1;
        end
    end
    
    dif = norm(V1-V);  
    V   = V1;   
    fprintf('iteration = %d, dif = %e\n', it, dif);          
end
toc

% Plots
figure
 hold on;
 plot(kGrid,V1(:,1),'-c','Linewidth',2);
 plot(kGrid,V1(:,n),'b*-','Linewidth',2);
 hold off;
 xlabel('k');
 ylabel('V(k)');
 legend('n=1','n=9', 'Location', 'best');
 title('Value Function');
 grid on;

figure
 hold on;
 plot(kGrid,K11(:,1),'-c','Linewidth',2);
 plot(kGrid,K11(:,n),'b*-','Linewidth',2)
 hold off;
 xlabel('k');
 ylabel('g(k)');
 legend('n=1','n=9', 'Location', 'best');
 title('Capital Policy Funtion');
 grid on;
 
% figure
%  hold on;
%  plot(kGrid,c,'-c','Linewidth',2);
%  plot(kGrid,c_actual,'b*-');
%  xlabel('k');
%  ylabel('c');
%  legend('VFI','True');
%  title('Consumption Policy Function');
%  grid on;
%   

function value = stoh_spline(kprime)

 global alpha beta delta k_init z_init kGrid Pi n V 

  g = interp1(kGrid,V,kprime,'spline');  % can be outside of the loop
     
  c = exp(z_init)*(k_init^alpha)+(1-delta)*k_init-kprime; % can be outside of the loop
 
 %for j=1:n    
%      g = interp1(kGrid,V,kprime,'spline');
%      
%      c = exp(z_init)*(k_init^alpha)+(1-delta)*k_init-kprime;    
   if  c<0
        value = -888888888888888888-800*abs(c);     % keeps it from going negative
   else
        for j=1:n
        value = log(c) + beta*(g*Pi(j,:)');
        end
   end
   value = -value;
%end
end
       
% Discretize the aggregate shock (z) using Tauchen method
function [Pi,z] = tauchen(n,mu,rho,sigma_e,w)

   z    = zeros(n,1);      % nodes for stohastic process (possible different states)
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
