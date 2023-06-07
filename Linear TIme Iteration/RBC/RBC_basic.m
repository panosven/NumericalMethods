%-----------------------------------------------------------------------
% Βasic RBC model (McKay's website)
% Solution Method: Pontus Rendahl's (2017) Linear Time Iteration (LTI) method
% Panagiotis Veneris, U of Liverpool
% Date: 26/10/2020
%------------------------------------------------------------------------

close all;
clear all;
clc;

% Define Parameters
alppha = 0.36;    % capital share
gamma  = 1;       % Risk aversion
delta  = 0.025;   % Depreciation rate
rho_z  = 0.95;    % shock persistence
beta   = 0.97;    % discount factor
it     = 0;       % iteration counter

% Define Steady State (only of the variables that appear in SS form in the loglin.system equations) 
R = 1/beta;

% Define Coefficients (of variables that appear in the loglinear system of equations)
coeff0 = (R-1+delta)/R;
coeff1 = (1/gamma) * coeff0;
coeff2 = ((alppha-1)/gamma) * coeff0;
coeff3 = R;
coeff4 = (R-1+delta)/alppha;
coeff5 = (R-1+(delta*(1-alppha)))/alppha;

%---------------------------------------------------
% Matrix form of Linear RE model
% model:A*E_t{X_t+1}+B*X_{t}+C*X_{t-1}+W*u_{t}=0, (instead of ε in the notes
% here I use W, and instead of ε_{t} I use u_{t})
% X = [ c k z]'
%---------------------------------------------------
A = [ 1 0 -coeff1;
    0 0 0;
    0 0 0];

B = [ -1 -coeff2 0;
    coeff5 1 -coeff4;
    0 0 1];

C = [ 0 0 0;
    0 -coeff3 0;
    0 0 -rho_z];

W = [ 0; 0; -1];    % for negative prod.shock remove (-)

%-----------------------------------------------------------------------------------------------
% Solver for Linear Time Iteration (Rendahl 2017)

% Algorithm Description:
% Step1 : Guess a value P_{0}
% Step 2: Calculate P_{1} from P_{1} = -[A*P_{0}+B]^(-1) * C
% Step 3: If P_{1} solves (approxim.) the matrix quadr.eq.: A*P_{1}^2 + B*P_{1}+C within some
% tolerance then stop--> P_{1} is the solvent
% Step 4: Otherwise iterate again to get P_{2} and follow the same steps
%-----------------------------------------------------------------------------------------------

Tol = 1e-13;        % tolerance level for time iteration algorithm, default: 1e-13

metric = 1;         
P = 0;              % initial guess for P

% Run the algorithm until convergence
% Find  x_t = P x_{t-1}+Q u_t, where X_t is the policy function (ex how C responds to a 1% techn.shock (u_t)).

tic
while metric>Tol
    it = it+1;   % next iteration
    P      = -(A*P+B)\C;                                  
    metric = max(max(abs(A*P*P+B*P+C)));                                                                                 
    Q      = -inv(A*P+B)*W;      % could be outside the loop as well           
    fprintf('iteration = %d, metric = %e\n', it, metric);                
end
toc
% Check for stable equilibrium (all eigenvalues of solvent P must be < 1)
eig_P = max(abs(eig(P)));
if eig_P < 1 
    disp('Stable solution exists')   % all eigenvalues must be < 1 in absolute value
end

% IRFs
nresp = 80;           % horizon of responses (IRFs)
x = zeros(3, nresp);  % creates a 3x80 matrix of x's, that is, the policy functions of consumption (1), capital (2) and productivity (3), for 80 periods                     
u = zeros(1, nresp);  % creates a 1x80 matrix of u's (errors/disturbances)

% Initialize u   
u(:,1)= 1;  
               
for t=2:nresp                         
    x(:,t) = P*x(:,t-1)+ Q*u(:,t-1);                                      
end                                                                      

% How to generate the IRFs of the variables not present in the reduced system (y,R,k,prod.) 
c_irf = x(1,2:end);             % IRF for Consumption
k_irf = x(2,2:end);             % IRF for Capital
z_irf = x(3,2:end);             % IRF for Productivity                                    
y_irf = z_irf + alppha*k_irf;   % IRF of Output (y) (check equil.condition for y in notes)
R_irf = coeff0*(y_irf - k_irf); % IRf of Real Interest Rate (R) (check equil.cond. for r in notes)

% Plotting
figure;
subplot(3,2,1); plot(c_irf, 'b','Linewidth',2);                        
title('Consumption');
grid on;

subplot(3,2,2); plot(y_irf, 'b', 'Linewidth',2);
title('Output');
grid on;

subplot(3,2,3); plot(x(3,2:end), 'b', 'Linewidth',2);
title('Productivity');
grid on;

%figure;
subplot(3,2,4); plot(x(2,2:end), 'b', 'Linewidth',2);
title('Capital');
grid on;

subplot(3,2,5); area(R_irf,'FaceColor',[.6 1 1]);
plot(R_irf, 'b', 'Linewidth',2);
title('Real Interest Rate');
grid on;
