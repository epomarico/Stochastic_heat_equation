% STOCHASTIC HEAT EQUATION
%
% This program solves numerically with an implicit Euler method
% the stochastic heat equation:
% 
% \partial_t u_t = alpha*\Delta u + sigma*\partial_t W (x,t) 
% + beta sin(u)
% with alpha, beta, sigma constants,
%      Dirichlet boundary conditions u(0,t) = u(1,t) = 0,
%      and u(x,0)=x(1?x) initial conditions 
%
% and produces a trajectory in space-time in the interval [0,1]


clear all

%Initialize random number generator
randn('state',100)

% Parameters of the equation
alpha=1; sigma=0.1; beta=1;
% Time parameters
tmax=1; Nt=100;
% Space parameters
L=1; Nx=100;
% Compute mesh in space and time 
dx = L/(Nx-1); dt = tmax/(Nt-1);

% Vectors x and t, and matrix u of size Nx x Nt 
x = linspace(0,L,Nx)';
t = linspace(0,tmax,Nt);
u = zeros(Nx,Nt); 

% Initial conditions
u(:,1) = x.*(1-x);    
% Dirichelet boundary conditions
u0 = 0;   uNx = 0;

% Coefficients of the system of equations
a = dt*(-alpha/dx^2)*ones(Nx,1); 
b = ones(Nx,1) - 2*a;
c = a;

% Dirichlet boundary conditions
b(1) = 1; c(1) = 0;       % to have  u_0^m = d(1) = 0. The value is put to zero in the for loop
a(end) = 0; b(end) = 1;   % to have  u_(Nx)^m = d(end) = 0.  The value is put to zero in the for loop

% Factorize tridiagonal matrix A into L and U matrices, characterized by
% vectors e and f that are returned by the LU_tridiag function
[e,f] = LU_tridiag(a,b,c);


for m=2:Nt  % Loop over time steps. A vector u is formed at each time step
    
  dW=sqrt(dt)*randn(Nx,1); % Brownian increments
  d = u(:,m-1) + (sigma*dW)/(sqrt(dx)) + beta * sin(u(:,m-1)) * dt; % d coefficients for the time-step in consideration
  d(1) = u0;  d(end) = uNx;     % boundary conditions
  
  u(:,m) = solve_Aud(d,a,e,f); % solves the system of equations and put the solution into a 
                               % a single column of the matrix u,
                               % correspondint to one time-step
end




% Plotting trajectory in space-time
close all
str = {'\fontsize{16}Stochastic heat equation','\alpha = 1, \sigma = 0.1, \beta = 10'};


figure(1);
imagesc(u(:,1:100))
dim = [0.35 0.65 0.26 0.25];
annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','Color','white');
xlabel('\fontsize{16} t')
ylabel('\fontsize{16} x')
colorbar;
