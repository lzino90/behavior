close all 
clear variables
clc
tic

global T_extinct T_cost

T_extinct = inf;
T_cost = inf;

%% Network Parameter setup

% n=5000;   %Number of individuals in the population
% A=WS(n,8,.125);   %Adjacency matrix for the influence layer

n=20000;   %Number of individuals in the population
k = 8;
A=WS(n,k,1/k);   %Adjacency matrix for the influence layer

%toc

m=13;      %Number of edges per activity

% a=.1;     %Activity rate constant

%Activity rate power law
a_min=0.1; %min 
a_max=1; %max
a=(a_max+(1/a_min-a_max)*rand(1,n)).^(1/(-2.09+1));
a1=mean(a); %mean
a2=mean(a.^2); %second moment
sigma=0.95;


%% Epidemic parameters
tau_I=4.1;
mu=1-exp(-1/tau_I);    %Recovery time
R = 2; 
lambda=R/(2*m*a1*tau_I);   %Infection probability
init=0.001;   %Initial fraction of infecteds 

%% Decision making setup
beta=6;      %Rationality of decision-making
c=.1;       %Base cost of social distancing for one time step
gamma = 0.9; %Discount rate for social distancing of past steps, scalar (code in YYY_game would need to be updated if vector)
u_in=[1, 0.4];       %Vector [Control law, Control effort]
f=[2, 3];        %Vector [Risk perception function case, scaling constant]
T=200;      %Number of time steps



%% Simulations

[z,zz,zzz,zzzz,R,u]= SIR_gameN(n,A,m,a, lambda, mu,beta,sigma,c,gamma,u_in, f,T, init);

toc