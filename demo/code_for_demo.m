close all 
clear variables
clc
tic

%% Parameter setup

% n=5000;   %Number of individuals in the population

n=10000;   %Number of individuals in the population
 A=WS(n,8,.125);   %Adjacency matrix for the influence layer


m=1;      %Number of edges per activity

%a=.2;     %Activity rate constant

%Activity rate power law
a_min=0.1; %min 
a_max=1; %max
a=(a_max+(1/a_min-a_max)*rand(1,n)).^(1/(-2.2+1));
a1=mean(a); %mean
a2=mean(a.^2); %second moment

%% Epidemic parameters
mu=1-exp(-1/7.8571);    %Recovery rate using Lancet paper
R = 1.6;
lambda=R/(2*m*a1*7.8571);   %Infection probability


%% Decision making dynamics
beta=6;      %Rationality of decision-making
c=.3;       %Base cost of social distancing for one time step
gamma = 0; %Discount rate for social distancing of past steps, scalar (code in YYY_game would need to be updated if vector)
u_in=[4,0];       %Vector [Control law, Control effort]
f=[3, 3];        %Vector [Risk perception function case, scaling constant]
T=600;      %Number of time steps
init=0.001;   %Initial fraction of infecteds 
sigma=0.99; %effectiveness of self-protections

% thres = (exp(beta)-beta*exp(-beta*c))/(m*(a1+sqrt(a2))*(exp(beta)+(1-beta)*beta*exp(-beta*c)));
% thres_ac = lambda/mu;



%% Simulations
[z,zz,zzz,zzzz,u,r]=SIS_game(n,A,m,a, lambda, mu,beta,sigma,c,gamma,u_in, f,T, init);

toc