close all 
clear variables
clc
tic

global T_extinct T_cost T1

T_extinct = inf;
T_cost = inf;
T1 = 35;
%% Network Parameter setup

% n=5000;   %Number of individuals in the population
% A=WS(n,8,.125);   %Adjacency matrix for the influence layer

n=10000;   %Number of individuals in the population
k = 8;
% A=zeros(n,n);
A=WS(n,k,1/k);   %Adjacency matrix for the influence layer
%  A=AB(n,k/2);

m=13;      %Number of edges per activity

% a=.1;     %Activity rate constant

%Activity rate power law
a_min=0.1; %min 
a_max=1; %max
a=(a_max+(1/a_min-a_max)*rand(1,n)).^(1/(-2.09+1));
a1=mean(a); %mean
a2=mean(a.^2); %second moment


%% Epidemic parameters
tau_E=6.4;
tau_I=5;
mu=1-exp(-1/tau_I);    %Recovery rate using Lancet paper
nu=1-exp(-1/tau_E);    %Transition from exposed to infected for SEIR model, using Lancet paper
R = 2.2;
lambda=R/(2*m*a1*tau_I);   %Infection probability
init=0.001;   %Initial fraction of infecteds 

%% Decision making setup
beta=6;      %Rationality of decision-making
c=.3;       %Base cost of social distancing for one time step
gamma =0.75; %Discount rate for social distancing of past steps, scalar (code in YYY_game would need to be updated if vector)
u_in=[5, 0.7];       %Vector [Control law, Control effort]
f=[2, 3];        %Vector [Risk perception function case, scaling constant]
T=400;      %Number of time steps


% thres = (exp(beta)-beta*exp(-beta*c))/(m*(a1+sqrt(a2))*(exp(beta)+(1-beta)*beta*exp(-beta*c)));
% thres_ac = lambda/mu;



%% Simulations

[z,zz,zzz,zzzz,R,u,expo,prev]= SEIR_game(n,A,m,a, lambda, nu, mu,beta,c,gamma,u_in, f,T, init);




% 28*u_in(2)+0.5*T1*u_in(2)
% (2(46.2-28*0.9))/0.9
toc