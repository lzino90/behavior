function [z,zz,zzz,zzzz,R,u,expo,prev] = SEIR_game_NI(n,A,m,a, lambda, nu, mu,beta,c,gamma,u_func, r_func,T, init)
%n: number of nodes
%A: adjacency matrix influence layer
%m: number of social interactions
%a: activity vector (if uniform, a single entry is ok)
%lambda: per-contact infection probability
%mu: recovery probability
%beta: rationality
%c: cost for self-protective behavior
%gamma: forgetting factor
%u_func: control function (see control_function.m)
%r_func: risk perception function (see risk_function.m)
%T: simulation time-horizon
%init: initial fraction of infected nodes, or initial state

global T_extinct T_cost Tcount

if length(init)==1
    x=zeros(1,n);
    x(randperm(n,round(n*init)))=1;
else
    x=init;
end
y=zeros(1,n);
if length(a)==1
    a=a*ones(1,n);
end
A=A*diag(sum(A,2).^(-1));   %Creates a column-stochastic adjacency matrix
z=zeros(1,T+1);     %The prevalence of the infected
expo=zeros(1,T+1);     %The prevalence of the exposed
prev = zeros(1,T+1);   %The prevalence of all individuals with disease (E+I)
R=zeros(1,T+1);     %Proportion of recovered individuals
zz=zeros(1,T+1);    %Proportion of individuals adopting self-protective behaviour
zzz = zeros(1,T+1); %Accumulation of the cost
zzzz = zeros(1,T+1); %Accumulation of the frustration
z(1)=sum(x)/n;
zz(1)=sum(y)/n;
f=zeros(1,n);   %Frustration function
u = zeros(1,T+1);
cost = zeros(1,T+1);
Tcount=-Inf;
for t=2:T+1
    %% Decision making
    f=f*gamma+gamma*c*y;   %Accumulation of the frustration
    r(t-1) = risk_function(z(t-1),r_func);   %Risk perception
    u(t-1) = control_function(z,u_func,t,u);
   %  pi0=(1-y)*A'-u(t-1);                %Payoff for not self-protecting
    % pi1=y*A'+r(t-1)-c-f;         %Payoff for adopting self-protection
        pi0=-u(t-1);                %Payoff for not self-protecting
    pi1=r(t-1)-c-f;         %Payoff for adopting self-protection
    pp=(exp(beta*pi1))./(exp(beta*pi1)+exp(beta*pi0));  %Probability of self-protection
    y=double(rand(1,n)<pp);        %Self-protection state update
    
    %% Epidemic
    
    W=ADN(n,m,a);
    ip=(1-(1-lambda).^(W*((x==2)')))';
    
    x=x+(rand(1,n)<=ip.*(x==0).*(1-y))+(rand(1,n)<nu*(x==1))+(rand(1,n)<mu*(x==2));
    z(t)=sum((x==2))/n;
    expo(t) = sum((x==1))/n;
    prev(t) = z(t)+expo(t);
    
    R(t)=sum((x==3))/n;
    zz(t)=sum(y)/n;
    
    zzz(t) = zzz(t-1)+c*sum(y);
    cost(t) = cost(t-1)+u(t-1)*sum(1-y);
    zzzz(t) = sum(c+f)/n;
    
    if prev(t) == 0
        T_extinct = min([t,T_extinct]);
        %        zzz(t) = zzz(t-1);
    end
    
    %     if (prev(t) == 0) && (zz(t) == 0)
    %         T_cost = min([t,T_cost]);
    %         zzz(t) = zzz(T_cost);
    %     end
end
figure
plot(0:T,zzz./n)
hold on
% plot(0:T,cost)
xlabel('Time Step, t')
ylabel('Accumulative Cost')
% legend('Cost')

figure
plot(0:T,u)
xlabel('Time Step, t')
ylabel('Control')



figure
hold on
yyaxis left
plot(0:T,prev)
plot(0:T,expo)
plot(0:T,z)
plot(0:T,R)
ylabel('Fraction of the Population')
ylim([0 max(R)+0.2])

yyaxis right
plot(0:T,zz)
ylim([0 1])
legend('Prevalence of Disease','Exposed','Infected','Recovered','Fraction Self-Isolating')
% plot(0:T,zzzz)
% legend('Prevalence of Disease','Fraction of Individuals Self-Isolating','Mean Individual Frustration')
title(['n = ',num2str(n)])
xlabel('Time Step, t')
end