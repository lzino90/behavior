function [z,zz,zzz,zzzz,u,r] = SIS_game_noSoc(n,m,a, lambda, mu,beta,sigma,c,gamma,u_func, r_func,T, init)
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
if length(init)==1
    x=zeros(1,n);
    x(randperm(n,round(n*init)))=1;
else
    x=init;
end
y=zeros(1,n);   %Initial self-protection is all zero
if length(a)==1
    a=a*ones(1,n);
end
z=zeros(1,T+1);     %The prevalence of the infected
zz=zeros(1,T+1);    %Proportion of individuals adopting self-protective behaviour
zzz = zeros(1,T+1); %Accumulation of the cost
zzzz = zeros(1,T+1); %Accumulation of the frustration
z(1)=sum(x)/n;
zz(1)=sum(y)/n;
f=zeros(1,n);   %Frustration function
u = zeros(1,T+1);
for t=2:T+1
    %% Decision making
    f=f*gamma+gamma*c*y;   %Accumulation of the frustration
    r(t-1) = risk_function(z(t-1),r_func);   %Risk perception
    u(t-1) = control_function(z,u_func,t,u);
    pi0=-u(t-1);                %Payoff for not self-protecting
    pi1=r(t-1)-c-f;         %Payoff for adopting self-protection
    pp=(exp(beta*pi1))./(exp(beta*pi1)+exp(beta*pi0));  %Probability of self-protection
    y=double(rand(1,n)<pp);        %Self-protection state update
    
    %% Epidemic
    W=ADN(n,m,a);          %Generate contact network for given time-step
    ip=(1-(1-lambda).^(W*x'))';    %Probability of individuals being infected
    %Infected state transition. First term captures transition to I, and
    %second term recovery to S.
    x=x+(rand(1,n)<=ip.*(1-x).*(1-sigma*y))-(rand(1,n)<mu*x);
    z(t)=sum(x)/n;
    zz(t)=sum(y)/n;
%         zzz(t) = (zzz(t-1)+c*sum(y))./(n*(t-1));
    zzz(t) = zzz(t-1)+c*sum(y);
    zzzz(t) = sum(c+f)/n;
end


%% Plotting
figure
plot(0:T,zzz)
xlabel('Time Step, t')
ylabel('Accumulative Cost')

% figure
% plot(0:T,zzzz)
% xlabel('Time Step, t')
% ylabel('Mean Individual Frustration')

figure(n)
plot(0:T,z)
hold on
plot(0:T,zz)
% plot(0:T,zzzz)
legend('Prevalence of Disease','Fraction of Individuals Self-Isolating')
xlabel('Time Step, t')
% title('n=5000')
mean(zz)

end