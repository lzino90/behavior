function [u] = control_function(z,func,t,p)
%CONTROL_FUNCTION Allows for different control functions
%INPUT:   func = Vector [The control function case, scaling constant]
%           z = prevalence of the disease
%           t = current time step
%           p = control input up to this point
%OUTPUT:  u = the control value

global Tcount T1

switch func (1)
    case 1   %Constant control
        u = func(2);
        
    case 2   %Constant if prevalence above certain threshold
        if z(t-1) > 0.01
            u = func(2);
        else
            u=0;
        end
        
    case 3   %Constant if prevalence above certain threshold at any given time, with reset
        
        if (z(t-1) <= 0.01) && (p(t-1) == 0)
            u = 0;
        else
            u = func(2);
        end
        
    case 4   %Constant if prevalence breaches certain threshold once
        
        if max(z) > 0.01 && z(t-1) ~= 0
            u = func(2);
        else
            u = 0;
        end
        
        
    case 5   %Constant after threshold is reached, and resets if after 14 days there is no new cases.
        T = 14;
        if z(t-1) >= 0.01
            u=func(2);
        else
            if t > 20
                if (sum(z(t-T-1:t-1)) == 0)
                    u = 0;
                else
                    u=func(2);
                end
            else
                u = 0;
            end
            
        end
        
        
    case 6   %Constant for T days after threshold is reached
        T = 28;
        if (t == 1) || (t==2)
            u = 0;
        elseif z(t-1) >= 0.01 && p(t-2)==0
            u= func(2);
            Tcount = t-1;
        elseif t-Tcount <=T+1
            u=func(2);
        else
            u=0;
        end
        
    case 7   %Constant after threshold is reached, then linear decay
        T = 28;
        if (t == 1) || (t==2)
            u = 0;
        elseif z(t-1) >= 0.01 && p(t-2)==0
            u= func(2);
            Tcount = t-1;
        elseif t-Tcount <=T+1
            u=func(2);
        else
            u=max(0,p(t-2)-func(2)/T1);
        end


  case 8   %Constant for T=56 days after threshold is reached
        T = 56;
        if (t == 1) || (t==2)
            u = 0;
        elseif z(t-1) >= 0.01 && p(t-2)==0
            u= func(2);
            Tcount = t-1;
        elseif t-Tcount <=T+1
            u=func(2);
        else
            u=0;
        end
end
end
