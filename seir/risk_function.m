function [r] = risk_function(z,func)
%RISK_FUNCTION Allows for different risk perception functions
%INPUT:   func = Vector [The risk function case, scaling constant]
%           z = prevalence of the disease
%OUTPUT:  r = the risk perception value
switch func (1)
    case 1   %Linear risk
        r = z*func(2);
    case 2   %Under reaction
        r = func(2)*(z^2);
    case 3   %Over reaction
        r = func(2)*(z^0.5);
end

end

