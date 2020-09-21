function pk = prob_audit(phi,profit,theta,k,n,pflag,pn_1,pn_2,pn_3,le,gamma,vi)
%------------------------- LEGEND ----------------------------------------%
%{
This function gives the prob of auditing. We allow for a number of
specifications, to be selected with "pflag".
INPUTS:
    phi: fraction of misreported business income
    profit: business income
    theta,k,n: production inputs
    pflag: integer variable to specify the argument of p(.)
    pn_1,pn_2,le,gamma,vi: extra parameters that need to be communicated

OUTPUT
    pk: probability of audit (a scalar)

NOTES:
    Updated by Alessandro Di Nola on September 16, 2020. 
%}

switch pflag
    
    case (1) %misreported income
        temp = phi*profit;
        
    case (2) %f(theta,k,n)
        temp = theta*prodfun(k,le,n,gamma,vi);
        
    case (3) %f(k,n)
        temp = prodfun(k,le,n,gamma,vi);
        
    case (4) %k (as it was in the previous version of the paper)
        temp = k;
        
    case (5) %n
        temp = n;
        
end

pk = pn_3 + (1-pn_3)./(1 + pn_1*exp(-pn_2*temp));

end %END FUNCTION "prob_audit"
