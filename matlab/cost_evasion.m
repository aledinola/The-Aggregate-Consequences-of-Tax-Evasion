function F = cost_evasion(phi,cc0,cc1,cc2)

%{
DESCRIPTION: 
This function returns the cost of evasion. It depends generally on three
parameters, c0,c1,c2 but in the published version of the paper we use only
c0, called as \kappa).
%}

if (phi>0) 
    F = cc0+cc1*phi+cc2*phi.^2;
else
    F = 0;
end
    
    
    
end %end function "cost_evasion"