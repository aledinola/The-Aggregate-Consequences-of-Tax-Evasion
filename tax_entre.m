function T = tax_entre( income,Parameters )
%{
INPUTS:
income = taxable income
y = income/average(income)

OUTPUT:
T = total taxes paid. The average tax rate is T/income

NOTES:
Impose normalization (for flat tax is not needed)
y = income/(w0*uncmean_eps);
%}

%Unpack parameters

taxfunc = Parameters.taxfunc;
b_entre = Parameters.b_entre;
p_entre = Parameters.p_entre;
s_entre = Parameters.s_entre;

tau_s   = Parameters.tau_s;

switch taxfunc
    
    case 1 % HSV
         T = y - lambda_entre*y.^(1-tau_entre) + tau_s;
       
    case 2 %'gouveia'
        
        T = income*(b_entre - b_entre*(s_entre*income^p_entre + 1)^(-1/p_entre)) + tau_s;
      
    case 3 % flat
        
         T = taxrate_entre*income + tau_s;
end


end %END FUNCTION "tax_entre"