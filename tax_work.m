function T = tax_work( income,Parameters )

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
b_work  = Parameters.b_work;
p_work  = Parameters.p_work;
s_work  = Parameters.s_work;

tau_s   = Parameters.tau_s;


switch taxfunc
    
    case 1 %'hsv' 
         T = y - lambda_work*y.^(1-tau_work) + tau_s;
       
    case 2 %'gouveia'
        
        T = income*(b_work - b_work*(s_work*income^p_work + 1)^(-1/p_work)) + tau_s;
        
    case 3 %'flat'
        
        T = taxrate_work*income + tau_s;
        
        
end

end