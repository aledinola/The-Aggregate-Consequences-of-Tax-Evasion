function F = prodfun(k,le,n,gamma,vi)
%{
DESCRIPTION: 
 This function returns f(k,le,n), where k,le,n
 are the production inputs used by the entrepreneur.
 
INPUTS
 k, le, n: capital, labor supply, labor hirings
 gamma, vi: additional parameters
 REMARK: production fun is defined without "theta"
        See also fortran function in "mod_baselib.m"

OUTPUTS
 F: production function f(k,le,n) (without ability "theta")

NOTES:
 Updated by Alessandro Di Nola on September 18, 2020.
%}

F = (k.^gamma.*(le+n).^(1-gamma)).^vi;

end %end function "prodfun"

