
"Hartree formula which returns the asymptotic coefficient C_κl of the atom."
function hartree_asymp_coeff(Z,Ip,l)
    if Z/sqrt(2*Ip) <= l
        @warn "[hartree_asymp_coeff] Z/κ<=l, unable to determine the asymptotic coefficient using Hartree formula. Using the l=0 case instead."
        l = 0
    end
    n = Z/sqrt(2*Ip)
    return 2^(n-1) / sqrt(n*gamma(n+l+1)*gamma(n-l))
end
