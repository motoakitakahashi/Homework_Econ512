function output = c(kappa, eta,  omega, l)
if omega<=l
    output = kappa * omega ^ (eta);
else
    output = kappa * l ^ (eta);
end
return output
end
