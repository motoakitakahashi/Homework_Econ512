function output = D(p1, p2)
% this is the demand for firm 1 (say)
% this function returns the L by L demand matrix for two L by L matrices of
% prices

global L rho eta kappa l v delta beta lambda c CRIT Omega;
vmat = v * ones(L, L);
output = exp(vmat - p1) ./ (ones(L, L) + exp(vmat - p1) + exp(vmat - p2));
end
