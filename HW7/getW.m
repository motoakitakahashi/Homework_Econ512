function output = getW(Vin, trans)
global L rho eta kappa l v delta beta lambda c CRIT Omega;

cont = ones(L, L, 3);
cont1 = trans(:, :, 1) * (Vin * (trans(:, :, 1)'));
cont2 = trans(:, :, 2) * (Vin * (trans(:, :, 1)'));
cont3 = trans(:, :, 1) * (Vin * (trans(:, :, 2)'));
cont(:, :, 1) = cont1;
cont(:, :, 2) = cont2;
cont(:, :, 3) = cont3;

output = cont;
end