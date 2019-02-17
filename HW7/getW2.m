function output = getW2(Vin)
global L rho eta kappa l v delta beta lambda c CRIT Omega;

cont = ones(L, L, 3);

% this function does not explicitly use the transition matrices

V = [Vin(1,:); Vin; Vin(L,:)];
% V = [V(1,:) ; V ; V(L,:)]

% cont1 is the continuation value W when both firms don't make a sale
cont1 = (ones(L, L) - ((1-delta)* ones(L, L)) .^ (Omega)) .* V(1:L, :) + (((1-delta)* ones(L, L)) .^ (Omega) .* V(2:L+1, :));
cont1 = [cont1(:,1), cont1, cont1(:,L)];

cont1 = ((ones(L, L) - ((1-delta)* ones(L, L)) .^ (Omega')) .* cont1(:, 1:L)) + (((1-delta)* ones(L, L)) .^ (Omega') .* cont1(:, 2:L+1));

cont(:, :, 1) = cont1;


% cont2 is the continuation value W when only firm 1 makes a sale

cont2 = ((ones(L, L) - ((1-delta)* ones(L, L)) .^ (Omega)) .* V(2:L+1, :)) + (((1-delta)* ones(L, L)) .^ (Omega) .* V(3:L+2, :));
cont2 = [cont2(:,1), cont2, cont2(:,L)];

cont2 = ((ones(L, L) - ((1-delta)* ones(L, L)) .^ (Omega')) .* cont2(:, 1:L)) + (((1-delta)* ones(L, L)) .^ (Omega') .* cont2(:, 2:L+1));

cont(:, :, 2) = cont2;

% cont3 is the continuation value W when only firm 2 makes a sale
cont3 = ((ones(L, L) - ((1-delta)* ones(L, L)) .^ (Omega)) .* V(1:L, :)) + (((1-delta)* ones(L, L)) .^ (Omega) .* V(2:L+1, :));
cont3 = [cont3(:,1), cont3, cont3(:,L)];

cont3 = ((ones(L, L) - ((1-delta)* ones(L, L)) .^ (Omega')) .* cont3(:, 2:L+1)) + (((1-delta)* ones(L, L)) .^ (Omega') .* cont3(:, 3:L+2));

cont(:, :, 3) = cont3;

output = cont;

end