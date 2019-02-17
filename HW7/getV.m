function output = getV(np, p, W)
global L rho eta kappa l v delta beta lambda c CRIT Omega;
cmat = kron(c, ones(1,L));
Done = D(np, p');
Dout = ones(L, L) - D(np, p') - D(p', np); % demand for the outside good
Dopp = D(p', np); % firm 2's demand (prob)
W1 = W(:,:,1);
W2 = W(:,:,2);
W3 = W(:,:,3);

output = Done .* (np - cmat) + beta * (Dout .* W1 + Done .* W2 + Dopp .* W3);
end