function output = getp(p, W)
global L rho eta kappa l v delta beta lambda c CRIT Omega;
%this function is from the FOC
cmat = kron(c, ones(1,L));
Done = D(p, p');
Dout = ones(L, L) - D(p, p') - D(p', p); % demand for the outside good
Dopp = D(p', p); % firm 2's demand (prob)

W1 = W(:,:,1);
W2 = W(:,:,2);
W3 = W(:,:,3);

output = cmat + ...
    (ones(L, L) - beta * W2 + beta * (Dout .* W1 + Done .* W2 + Dopp .* W3)) ./ (ones(L, L) - Done);
end