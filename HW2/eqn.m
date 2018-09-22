function output = eqn(p,q)
output = ones(2, 1)-diag(demand(p,q).*(ones(2,1)-demand(p,q)))*p
end