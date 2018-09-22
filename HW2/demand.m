function d = demand(p, q)
d=exp(q-p)/sum(1+sum(exp(q-p)))
end