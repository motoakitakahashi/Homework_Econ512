function p = sim_pi2(x)
x2 = sqrt(1-x.^2);
p = 4*sum(x2)/length(x);
end