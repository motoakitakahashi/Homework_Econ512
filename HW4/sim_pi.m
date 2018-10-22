function p = sim_pi(x)
x = x.^2;
x2 = sum(x, 2);
x2 = x2(x2<=1);
p = 4*length(x2)/size(x, 1);
end