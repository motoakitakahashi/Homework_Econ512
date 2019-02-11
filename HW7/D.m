function output = D(v, p1, p2)
output = exp(v-p1) / (1 + exp(v-p1) + exp(v-p2));
end
