function output = logit(x)
output = (ones(size(x))+exp(-x)).^(-1);
end