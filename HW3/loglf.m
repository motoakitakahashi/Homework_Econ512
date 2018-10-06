function l = loglf(x, y, beta)
l = -sum(exp(x*beta))+sum(y.*(x*beta))-sum(log(factorial(y)));
end
