function [l,grad] = loglfgrad(beta)
load('hw3.mat');
l = -(-sum(exp(X*beta))+sum(y.*(X*beta))-sum(log(factorial(y))));
grad = -(-transpose(X)*exp(X*beta)+transpose(X)*y);
end
