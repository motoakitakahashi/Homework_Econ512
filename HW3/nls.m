function n = nls(X,y,beta)
n=sum((y-exp(X*beta)).^2);
end