function output = withoutu (X, Y, Z, par, node, method, N, T)
% this function yields negative log likelihoods to be minimized

    gamma = par(1);
    beta0 = par(2);
    sigmab = par(3); % variance (not sd) of the normal in rc
if method == 1
    % Gaussian quadrature
    [rc, w] = qnwnorm(node, beta0, sigmab); % 3rd argument being the variance
    integrand = exp(1)*ones(length(rc), N);
    for i = 1:length(rc)
        betai = rc(i, 1);
        integrand(i,:) = prod(((logit(betai*X+gamma*Z)).^(Y)).*(ones(size(X))-logit(betai*X+gamma*Z)).^(ones(size(X))-Y), 1);
    end
    
    output = -sum(log((w.')*integrand)); % we get the negative log-likelihood
    clear integrand;
    
elseif method == 2
    % MC
    hNorm = haltonNormShuffle(node, 1, 6);
    rc = repmat(beta0, node, 1) + sigmab * hNorm';
    integrand = exp(1)*ones(node, N);
    for i = 1:node
        betai = rc(i, 1);
        integrand(i, :) = prod(((logit(betai*X+gamma*Z)).^(Y)).*(ones(size(X))-logit(betai*X+gamma*Z)).^(ones(size(X))-Y), 1);
    end
    
    output = -sum(log(mean(integrand, 1))); % we get the negative log-likelihood
    clear integrand;
end % end for if
end % end for function
    
        
        
