function output = withu (X, Y, Z, par, node, N, T)
% this function yields negative log likelihoods to be minimized

    gamma = par(1);
    beta0 = par(2);
    sigmab = par(3); % variance (not sd) of the normal in rc
    u0 = par(4);
    sigmau = par(5);
    sigmabu = par(6);
    mu = [beta0; u0]; % mean vector
    Sigma = [sigmab, sigmabu; sigmabu, sigmau]; % variance matrix
    %Cholesky decomposition
    U = chol(Sigma); %upper triangular matrix

    % MC
    hNorm = haltonNormShuffle(node, 2, 6); % the last argument being a seed
    rc = repmat(mu, 1, node) + U' * hNorm;
    integrand = exp(1)*ones(node, N);
    for i = 1:node
        betai = rc(1, i);
        ui = rc(2, i);
        integrand(i, :) = prod(((logit(betai*X+gamma*Z+repmat(ui, T, N))).^(Y)).*(ones(size(X))-logit(betai*X+gamma*Z+repmat(ui, T, N))).^(ones(size(Y))-Y), 1);
    end
    
    output = -sum(log(mean(integrand, 1))); % we get the negative log-likelihood
    clear integrand;

end % end for function