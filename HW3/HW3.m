% Motoaki Takahashi
% HW3 for Econ 512 Empirical Method

clear
diary hw3.out

% x is the n by 6 matrix (explanatory), y is the n by 1 vector (explained)

%% Question 1

load('hw3.mat');
% we have X and y

% write the negative likelihood function as a function of a parameter
loglf_beta = @(beta) -loglf(X,y,beta)

% set the initial guess
beta = [log(mean(y)); zeros(5,1)];

[est_nm, valnm] = fminsearch(loglf_beta, beta)

%% Question 2

% I use the BFGS method.
% the function loglfgrad contains the objective function and the gradient

options = optimoptions('fminunc','Algorithm','quasi-newton',...
          'SpecifyObjectiveGradient',true, 'Display','iter', 'MaxFunctionEvaluations', 30000, 'MaxIterations', 10000);
[est_BFGS, valBFGS] = fminunc('loglfgrad', beta, options)
disp(est_BFGS);

% I calculate the BFGS outcome without the analytical gradient as well

[est_BFGS_wog, val_BFGS_wog] = fminunc(loglf_beta, beta)

% See the difference b/w the outcomes with and without the gradient
est_BFGS-est_BFGS_wog

%% Question 3
nls_beta=@(beta) nls(X,y,beta);
options1 = optimoptions(@lsqnonlin, 'MaxFunctionEvaluations', 30000, 'MaxIterations', 10000)
nls_com = lsqnonlin(nls_beta, beta, -Inf, +Inf, options1)
% why is it so different from the other estimates? it should have rung
% alarms in you
%% Question 4
nls_nm = fminsearch(nls_beta, beta)

%% Question 5

% See what happens if we move the initial guess for the 3rd element (the coefficient for the years married) from
% -10 to 10. The initial guess for beta_0 is log(mean(y)). The others are kept 0.

grid = [-10:0.5:10];

beta_mat = zeros(6, length(grid));

% the initial guess for beta_0 is log(mean(y))
beta_mat(1, [1:length(grid)]) = log(mean(y))*ones(1, length(grid));

beta_mat(3,[1:length(grid)]) = grid;
% 
% %Nelder-Mead for ML
nm_mat = zeros(6, length(grid));
tic
for n=1:length(grid);
    nm_mat([1:6],n) =  fminsearch(loglf_beta, beta_mat([1:6],n));
end
toc
% 
% %BFGS for ML
BFGS_mat = zeros(6, length(grid));
tic
for n=1:length(grid);
    BFGS_mat([1:6],n) =  fminunc('loglfgrad', beta_mat([1:6],n), options);
end
toc
% % BFGS for ML without the analytical gradient
BFGS_wog_mat =  zeros(6, length(grid));
tic
for n=1:length(grid);
    BFGS_wog_mat([1:6],n) =  fminunc(loglf_beta, beta_mat([1:6],n));
end
toc



% % %lsqnonlin
nls1_mat = zeros(6, length(grid));
tic
for n=1:length(grid);
    nls1_mat([1:6],n) =  lsqnonlin(nls_beta, beta_mat([1:6],n), -Inf, +Inf, options1);
end
toc
% % 
% % %Nelder-Mead for NLS
nls2_mat = zeros(6, length(grid));
tic
for n=1:length(grid);
    nls2_mat([1:6],n) =  fminsearch(nls_beta, beta_mat([1:6],n));
end
toc

% Draw a figure
% I plot the estimates of beta_3 for each initial guess 
plot(grid, nm_mat(3, [1:length(grid)]), grid, BFGS_mat(3, [1:length(grid)]), grid, BFGS_wog_mat(3, [1:length(grid)]), grid, nls1_mat(3, [1:length(grid)]), grid, nls2_mat(3, [1:length(grid)]))

legend('Nelder-Mead for ML', 'BFGS with the grad for ML', 'BFGS without the grad for ML', 'Built-in for NLS', 'Nelder-Mead for NLS')

% Since the nelder-mead for NLS is exceptionally bad, so plot the others
% separately
plot(grid, nm_mat(3, [1:length(grid)]), grid, BFGS_mat(3, [1:length(grid)]), grid, BFGS_wog_mat(3, [1:length(grid)]), grid, nls1_mat(3, [1:length(grid)]))

legend('Nelder-Mead for ML', 'BFGS with the grad for ML', 'BFGS without the grad for ML', 'Built-in for NLS')

% Plot the first three

plot(grid, nm_mat(3, [1:length(grid)]), grid, BFGS_mat(3, [1:length(grid)]), grid, BFGS_wog_mat(3, [1:length(grid)]))

legend('Nelder-Mead for ML', 'BFGS with the grad for ML', 'BFGS without the grad for ML')

% Plot the first two

plot(grid, nm_mat(3, [1:length(grid)]), grid, BFGS_mat(3, [1:length(grid)]))

legend('Nelder-Mead for ML', 'BFGS with the grad for ML')

diary off