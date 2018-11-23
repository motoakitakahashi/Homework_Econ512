% Motoaki Takahashi
% HW5 for Econ 512 Empirical Method

clear
delete HW5log.txt
diary('HW5log.txt')
diary on
load('hw5.mat')
addpath('../CEtools/');
X = data.X;
Y = data.Y;
Z = data.Z;
N = 100;
T = 20;

%% Q1
disp('Question 1')
par = [0, 0.1, 1];
-withoutu(X, Y, Z, par, 20, 1, N, T)

%% Q2
disp('Question 2')
-withoutu(X, Y, Z, par, 100, 2, N, T)

%% Q3
disp('Question 3')
disp('Gaussian Quadrature')
% restrict the arguments to only par
withoutu_min = @(par) withoutu(X, Y, Z, par, 20, 1, N, T);
par = ones(1,3); % When I started with pi*ones(1,3), I didn't get the result.
[x, fval] = fminsearch(withoutu_min, par);
disp('The minimizer is')
disp(x)
disp('The value of the negative log-likelihood is')
disp(fval)

disp('Monte Carlo')
% restrict the arguments to only par
withoutu_min = @(par) withoutu(X, Y, Z, par, 100, 2, N, T);
par = ones(1,3);
[x, fval] = fminsearch(withoutu_min, par);
disp('The minimizer is')
disp(x)
disp('The value of the negative log-likelihood is')
disp(fval)


%% Q4

X = data.X;
Y = data.Y;
Z = data.Z;
N = 100;
T = 20;


disp('Question 4')
% restrict the arguments to only par
withu_min = @(par) withu(X, Y, Z, par, 100, N, T);
par = [ 1,1,1,1,1,0.3 ]; %cholesky decomposition needs a pd matrix
[x, fval] = fminsearch(withu_min, par);
disp('The minimizer is')
disp('   gamma      beta0    sigmab    u0    sigmau   sigmaub')
disp(x)
disp('The value of the negative log-likelihood is')
disp(fval)

diary off
