% HW2 for Econ 512 Empirical Methods
% Motoaki Takahashi

clear

%% Question 2
disp('Question 2')

% Define the quality vector (2 by 1)
q = [2; 2];

% Refer to the function demand, demand.m, which returns the 2 by 1 demand vector for
% 2 by 1 price vector and 2 by 1 quality vector

% eqn.m defines the "left-hand side" of the equation of interest.
% That is, we want to get the value p that solves eqn(p, q)=0 for given
% quality q.

eqn([1;1],q)


