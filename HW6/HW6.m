% Motoaki Takahashi
% HW6 Econ 512

clear all
delta = 0.95;
p0 = 0.5;
rho = 0.5;
N = 1000;
k0 = 100; % initial stock of lumber
k = (k0/N):k0/N:k0;

sigmau = 0.1;

%% Question 2
Z = 21; % number of grid points
% Z = 5 % for coarse grid

 
 [prob,grid]=tauchen(Z,p0,rho,sigmau);
disp(['The dimensions of prob are ' num2str(size(prob)) ])
disp(['The dimensions of grid are ' num2str(size(grid)) ])

%% Question 3


v = zeros(N, Z); % initial guess for value function
decision = zeros(N,Z); % this will contain the firm's policy
newv = zeros(N,Z); % this will contain 

%% value function iteration

dif = 1;
tol = 1E-4;
while dif > tol
    EV = v * prob';
    for i = 1:N
        prof = kron(grid, k(i)*ones(N, 1)-k');
        prof(prof < 0) = -1E5; % punish a negative stock of lumber
        inv = k(i)*ones(N, 1)-k';
        inv(inv<0) = 0; % avoid generating an imaginary number
        inv = kron(ones(1, Z), inv);
        prof = prof - 0.2 * inv .^ (1.5); % subtract inv costs from the gross profits
        [vnew(i,:),decision(i,:)]=max(prof + delta * EV);

    end
   dif=norm(vnew-v)/norm(vnew);
   disp(dif)
   v=vnew;
end

%% 
plot(k, v(:, 8), k, v(:, 11), k, v(:, 14)) % for grid Z = 21
% plot(k, v(:, 2), k, v(:, 3), k, v(:, 4)) % for coarse grid Z = 5
title('Value function')
xlabel('Stock of lumber')
ylabel('Value')

%% Question 4

plot(grid, decision(1,:), grid, decision(101,:), grid, decision(201,:), grid, decision(301,:), grid, decision(401,:), grid, decision(501,:), grid, decision(601,:),grid, decision(701,:),grid, decision(801,:), grid, decision(901,:))

title('Next period stock vs lumber prices')
xlabel('Lumber price')
ylabel('Next period stock of lumber')

%% Question 5

% simulation
% let's run 1000 paths and numerically get the means and the confidence
% interval of stocks

p = 1; % the initial value of lumber
k0 = 100; % the initial stock of lumber

% we will get n = 1000 simulations
n = 1000;
numSteps = 21; % initial state is p=1 and 20 periods ahead

simu_price = zeros(n, numSteps); % this will contain the simulated paths of prices
simu_price(:, 1) = ((Z+1)/2) * ones(n, 1);

% make the cumulative version of the transition matrix
cumu = cumsum(prob, 2);

for s = 1:n
    for step = 2:21
        r = rand;
        simu_price(s, step) = find(cumu(simu_price(s, step-1),:) > r, 1);
    end
end

% Get the simulated paths of lumber stocks associated with the prices

simu_ind = zeros(n, numSteps);
simu_ind(:, 1) = N * ones(n, 1);
simu_cap = zeros(n, numSteps);
simu_cap(:, 1) = k0 * ones(n, 1);

for s = 1:n
    for step = 2:21
        simu_ind(s, step) = decision(simu_ind(s, step-1), simu_price(s, step-1));
        simu_cap(s, step) = k(simu_ind(s, step));
    end
end

meancap = mean(simu_cap);
secap = std(simu_cap) ./ sqrt(n);
ts = tinv([0.05  0.95], n-1); 

cicap = kron(meancap, ones(2, 1)) +  kron(ts', secap);
memori = 1:numSteps;
plot(memori, meancap, memori, cicap)

title('Mean and 90 percent CI for lumber stocks')
xlabel('Period')
ylabel('Stock of lumber')

%% Question 6

% for question 6, redo with Z = 5



