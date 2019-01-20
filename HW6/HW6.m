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
        inv(inv<0) = 0; % avoid generating a imaginary number
        inv = kron(ones(1, Z), inv);
        prof = prof - 0.2 * inv .^ (1.5); % subtract inv costs from the gross profits
        [vnew(i,:),decision(i,:)]=max(prof + delta * EV);

    end
   dif=norm(vnew-v)/norm(vnew);
   disp(dif)
   v=vnew;
end

%% 
plot(k, v(:, 8), k, v(:, 11), k, v(:, 14))
title('Value function')
xlabel('Stock of lumber')
ylabel('Value')

%% Question 4

plot(grid, decision(1,:), grid, decision(101,:), grid, decision(201,:), grid, decision(301,:), grid, decision(401,:), grid, decision(501,:), grid, decision(601,:),grid, decision(701,:),grid, decision(801,:), grid, decision(901,:))

title('Next period stock vs lumber prices')
xlabel('Lumber price')
ylabel('Next period stock of lumber')

%% Question 5







