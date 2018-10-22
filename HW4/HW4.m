% Motoaki Takahashi
% HW4 for Econ 512 Empirical Method

clear
%% Question 1

% I draw 1000 points in the unit square from Halton sequence.
n = 100^2;
h = haltonseq(n, 2);
hsq = h.^2;
hsq = sum(hsq, 2);
hsq = hsq(hsq<=1);
pi1=4*length(hsq)/n

%% Question 2
% weights are 1/100^2 where 100^2 is the # of draws (points)
x = transpose(0.01:0.01:1); % 100 by 1 vector running from 0.01 to 1
y = transpose(0.01:0.01:1); % 100 by 1 vector running from 0.01 to 1
grid =[kron(x, ones(100,1)), kron(ones(100,1), y)];
grid = grid.^2;
grid = sum(grid, 2);
grid = grid(grid<=1);
pi2 = 4*length(grid)/(100^2)

%% Question 3









