% Motoaki Takahashi
% HW7 Econ 512 Empirical Method

clear

global L rho eta kappa l v delta beta lambda c CRIT Omega trans;
L = 30;
rho = 0.85;
eta = log(rho)/log(2);
kappa = 10;
l = 15; % at which the learning curve flattens out
v = 10;
delta = 0.03;
beta = 1/1.05;
lambda = 0.1; % for dampning
% c is the cost
c = zeros(L, 1);
c(1:l) = kappa * [1:l] .^ (eta);
c(l+1:L) = kappa * (l * ones(L-l, 1)) .^ (eta);
CRIT = 1e-5;
Omega = kron((1:L)', ones(1, L));
trans = zeros(L, L, 2); % this will contain transition matrices for q=0, 1
trans(:, :, 1) = diag([1, (1-delta) .^ (2:L)]) + diag(ones(1, L-1) - (1-delta) .^ (2:L), -1);
trans(:, :, 2) = diag([ones(1, L-1) - (1-delta) .^ (1:L-1), 1]) + diag((1-delta) .^ (1:L-1) ,1);







%% Question 1

p_initial = ones(L, L);
V_initial = ones(L, L);

check = 1;
iter = 0;

p = p_initial;
V = kron((1:30)', ones(1, 30))';


while check > CRIT && iter < 10000
    W = getW(V, trans);
    np = getp(p, W);
    nV = getV(np, p, W);
    check = max( max(max(abs((nV-V)./(1+nV)))),  max(max(abs((np-p)./(1+np)))))
    V = lambda * nV + (1-lambda)*V;
    p = lambda * np + (1-lambda)*p;

    iter = iter+1;
end

figure(1);
mesh(V);
title('Value Function');

figure(2);
mesh(p);
title('Policy');

%% Question 2

% construct a 900 by 900 transition matrix
% elements are like (omega1, omega2) = (1,1), (1,2), ..., (1,30), ...,
% (30,30)

big_trans = zeros(L*L, L*L);
Dmat = D(p, p');
for omega1 = 2:L-1
    for omega2 = 2:L-1
        D1 = Dmat(omega1, omega2);
        D2 = Dmat(omega2, omega1);
        D0 = 1 - D1 - D2;
        
        big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2-1) = D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
        big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2) = D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2);
        big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2+1) = D2 * (1-(1-delta)^omega1) * ((1-delta)^omega2);
        
        big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2-1) = D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
        big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2) = D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
        big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2+1) = D2 * ((1-delta)^omega1) * ((1-delta)^omega2);
        
        big_trans(L*(omega1-1)+omega2, L*(omega1)+omega2-1) = D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
        big_trans(L*(omega1-1)+omega2, L*(omega1)+omega2) = D1 * ((1-delta)^omega1) * ((1-delta)^omega2);
    end
end

% the cases where omega1 or omega2 is 1 or L
% omega1 = 1 and omega2 = 1
Dmat = D(p, p');
D1 = Dmat(1, 1);
D2 = Dmat(1, 1);
D0 = 1 - D1 - D2;
big_trans(1, 2) = D2 * (1-delta) ^ 1;
big_trans(1, 30+1) = D1 * (1-delta) ^ 1;
big_trans(1, 1) = 1 - D2 * (1-delta) ^ 1 - D1 * (1-delta) ^ 1;

% omega1 = 1 and 2=<omega2=<L-1
omega1 = 1;
Dmat = D(p, p');
for omega2 = 2:L-1
    D1 = Dmat(omega1, omega2);
    D2 = Dmat(omega2, omega1);
    D0 = 1 - D1 - D2;
    
    
    big_trans(omega2, omega2-1) = D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
        D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
        D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
    big_trans(omega2, omega2) = D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2)+...
            D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
    big_trans(omega2, omega2+1) = D2 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
        D2 * ((1-delta)^omega1) * ((1-delta)^omega2);
    
    big_trans(omega2, L+omega2-1) = D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
    big_trans(omega2, L+omega2) = D1 * ((1-delta)^omega1) * ((1-delta)^omega2);
end

% 2=<omega1=<L-1, omega2=1
omega2 = 1;
for omega1 = 2:L-1
    D1 = Dmat(omega1, omega2);
    D2 = Dmat(omega2, omega1);
    D0 = 1 - D1 - D2;
    
    
    big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2) = D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
        D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
        D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2);
    big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2+1) = D2 * (1-(1-delta)^omega1) * ((1-delta)^omega2);
    
    big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2) = D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
    big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2+1) = D2 * ((1-delta)^omega1) * ((1-delta)^omega2);
    
    big_trans(L*(omega1-1)+omega2, L*(omega1)+omega2) = D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
        D1 * ((1-delta)^omega1) * ((1-delta)^omega2);
end

% omega1=L, 2=<omega2=<L-1
omega1 = L;
Dmat = D(p, p');
for omega2 = 2:L-1
    D1 = Dmat(omega1, omega2);
    D2 = Dmat(omega2, omega1);
    D0 = 1 - D1 - D2;
    
    big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2-1) = D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
    big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2) = D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
           D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2);
    big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2+1) = D2 * (1-(1-delta)^omega1) * ((1-delta)^omega2);
    
    big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2-1) = D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
    big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2) = D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * ((1-delta)^omega1) * ((1-delta)^omega2);
    big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2+1) = D2 * ((1-delta)^omega1) * ((1-delta)^omega2);
end

% omega2=L, 2=<omega1=<L-1
omega2 = L;
Dmat = D(p, p');
for omega1 = 2:L-1
    D1 = Dmat(omega1, omega2);
    D2 = Dmat(omega2, omega1);
    D0 = 1 - D1 - D2;
    
       big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2-1) = D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
       big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2) = D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
           D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2)+...
           D2 * (1-(1-delta)^omega1) * ((1-delta)^omega2);
       
       big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2-1) = D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
       big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2) = D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * ((1-delta)^omega2);
        
       big_trans(L*(omega1-1)+omega2, L*(omega1)+omega2-1) = D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
       big_trans(L*(omega1-1)+omega2, L*(omega1)+omega2) = D1 * ((1-delta)^omega1) * ((1-delta)^omega2);

end

% omega1=1, omega2=L
omega1=1;
omega2=L;
D1 = Dmat(omega1, omega2);
D2 = Dmat(omega2, omega1);
D0 = 1 - D1 - D2;

big_trans(omega2, omega2-1)= D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
    D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
    D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
big_trans(omega2, omega2) = D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2)+...
            D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * ((1-delta)^omega2);
big_trans(omega2, L*(omega1-1+1)+omega2-1) = D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
big_trans(omega2, L*(omega1-1+1)+omega2) = D1 * ((1-delta)^omega1) * ((1-delta)^omega2);

% omega1=L, omega2=1
omega1=L;
omega2=1;
D1 = Dmat(omega1, omega2);
D2 = Dmat(omega2, omega1);
D0 = 1 - D1 - D2;

big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2) = D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
    D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
    D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2);
big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2+1) = D2 * (1-(1-delta)^omega1) * ((1-delta)^omega2);

big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2) = D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * ((1-delta)^omega1) * ((1-delta)^omega2);
big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2+1) = D2 * ((1-delta)^omega1) * ((1-delta)^omega2);

% omega1=L, omega2=L
omega1=L;
omega2=L;
D1 = Dmat(omega1, omega2);
D2 = Dmat(omega2, omega1);
D0 = 1 - D1 - D2;

big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2-1) = D0 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2);
big_trans(L*(omega1-1)+omega2, L*(omega1-1-1)+omega2) = D0 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2*(1-(1-delta)^omega1)*(1-(1-delta)^omega2)+...
            D2 * (1-(1-delta)^omega1) * ((1-delta)^omega2);
big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2-1) = D0 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D1 * ((1-delta)^omega1) * (1-(1-delta)^omega2);
big_trans(L*(omega1-1)+omega2, L*(omega1-1)+omega2) = D0 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * (1-(1-delta)^omega1) * ((1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * (1-(1-delta)^omega2)+...
            D2 * ((1-delta)^omega1) * ((1-delta)^omega2)+...
            D1 * ((1-delta)^omega1) * ((1-delta)^omega2);

 max(sum(big_trans, 2)-ones(L*L, 1))

% draw the distribution of states after 10, 20, 30 periods, starting with
% (1,1)

% initial state
state_int = zeros(1,L*L);
state_int(1,1) = 1;

% 10 period
Dstrbn10 = state_int * (big_trans)^10;
Dstrbn10 = reshape(Dstrbn10, [30, 30]);



% 20 period
Dstrbn20 = state_int * (big_trans)^20;
Dstrbn20 = reshape(Dstrbn20, [30, 30]);


% 30 period
Dstrbn30 = state_int * (big_trans)^30;
Dstrbn30 = reshape(Dstrbn30, [30, 30]);

figure(3);
mesh(Dstrbn10);
title('Distribution of states after 10 periods');

figure(4);
mesh(Dstrbn20);
title('Distribution of states after 20 periods');

figure(5);
mesh(Dstrbn30);
title('Distribution of states after 30 periods');

%% Question 3

% Iterate the transition until the dist of states converges

state = state_int;
diff = 1;
while diff > CRIT
    state_new = state * big_trans;
    diff = max(abs(state - state_new))
    state = state_new;
end

stn = reshape(state, [30, 30])

figure(6);
mesh(stn);
title('Stationary Distribution');