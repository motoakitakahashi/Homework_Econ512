% HW2 for Econ 512 Empirical Methods
% Motoaki Takahashi

clear

%% Question 2
disp('Question 2')

% Define the quality vector (2 by 1)
q = [2; 2];

% Refer to the function demand, demand.m, which returns the 2 by 1 demand vector for
% a 2 by 1 price vector and a 2 by 1 quality vector

% eqn.m defines the "left-hand side" of the equation of interest.
% That is, we want to get the value p that solves eqn(p, q)=0 for given
% quality q.

% We are interested in eqn(p,q)=0 for q = [2; 2]. Name this function
% eqn_22.
eqn_22 = @(p) eqn(p,q);

% I follow and quote the code for Broyden method in the lecture note.
% myJac.m is taken from the lecture notes.

% the initial values for the equlibrium prices
p = [1; 1];

fVal = eqn_22(p)

%iJac = inv(myJac(eqn_22, p))
iJac = eye(size(p,1));
% the Broyden iterations

maxit = 100;
tol = 1e-6;
for iter = 1:maxit
    fnorm = norm(fVal);
    fprintf('iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n', iter, p(1), p(2), norm(fVal));
    if norm(fVal) < tol
        break
    end
    d = - (iJac * fVal);
    p = p+d;
    fOld = fVal;
    fVal = eqn_22(p);
    u = iJac*(fVal - fOld);
    iJac = iJac + ( (d - u) * (d'*iJac) )/ (d'*u);
end

%% Question 2

p = [1; 1];
pOld = [0;0];

maxit = 100;
tol = 1e-6;

for iter = 1:maxit;
    if eqn_22(p)<tol
        fprintf('Converged: iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n', iter, p(1), p(2), norm(eqn_22(p)));
        break
    end
    
    % Secant iteration
    
    eqn_A = @(pa) eqn_f([pa; p(2,1)], q, 1);
    
    % Secant for firm A
    pa = p(1, 1);
    paOld = pOld(1,1);
    fOld = eqn_A(paOld);
    
    tolA = 1e-8;
    for iter =1:maxit
    fVal = eqn_A(pa);
    fprintf('iter %d: pa = %.8f, f(x) = %.8f\n', iter, pa, fVal);
    if abs(fVal) < tolA
        break
    else
        paNew = pa - ( (pa - paOld) / (fVal - fOld) )* fVal;
        paOld = pa;
        pa = paNew;
        fOld = fVal;
    end
    end
    
    % Secant for firm B
    eqn_B = @(pb) eqn_f([pa; pb], q, 2);
    pb = p(2, 1);
    pbOld =  pOld(2,1);
    fOld = eqn_B(pbOld);
    
    for iter =1:maxit
    fVal = eqn_B(pb);
    fprintf('iter %d: pb = %.8f, f(x) = %.8f\n', iter, pb, fVal);
    if abs(fVal) < tolA
        break
    else
        pbNew = pb - ( (pb - pbOld) / (fVal - fOld) )* fVal;
        pbOld = pb;
        pb = pbNew;
        fOld = fVal;
    end
    p=[pa; pb]
    end
    
    
end

%% Question 4

%initial guess for p
p = [1; 1];

maxit = 100;
tol = 1e-6;
for iter = 1:maxit
    pnew = ones(2,1)./(ones(2,1)-demand(p, q));
    fprintf('iter %d: p(1) = %f, p(2) = %f, norm(p-pnew) = %.8f\n', iter, p(1), p(2), norm(p-pnew));
    if norm(p-pnew) < tol
        fprintf('Converged: iter %d: p(1) = %f, p(2) = %f, norm(p-pnew) = %.8f\n', iter, p(1), p(2), norm(p-pnew));
        break
    end
    p = pnew;
end
p


%% Question 5



