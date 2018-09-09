% Motoaki Takahashi
% HW1 for Econ 512 Empirical Methods
clear

diary hw1.out


%% Question 1
disp('Question 1')

X=[1, 1.5, 3, 4, 5, 7, 9, 10];
Y1=-2+0.5*X
Y2=-2+0.5*(X.^2)
Y=[Y1; Y2]
plot(X, Y)

%% Question 2
disp('Question 2')

X=linspace(-10, 20, 200)

%% Question 3
disp('Question 3')

A=[2, 4, 6;
   1, 7, 5;
   3, 12, 4];

b=[-2; 3; 10];

C=A.'*b
D=inv(A.'*A)*b
E=sum(b.'*A)
F=A([1,3],[1,2])

% solve the system of linear equations for the vector x
x=inv(A)*b

% for use on older versions of matlab it is faster to write A\b, it does
% the same thing, but it's much faster

%% Question 4
disp('Question 4')

B=kron(eye(5,5), A)

%% Question 5
disp('Question 5')

% 5X3 matrix whose elements are random draws from N(10, 5^2)
A=normrnd(10, 5, [5,3])
% B has an element 1 if the correspondent in A is bigger than or equal to 10
B=A>=10

%% Question 6
disp('Question 6')

% read the csv file
data='datahw1.csv';

% csvread replaces NaNs with zeros and that makes your estimates biased
data=csvread(data);
% X is a matrix of explanatory variables
X=data(1:size(data,1),[3,4,6]);
X=[ones(size(data,1),1), X];

% y is a vector of a explained variable
y=data(1:size(data,1), 5);

% OLS estimate
disp('OLS estimate')
b=inv(X.'*X)*(X.'*y)

% residual
e= y -X*b;
% the estimate for the variance of beta-hat (here b)
Var_hat=((e'*e)/(size(data,1)-4))*inv(X'*X);
disp('Standard Errors')
se=diag(Var_hat.^(1/2))

diary off