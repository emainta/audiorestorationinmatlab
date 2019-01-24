%function [x_i] = x_estimator(a_i, t_pos, sig)
%X_ESTIMATOR Given an estimation of 'a', an estimation of 'x' is given
%   Q(a,x) is minimized as function of x to obtain x. (Section III.B)


%% Prova
clc
clear all
p_i =14;
a_i= [1; rand(p_i,1)];
sig = [ rand(1,20) 0 0 0 0 rand(1,20) ];
t_pos = [ 21 22 23 24];
x_i= sig(t_pos(1):t_pos(4))';


%%  Var
p_i = length(a_i)-1; %length(a) = p + 1  
b = zeros(1, 2*p_i+1);
m_i = length(t_pos);  %length(t_pos) = m ?
B = zeros(m_i);
z = zeros(m_i,1);

a_i = [zeros(p_i,1); a_i; zeros(2*p_i, 1)];

%%  definition of B

%Mi serve realmente calcolarmi tutti i b(l) ?

for l = (-p_i):p_i
    for k = (p_i+1):(2*p_i+1)
        b(l+p_i+1) = b(l+p_i+1) + a_i(k)*a_i(k+l);
    end
end

%   B is an m x m matrix

for i = 1:m_i
    for j=1:m_i
        B(i,j) = b(t_pos(i)) - b(t_pos(j));
    end
end

chol(B);

%%  definition of z

for i=1:m_i
    for k=1:(2*p_i+1)
%       z(i)= z(i) + b(k)*sig(t_pos(i) - k+p_i+1);
    end
end
    

%%  W/ Cholesky decomposition  
%   (II.14) : B(a_).*x_ = -z(a_)
%
%   B(a) is positive definite => We can use Cholesky decomposition of B
%   for solving x in O(m^3) operations.

%L = chol(B,'lower');

%%  Simpler alternative ( To improve )
%   Per adesso ho deciso usare l'operatore mldivide (\)
%   x = C\E solves the system of linear equations C*x = E

x_i = B\(-z);
%end

