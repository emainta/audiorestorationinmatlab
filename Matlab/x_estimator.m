function [sig] = x_estimator(a_i, t_pos, sig)
%   X_ESTIMATOR: Given an estimation of 'a', an estimation of 'x' is
%   returned
%   Q(a,x) is minimized as function of x (to obtain x). (Section III.B)


%% Prova

%{
clc
clear all
sig_or = sin(0.5*[1:1024]);
%sig = [ rand(1,20) 0 0 0 0 rand(1,20) ];
p_i =20;
a_i = arcov(sig_or,p_i)';
t_pos = linspace(100,120,21);
sig = sig_or;
sig(t_pos) = 0;
x_i= sig(t_pos(1):t_pos(end))';

figure();
subplot(3,1,1)
stem(sig([80:130]))
hold on
stem(sig_or(80:130)+0.001)
hold on
%}

%%  Var
p_i = length(a_i)-1; %length(a) = p + 1
m_i = length(t_pos);  %length(t_pos) = m ?
x_i = sig(t_pos(1):t_pos(m_i))';

% b(l) = b'(l+p+1)
b = zeros(1, 2*p_i+1); 

B = zeros(m_i);
z = zeros(m_i,1);

%a traslato di (p+1) e circondato di 0
a_tmp = [zeros(p_i,1); a_i; zeros(2*p_i, 1)]; 

%%  definition of B


for l = (-p_i):p_i
    for k = 1:(p_i+1)
        b(l+p_i+1) = b(l+p_i+1) + a_tmp(k+p_i)*a_tmp(k+p_i+l);
    end
end

%   B is an m x m matrix

for i = 1:m_i
    for j=1:m_i
        B(i,j) = b(t_pos(i) - t_pos(j) + p_i + 1);
    end
end

chol(B);

%%  definition of z

for i=1:m_i
    for k=1:(2*p_i+1)
        z(i)= z(i) + b(k)*sig(t_pos(i) - (k - p_i - 1));
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

%b(p_i+1)
x_i = B\(-z);
sig(t_pos) = x_i;




%% Prova - Un po' di plot
%{
sig_fin = sig;
sig_fin(t_pos) = x_i;

hold on
subplot(3,1,2)
stem(sig_or(80:130))
hold on
subplot(3,1,3)
stem(sig_fin(80:130))

err = max(sig_or(80:130) - sig_fin(80:130) );

fprintf("\n Errore massimo : " + err + "\n")


%}
end

