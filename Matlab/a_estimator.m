function [a_i] = a_estimator(sig, p_i)
%A_ESTIMATOR Given an estimation of 'x', an estimation of 'a' is given
%   Q(a,x) is minimized as function of a to obtain a. (Section III.A)
%   We'll use the AUTOCORRELATION METHOD
%   R.*a = -r

%   Also Try estimate : Estimate ARIMA or ARIMAX model parameters
%   Also try ar: Estimate parameters of AR model for scalar time series
%
%   a_i = a(2) ... (p+1);
%%  PROVA
%{
clc
clear all
sig = sin(0.23.*pi.*[1:100]);
t_pos = [ 31 32 33 ];
p_i = 8;

figure();
subplot(2,1,1)
stem(sig)
hold on
%}

%%  Var
N_i= length(sig);
r_tmp = zeros(2*p_i+1,1);
r = zeros(p_i,1);
R = zeros(p_i);

%%  definition of r
%   r(j): a (biased) estimate for the j-th lag of s

for j= (-p_i):p_i
    for k=1:(N_i - abs(j))
        r_tmp(p_i+j+1)= r_tmp(p_i+j+1) + sig(k)*sig(k+abs(j));
    end
end
r_tmp = r_tmp./N_i;
r = r_tmp((p_i+2):end); % r=[r(1),...,r(p)]

%%  definition of R (Autocorrelation matrix)

for j=1:p_i
    for i = 1:p_i
        R(i,j) = r_tmp(p_i+1+(i-j));
    end
end

%Nota: R è simmetrica

%%  Solving the system obtained 
%   R.*a = -r   
%       (Section III.A)
a_i = R\(-r);

%{
a_i = [1;a_i];
a_est = arcov(sig,p_i);
ahi= zeros(p_i+1,2);
ahi(:,1) = a_i;
ahi(:,2) = a_est';

%}
end

