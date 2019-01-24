function [a_i] = a_estimator(x_i,t_pos, sig, p_i)
%A_ESTIMATOR Given an estimation of 'x', an estimation of 'a' is given
%   Q(a,x) is minimized as function of a to obtain a. (Section III.A)
%   We'll use the AUTOCORRELATION METHOD
%   R.*a = -r

%   Also Try estimate : Estimate ARIMA or ARIMAX model parameters
%   Also try ar: Estimate parameters of AR model for scalar time series
%%  PROVA
%clc
%clear all
%sig = [ rand(1,10) 0 0 rand(1,10)];
%x_i = rand(1,2);
%t_pos = [ 11 12 ];
%p_i = 8;

%%  Var
N_i= length(sig);
sig(t_pos) = x_i; %Estimated signal
r_tmp = zeros(2*p_i+1,1);
r = zeros(p_i,1);
R = zeros(p_i, p_i);

%%  definition of r
%   r(j): a (biased) estimate for the j-th lag of s

for j= (-p_i):p_i
    for k=1:(N_i - abs(j))
        r_tmp(p_i+j+1)= r_tmp(p_i+j+1) + sig(k)*sig(k+abs(j));
    end
end

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

end

