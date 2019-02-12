function [a_i] = a_estimator(sig, p_i, met)
%A_ESTIMATOR : Dato un segnale 'sig', questa funzione ritorna una stima di
%   'a' , vettore dei coefficienti del modello AR
%
%   IMPORTANTE: Passare solo a_i = a(2) ... (p+1)
%   a(1) è sempre 1. 

%   Alla chiamata di questa funzione è possibile usare uno tra i due metodi
%   descritti nel dettaglio nella relazione:
%   - se (met == 'acov') verrà usato il METODO DELL'AUTO-COVARIANZA
%   - se (met == 'acor') verrà usato il METODO DELL'AUTO-CORRELAZIONE

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

%% Check di met - METODO DELL'AUTOCORRELAZIONE
%   R.*a = -r

if (met == 'acor') 
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

%% Check di met - METODO DELL'AUTOCOVARIANZA
elseif (met == 'acov') 
    a_est = arcov(sig,p_i);
    a_i = a_est(2:end)';

else
    fprintf("Errore nella richiesta del metodo di stima di 'a_i'")
end
%{
a_i = [1;a_i];
a_est = arcov(sig,p_i);
ahi= zeros(p_i+1,2);
ahi(:,1) = a_i;
ahi(:,2) = a_est';

%}
end

