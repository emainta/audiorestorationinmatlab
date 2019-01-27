% ---------------------------------------------------
% | Adaptive Interpolation of Discrete-Time Signals | 
% | That Can Be Modeled as Autoregressive Processes |
% ---------------------------------------------------
clc
clear all
close all

%% Richieste a utente
prompt = 'N iterazioni? ';
n_it = input(prompt);

%   Scelta del metodo di stima dei parametri del modello AR
fprintf('Scelta del metodo di stima dei parametri del modello AR: \n')
fprintf('[1] : METODO DELL''AUTOCORRELAZIONE \n')
fprintf('[2] : METODO DELL''AUTOCOCOVARIANZA \n')
prompt = 'Digita [1] o [2]: \n';

tmpA = input(prompt);
a_met = met_choice(tmpA);

%%  Prova
%s_tmp = sin(0.23.*pi.*[1:512]);

[s_tmp, fs] = audioread('d.wav');

sig = s_tmp;

figure();
title("Questo è solo un test")
subplot(3,1,1)
ylabel("Segnale Originale")
plot(sig(90:150))
hold on

%%  Given
%t: vect of the unknown samples indexes position
t = linspace(100,110,11);
%t = [100 101 102 103 ]; 

sig(t) = 0 ;
hold on
subplot(3,1,2)
ylabel("Segnale Compromesso")
plot(sig(90:150))

%% Variables
N = length(sig); %number of samples of the available data
m = length(t'); %m: n' of unknown samples
x = zeros(1,m); %x: vect of the unknown samples
p = 3*m+2; %p: order of the AR process
a = [1 zeros(1,p)].' ; % a: col vect of the prediction coeff., a(1)=1 
                               % Remember : length(a) = p+1
%% Prova

%% Check the input

%% Sub-optimal approach
tic

for i=1:n_it
%% Estimation of a
%[a_i] = a_estimator(sig, p_i, met)
a(2:end) = a_estimator(sig, p, a_met);

%% Estimation of x
%[sig] = x_estimator(a_i, t_pos, sig)
sig = x_estimator(a,t,sig);

sig_ = x_estimator(a_,t,sig);
end

subplot(3,1,3)
ylabel("Segnale Ricostruito")
plot(sig(90:150))

toc
fprintf('Done!')

