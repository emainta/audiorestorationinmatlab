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

%%  Test Signals
fprintf('[1] : Sinusoide (f= 500Hz), sequenza di 44100 sample, m=4. \n')
fprintf('[2] : File audio percussivo, fs=44100, m=8. \n')


prompt = 'Scegli il tuo segnale di prova: \n';
tmpB = input(prompt);

switch tmpB
    case 1
    s_tmp = sin(0.23*[1:44100]);
    t = [1000 1001 1002 1003 ]; 

    case 2
    [s_tmp, fs] = audioread('d.wav');
    t = linspace(120,128,9);
    s_tmp = s_tmp';

    otherwise
        fprintf('Nessun segnale di test corrisponde al codice inserito\n')
    
end



%%

sig = s_tmp;

figure();
title("Questo è solo un test")
subplot(3,1,1)
stem(sig((t(1)-20):(t(end)+20)))
title("Segnale Originale")
hold on

%%  Given
%t: vect of the unknown samples indexes position
%t = linspace(100,110,11);
%t = [100 101 102 103 ]; 

sig(t) = 0 ;
hold on
subplot(3,1,2)
stem(sig((t(1)-20):(t(end)+20)))
hold on

stem([20:(20+length(t)-1)], sig((t(1)):(t(end))))
title("Segnale Compromesso")

%% Variables
N = length(sig); %number of samples of the available data
m = length(t'); %m: n' of unknown samples
x = zeros(1,m); %x: vect of the unknown samples
p = 3*m+2; %p: order of the AR process
a = [1 zeros(1,p)].' ; % a: col vect of the prediction coeff., a(1)=1 
                               % Remember : length(a) = p+1
%% Prova
lista= zeros(n_it,length(sig));
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

lista(i,:)=sig;
end

subplot(3,1,3)
stem(sig((t(1)-20):(t(end)+20)))
hold on
stem(s_tmp((t(1)-20):(t(end)+20)))
legend('Ricostruito','Originale')
title("Segnale Ricostruito")

toc
fprintf('Done! \n')

