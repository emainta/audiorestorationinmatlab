% ---------------------------------------------------
% | Adaptive Interpolation of Discrete-Time Signals | 
% | That Can Be Modeled as Autoregressive Processes |
% ---------------------------------------------------
clc
clear all
close all

%% Richieste a utente
fprintf('[6] : File audio musicale, musica classica. \n')
prompt = 'N iterazioni? ';
n_it = input(prompt);

%   Scelta del metodo di stima dei parametri del modello AR
fprintf('Scelta del metodo di stima dei parametri del modello AR: \n')
fprintf('[1] : METODO DELL''AUTOCORRELAZIONE \n')
fprintf('[2] : METODO DELL''AUTOCOVARIANZA \n')
prompt = 'Digita [1] o [2]: \n';

tmpA = input(prompt);
a_met = met_choice(tmpA);

%%  Test Signals

    fprintf('File audio musicale:  \n')
    fprintf('Ludwig van Beethoven - Ode to Joy - Symphony No. 9 in D minor Choral Op. 125 \n')
    fprintf('Durata: 17 secondi\n')
    fprintf('Fs = 44100 Hz \n')
    fprintf('Bursts di m=16. \n')
    fprintf('Rate dei bursts: 1/10 s \n')
    
    [s_tmp, fs] = audioread('beethoven.wav'); 
    %Definizione di t particolare per questo caso
    i_6 = 5*fs;
    k_6 = 1;
    t(1) =0;
    fprintf('Richiede tempo... \n');
        while ( i_6 < 10*fs )
           t(k_6:k_6+16) = [i_6:(i_6+16)];
          k_6 = k_6 + 16;
           i_6 = i_6 + 0.1*fs;
         end
    
    
    clear i_6 k_6
    s_tmp = s_tmp';
    fprintf('NOW PLAYING: File originale\n')
    sound(s_tmp, fs);    
           
%% SCELTA P

prompt = 'Che ''p'' usare? \n';
fprintf('[1] : AUTO ( p=3m+2 ) \n')
fprintf('[2] : p = 50 \n')
fprintf('[3] : p = min(3m+2,50) \n')
tmpP = input(prompt);
%%

sig = s_tmp;

figure('Name','Experiment','NumberTitle','off');
subplot(3,1,1)
stem([(t(1)-20):(t(end)+20)],sig((t(1)-20):(t(end)+20)), ...
    'LineStyle',':','Marker','.','MarkerSize',7)
title("Segnale Originale")
xlabel("n-esimo Sample",'FontSize',6)
hold on

%%  Given
%t: vect of the unknown samples indexes position

sig(t) = 0 ;
s_comp = sig;
hold on
subplot(3,1,2)
stem([(t(1)-20):(t(end)+20)], sig((t(1)-20):(t(end)+20)), ...
        'LineStyle',':','Marker','.','MarkerSize',7) 
hold on
stem(t, sig(t), ...
        'LineStyle',':','Marker','.','MarkerSize',7) 
%stem([21:(21+length(t)-1)], sig((t(1)):(t(end))))
title("Segnale Compromesso")
xlabel("n-esimo Sample",'FontSize',6)

%% Variables
N = length(sig); %number of samples of the available data
m = length(t'); %m: number of unknown samples
x = zeros(1,m); %x: vect of the unknown samples

switch tmpP
    case 1
        p = 3*m+2; %p: order of the AR process
    case 2 
        p = 50;
    otherwise
        p = min(3*m+2 , 50);
end

a = [1 zeros(1,p)].' ;  % a: col vect of the prediction coeff., a(1)=1 
                        % Remember : length(a) = p+1
%% Prova
%lista= zeros(n_it,length(a));

%% Check the input

if ( t(1) >= (p+1) && ...
        t(end) <= (N-p) )
    fprintf("Input idoneo. \n");
else
        error("Input non idoneo.");
    end
    
%% Sub-optimal approach
tic

for i=1:n_it
%% Estimation of a
%[a_i] = a_estimator(sig, p_i, met)
a(2:end) = a_estimator(sig, p, a_met);

%% Estimation of x
%[sig] = x_estimator(a_i, t_pos, sig)
sig = x_estimator(a,t,sig);

end

subplot(3,1,3)
stem([(t(1)-20):(t(end)+20)],sig((t(1)-20):(t(end)+20)), ...
        'LineStyle',':','Marker','.','MarkerSize',8)
hold on
stem([(t(1)-20):(t(end)+20)],s_tmp((t(1)-20):(t(end)+20)), ...
        'LineStyle',':','Marker','o','MarkerSize',4)
legend('Ricostruito','Originale')
title("Segnale Ricostruito")
xlabel("n-esimo Sample",'FontSize',6)

toc

%% Calcolo ERORRE QUADRATICO RELATIVO MEDIO di interpolazione
% [e] = eqrm(s_or, s_ric, t_poz)
e_i = eqrm(s_tmp, sig, t);

fprintf("EQRM: " + e_i + " \n")
%% Calcolo di var(e)^2
%{
sigma_e_sq = sigma_calc(N,p,m,a,sig);
[a_, sigma_e_sq2] = arcov(sig,p);
fprintf("Sigma: " + sigma_e_sq + " \n")
fprintf("Sigma: " + sigma_e_sq2 + " \n")
%}

prompt = 'Calcolare Sigma(e)^2?    [ y/n ] \n';
confirm = input(prompt,'s');
if (confirm == 'y')
prompt = 'Calcolo di sigma(e)^2 \n';
fprintf("[1] : Algoritmo suggerito (Eq. II.15) \n");
fprintf("[2] : Funzione arcov \n");
e_met = input(prompt);

switch e_met
    case 1
        sigma_e_sq = sigma_calc(N,p,m,a,sig);
    case 2
        [a_, sigma_e_sq] = arcov(sig,p);
        clear a_
    otherwise
        fprintf("Comando non riconosciuto \n")
end

fprintf("Sigma: " + sigma_e_sq + " \n")

end

%%  Scostamento massimo dal segnale originale
delta = max(abs(s_tmp-sig));
delta_perc = delta/max(sig);
delta_perc = 100*round(delta_perc,3);
fprintf("Scostamento massimo: " + delta_perc + "%% \n")

%% Riproduzione per i Test Audio

        fprintf('Questo è il segnale compromesso \n');
        sound(s_comp, fs);
        fprintf('Premi un tasto per proseguire \n')
        pause
        fprintf('Questo è il segnale ricostruito \n')
        sound(sig, fs);


fprintf('Done! \n')

