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
fprintf('[2] : METODO DELL''AUTOCOVARIANZA \n')
prompt = 'Digita [1] o [2]: \n';

tmpA = input(prompt);
a_met = met_choice(tmpA);

%%  Test Signals
fprintf('[1] : Unica sinusoide. \n')
fprintf('[2] : File audio percussivo, campionamento di un rullante. \n')
fprintf('[3] : Multiple sinusoidi. \n')
fprintf('[4] : Realizzazione di un processo AR di ordine 10. \n')
fprintf('[5] : Due sinusoidi corrotte da rumore bianco. \n')
fprintf('[6] : File audio musicale, musica classica. \n')

prompt = 'Scegli il tuo segnale di prova: \n';
tmpB = input(prompt);

switch tmpB
    case 1
    fprintf('Sinusoide @f= 500Hz, sequenza di 44100 sample, unico burst,m=4. \n')
    s_tmp = sin(2*pi*(500/44100).*[1:44100]);
    t = [1000 1001 1002 1003 ]; 

    case 2
    fprintf('File audio percussivo: fs=44100, unico burst di m=8. \n')
    [s_tmp, fs] = audioread('d.wav');
    t = linspace(121,128,8);
    s_tmp = s_tmp';
    sound(s_tmp, fs);
    
    case 3
        fprintf('Due sinusoidi, sequenza di 512 sample, con 3 bursts da m=4. \n')
        s_tmp = sin(0.23*pi*[1:512] +0.3*pi) + 0.6*sin(0.4*pi*[1:512] +0.3*pi);
        t = [100 101 102 103 221 222 223 224 341 342 343 344 ];
    case 4
        fprintf("Realizzazione, artificialmente generata di un processo AR")
        fprintf(" di ordine 10 con spettro generico\n")
        fprintf("Sequenze di 512 sample. Pattern dei burst con m=8. \n")
        
        mdl = regARIMA(10,0,0);
        mdl.Intercept =0 ;
        mdl.Variance = 1 ;
        mdl.AR = {1 -0.7185 0.5502 0.7970e-1 0.1586e-2 0.3802e-1 0.5296e-1 ...
            0.2792e-1 -0.11538e-1 -0.54262e-1 -0.5943e-2 };
        rng('default')
        s_tmp = simulate(mdl,512);
        s_tmp = s_tmp';
        t = [123 124 125 126 127 128 129 130 ...
             201 202 203 204 205 206 207 208 ];
         
     case 5
        fprintf('Due sinusoidi corrotte da rumore bianco. \n')
        fprintf('Sequenza di 512 sample, 3 bursts da m=4. \n')
        fprintf('SNR = 40dB \n')
        s_tmp = sin(0.23*pi*[1:512] +0.3*pi) + ...
            0.6*sin(0.4*pi*[1:512] +0.3*pi) + 0.026*(rand(1,512)-0.5);
        t = [100 101 102 103 221 222 223 224 341 342 343 344 ];
    case 6
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
    otherwise
        error('Nessun segnale di test corrisponde al codice inserito')
    
end



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
p = 3*m+2; %p: order of the AR process
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
%lista(i,:)=a;
%% Estimation of x
%[sig] = x_estimator(a_i, t_pos, sig)
sig = x_estimator(a,t,sig);

%lista(i,:)=sig;
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
if (tmpB == 2 || tmpB == 6) 
        fprintf('Questo è il segnale compromesso \n');
        sound(s_comp, fs);
        fprintf('Premi un tasto per proseguire \n')
        pause
        fprintf('Questo è il segnale ricostruito \n')
        sound(sig, fs);
end

fprintf('Done! \n')

