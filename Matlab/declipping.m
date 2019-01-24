% ---------------------------------------------------
%   Adaptive Interpolation of Discrete-Time Signals 
%   That Can Be Modeled as Autoregressive Processes
% ---------------------------------------------------
clc
clear all
close all

%%  Prova
s_tmp = sin(0.23.*pi.*[1:512]);

%%  Given
t = [100 101 102 103]; %t: vect of the unknown samples indexes position

%% Variables
N = length(s_tmp); %number of samples of the available data
m = length(t); %m: n' of unknown samples
x = zeros(1,m); %x: vect of the unknown samples
p = 3*m+2; %p: order of the AR process
a = [1 zeros(1,length(p))].' ; %a: col vect of the prediction coeff., a(1)=1 
%
fprintf('%i ', s_tmp(t))
fprintf('\n')
%% Sub-optimal approach

%% Estimation of a
a = a_estimator(x,t,s_tmp,p);

%% Estimation of x
x = x_estimator(a,t,s_tmp);

fprintf('%i ', x)
fprintf('\n')
