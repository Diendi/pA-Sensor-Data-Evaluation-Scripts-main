%characterization of various matlab lowpass implementations
%
%Matthias Musil
%FH Wels RING
%05.02.2024

clear
close all

%% setup
fs = 12.5e3; %sampling rate
f_stop = 100; %Hz, stop/cutoff frequency

%lowpass filter variables
filt_type = 'fir';
filt_steep = 0.95;

%filter design variables
Ap = 0.1; %dB, pass band attenuation
Ast = 60; %dB, stop band attenuation

%% test signal

%% lowpass filter as used so far in pA_VIchar.m

%load original filter objetc
load filter.mat

fvt1 = fvtool(filt_obj,'Fs',fs,'Color','white');
%--> f_pass = 10 Hz, but f_stop almost 400 Hz, unsuitable
%design my own filter here

%% FIR design
%doc: https://de.mathworks.com/help/dsp/ref/fdesign.html

%Fpass = f_stop/20;
%Fstop = f_stop;
Fpass = 10;
Fstop = 100;
spec = fdesign.lowpass(Fpass,Fstop,Ap,Ast,fs);
% Give more weight to passband
firFiltObj = design(spec,'equiripple','SystemObject',true);

fvt =  fvtool(firFiltObj,'Fs',fs,'Color','white');
%looks pretty good, try this one

save('fir_filt_10_100.mat','firFiltObj');
