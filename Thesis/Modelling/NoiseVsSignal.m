%Thesis, Modelling the noise unburying possible with a know relationship
%between electrical signal and vestibular signal

clear all;
close all;
clc;

%% Constants

theta = pi/2; %heading of the shark in radians
B = 40*10^-6; %strength of the magnetic field in Tesla
veloc = 1; %forward velocity of the shark in m/s
omega = 1; %2pi*vestibular frequency in Hz
headAmp = 0.1; %max amplitude of the head movement in m
lambdaV = 0.8; %proportionality constant related to sensitivity for V
lambdaM = 0.99; %proportionality constant related to sensitivity for M
restFire = 34; %resting firing rate of the neurons in Hz
sigmaV = 1; %std dev of noise in Voltage
sigmaM = 1; %std dev of noise in Movement
modNoiseV = 10^-7; %the order of magnitude of noiseV
modNoiseM = 10^-3; %the order of magnitude of noiseM
dt = 0.1; %sampling time for the neuron
length = 1000; %length of sampling period

%Vector quantities
time = (0:length-1)*dt; %time vector in s
accel = -(headAmp*(omega^2))*sin(omega*time); %head acceleration
psi = -(1/(omega^2))*accel + theta; %psi gives the angle between v and B
vElec = lambdaV*veloc*B*sin(psi); %the potential the shark measures minus noise

noiseV = modNoiseV*gaussNoise(sigmaV, 0, time); %the noise in vSignal
noiseM = modNoiseM*gaussNoise(sigmaM, 0, time); %the noise in mSignal


%Signals reaching the shark
Vsignal = noiseV + vElec; %what the shark measures
Msignal = noiseM + lambdaM*accel;


%% Plotting initial functions
% plot(time, Msignal); %plotting the vesitbular signal + noise
% xlabel('time (s)'); 
% ylabel('amplitude');
% title('Vestibular Signal vs Time');
figure;
plot(time, Vsignal);
xlabel('time (s)'); 
ylabel('Measured Voltage');
title('Voltage Signal vs Time');
figure;
plot(time, noiseV);
figure;
plot(time, vElec);


%% Comparing Vsignal against the vestibular signal
%Going to try and take the FFT of Msignal and Vsignal and then get the main
%component of the frequency of F(Msignal) and then just take info from
%Vsignal that has frequency in integer multiples of the main component of
%F(Msignal).

NFFT = 2^nextpow2(length); %next power of 2 from length of time
F_Msignal = fft(Msignal, NFFT)/length; %FT of Msignal
F_Vsignal = fft(Vsignal,NFFT)/length; %FT of Vsignal
freq = 1/(2*dt)*linspace(0,1,NFFT/2+1); %frequency vector

%Plotting the FT's of Msignal and Vsignal
% figure;
% plot(freq,2*abs(F_Msignal(1:NFFT/2+1)), 'r'); 
% figure;
% plot(freq,2*abs(F_Vsignal(1:NFFT/2+1)), 'k');

correlated = ifft((modNoiseV/modNoiseM)*F_Msignal.*F_Vsignal);






