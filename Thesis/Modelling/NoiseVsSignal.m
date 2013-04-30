%Thesis, Modelling the noise unburying possible with a know relationship
%between electrical signal and vestibular signal

clear all;
close all;
clc;

%% Constants

theta = 0; %heading of the shark in radians
B = 40*10^-6; %strength of the magnetic field in Tesla
veloc = 1; %forward velocity of the shark in m/s
omega = 2*pi; %2pi*vestibular frequency in Hz
headAmp = 0.2; %max amplitude of the head movement in m
lambdaV = 0.8; %proportionality constant related to sensitivity for V
lambdaM = 0.99; %proportionality constant related to sensitivity for M
sigmaV = 0.001; %std dev of noise in Voltage
sigmaM = 1; %std dev of noise in Movement
modNoiseM = 10^-3; %the order of magnitude of noiseM
dt = 0.001; %sampling time for the neuron
length = 8192; %length of sampling period
numNeurons = 2000; %number of neurons
relaxation = 0; %relaxation time in s

%Vector quantities
time = (0:length-1)*dt; %time vector in s

accel = -(headAmp*(omega^2))*sin(omega*time); %head acceleration
psi = -(1/(omega^2))*accel + theta; %psi gives the angle between v and B
vElec = lambdaV*veloc*B*sin(psi); %the potential the shark measures minus noise

%had to change the noiseV scaling as the size of V changes significantly
%for different theta
modNoiseV = (max(vElec) - min(vElec)); %the order of magnitude of noiseV
noiseV = modNoiseV*gaussNoise(sigmaV, 0, time); %the noise in vSignal
noiseM = modNoiseM*gaussNoise(sigmaM, 0, time); %the noise in mSignal


%Signals reaching the shark
Vsignal = noiseV + vElec; %what the shark measures
Msignal = noiseM + lambdaM*accel; %what the shark 'feels'


%% Plotting initial functions
% plot(time, Msignal); %plotting the vesitbular signal + noise
% xlabel('time (s)'); 
% ylabel('amplitude');
% title('Vestibular Signal vs Time');
% figure;
% plot(time, Vsignal);
% xlabel('time (s)'); 
% ylabel('Measured Voltage');
% title('Voltage Signal vs Time');
%figure;
% plot(time, noiseV);
% xlabel('time (s)'); 
% ylabel('Noise Voltage');
% title('Noise Voltage vs Time');
% figure;
% plot(time, vElec);
% xlabel('time (s)'); 
% ylabel('Induced Voltage');
% title('Induced Voltage vs Time');

%% Converting Vsignal into an actionPotential and back
actionpotentials = zeros(numNeurons, length); %starting matrix for all neuron inputs

%the loop below creates an mxlength matrix with entries from all neurons
for m = 1:numNeurons
    actionpotentials(m, :) = actionPot(vElec, relaxation, dt); %each row is one neuron
end

actionSignal = invActPot(actionpotentials, dt); %what the sharks brain resolves
% figure;
% plot(time, actionSignal)
% xlabel('Time (s)');
% ylabel('Neuron Signal at brain');
% title('All inputs decoded');


%% Comparing actionSignalignal against the vestibular signal
% Going to try and take the FFT of Msignal and Vsignal and then get the main
% component of the frequency of F(Msignal) and then just take info from
% Vsignal that has frequency in integer multiples of the main component of
% F(Msignal).

NFFT = 2^nextpow2(length); %next power of 2 from length of time
F_Msignal = fft(Msignal, NFFT)/length; %FT of Msignal
F_actSignal = fft(actionSignal,NFFT)/length; %FT of Vsignal
freq = 1/(2*dt)*linspace(0,1,NFFT/2+1); %frequency vector
modF_Msignal = abs((1/max(2*abs(F_Msignal(1:NFFT/2+1)))))*(F_Msignal); 
                                        %F_Msignal with max = 1
peak = find(abs(modF_Msignal((1:NFFT/2+1))) == max(abs(modF_Msignal((1:NFFT/2+1)))));
                                        %finds the peak of modF_Msignal
peakvec = modF_Msignal(1:2*peak);                                        
peakvec = [zeros(1,peak) peakvec zeros(1,size(F_Msignal,2))];
filtervec = 2*abs(modF_Msignal(1:NFFT/2+1))+ 2*abs(peakvec(1:NFFT/2+1)); %what picks out the
                                                %frequencies at omega and
                                                %2*omega
P_filtervec = [filtervec, zeros(1, size(F_Msignal,2) - size(filtervec,2))];
                                    %filtervec padded with zeros
filtered_V = P_filtervec.*(F_actSignal); %filters Vsignal
% Plotting the FT's of Msignal and Vsignal
figure;
plot(freq,2*abs(F_Msignal(1:NFFT/2+1)), 'r');
xlabel('Frequency (f)'); 
ylabel('Amplitude of FT of Msignal');
title('FT of Msignal vs Frequency');
figure;
plot(freq,2*abs(F_actSignal(1:NFFT/2+1)), 'k');
xlabel('Frequency (f)'); 
ylabel('Amplitude of FT of actSignal');
title('FT of actSignal vs Frequency');
figure;
plot(freq,filtervec, 'k');
xlabel('Frequency (f)'); 
ylabel('Amplitude of filtervec');
title('Filtervec vs Frequency');
figure;
plot(freq, 2*abs(filtered_V(1:NFFT/2+1)));
xlabel('Frequency (f)'); 
ylabel('Amplitude of filtered FT of actSignal');
title('Filtered FT of Vsignal vs Frequency');

correlated = (length^2)*ifft((modNoiseV/modNoiseM)*(filtered_V));
%here I found a strange problem with the correlated signal, it doesn't
%always sync up with the start i.e. it doesn't start at zero as Vsignal
%does

ampCor = 0.5*(max(correlated) - min(correlated));
ampV = 0.5*(max(Vsignal) - min(Vsignal));
correlated = correlated*(ampV/ampCor); %to get the amplitudes of the two the same
offset = mean(correlated) - mean(Vsignal); 
correlated = correlated - offset; %to get the means alligned

figure;
plot(time, correlated, 'r', time, Vsignal, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('Resolved image vs time');
legend('Correlated signal', 'Vsignal');
