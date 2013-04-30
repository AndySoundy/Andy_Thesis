%Thesis, Just testing that the matrix form of the action potentials of all
%the neurons combine and work in some vague sort of way
close all;
clear all;
clc;

%% Constants and vector initial setup
%Constants
dt = 0.0005; %sampling time for the neuron
length = 8192; %length of sampling period
lambdaV = 0.8; %proportionality constant related to sensitivity for V
veloc = 1; %forward velocity of the shark in m/s
omega = 2*pi; %2pi*vestibular frequency in Hz
headAmp = 0.05; %max amplitude of the head movement in m
relaxation = 0.001; %relaxation time in s
B = 40*10^-6; %strength of the magnetic field in Tesla
numNeurons = 30000; %number of neurons
sigmaV = 0.1; %std dev of noise in Voltage
bandwidth = 0.05; %bandwidth in Hz


%Vectors
time = (0:length-1)*dt; %time vector in s
accel = -(headAmp*(omega^2))*sin(omega*time); %head acceleration
psi = -(1/(omega^2))*accel; %psi gives the angle between v and B
vElec = lambdaV*veloc*B*sin(psi); %the potential the shark measures w/ no noise
accel = -(headAmp*(omega^2))*sin(omega*time); %head acceleration
modNoiseV = (max(vElec) - min(vElec)); %the order of magnitude of noiseV
noiseV = modNoiseV*gaussNoise(sigmaV, 0, time); %the noise in vSignal
Vsignal = noiseV + vElec; %what the shark measures

%% Getting the action potential matrix of many neurons
actionpotentials = zeros(numNeurons, length);

for m = 1:numNeurons
    actionpotentials(m, :) = actionPot(vElec, relaxation, dt); %each row is one neuron
end

figure(1)
actionSignal = invActPot(actionpotentials, dt);
subplot(3,1,1);
plot(time, vElec, 'r')
xlabel('Time (s)');
ylabel('Voltage');
title('Original signal minus noise')
subplot(3,1,2);
plot(time, Vsignal, 'b');
xlabel('Time (s)');
ylabel('Voltage');
title('Original signal with noise')
subplot(3,1,3);
plot(time, actionSignal, 'k');
xlabel('Time (s)');
ylabel('Voltage');
title('actionSignal')

%% Filter Testing Area
%Here can test different types of filtervec to see which one will give the
%better fit to Vsignal

F_actionSignal = fft(actionSignal)/length;
F_accel = fft(accel)/length;
peak = find(abs(F_accel) == max(abs(F_accel)),1); %index of peak of accel
bandwidth = 10; %at the moment in df but could upgrade to Hz relatively easily

%here we're going with the square filter but could replace that later if
%neccessary
Square = squareFilter(peak, 1, length, dt);
Filter = Square.*F_accel; %just leaves you with a vector with a spike at omega
                          %and zeros elsewhere
FilteredSig = Filter.*F_actionSignal;
correlated = ifft(FilteredSig)*length;

%to get the amplitudes of the two the same
ampCor = 0.5*(max(correlated) - min(correlated));
ampV = 0.5*(max(Vsignal) - min(Vsignal));
correlated = correlated*(ampV/ampCor)*lambdaV;

%to get the start to more or less line up with zero


figure(3)
subplot(3,1,1);
plot(time, vElec, 'r')
xlabel('Time (s)');
ylabel('Voltage');
title('Original signal minus noise')
subplot(3,1,2);
plot(time, Vsignal, 'b');
xlabel('Time (s)');
ylabel('Voltage');
title('Original signal with noise')
subplot(3,1,3);
plot(time, correlated, 'k');
xlabel('Time (s)');
ylabel('Voltage');
title('Correlated Signal')


%% Getting the parameters in the end
% 
% params = sinefit(actionSignal, time, omega/(2*pi));
% constantEst = real(params(1,1)); %offset from y = 0
% amplitudeEst = params(1,2);
% freqEst = params(1,3); %in Hz
% thetaEst = params(1,4); %in rad, from theta = 0 NOT THE SHARK HEADING!!!
% 
% sineEst = constantEst + amplitudeEst*sin(2*pi*freqEst*time + thetaEst);
% 
% plot(time, sineEst, 'b') 
% hold off


%% Graph of fit vs number of neurons
%Will get the sum of least squares for the correlated fit vs vElec (not
%Vsignal) and then graph that against the number of neurons required for
%that fit

% startNum = 50; %these three parameters refer to number of neurons in test
% stepNum = 50;
% maxNum = 2000;
% sumRsquares = zeros(1,maxNum/stepNum);
% Numvec = startNum:stepNum:maxNum;
% 
% for num = startNum:stepNum:maxNum;
%     sumRsquares(1,num/stepNum) = sumLeastSquares(num, Vsignal,accel,dt,vElec,relaxation);
% end
% 
% plot(Numvec, sumRsquares)





