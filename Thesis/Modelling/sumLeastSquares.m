function [ sumRsquared ] = sumLeastSquares( numNeurons, Vsignal, Msignal, dt,...
    vElec, relaxation)
%leastSquares. This fntn file is designed to take the inputs, spit out a
%correlated signal and then provide sqrt(residual^2), i.e. square root of
%the sum of least squares. Hopefully this fntn file can be used with some
%loop to give a plot of fit versus number of neurons
%   Detailed explanation goes here
length = size(Vsignal,2);
ampV = 0.5*(max(Vsignal) - min(Vsignal));
ampM = 0.5*(max(Msignal) - min(Msignal));

%% Converting Vsignal into an actionPotential and back
actionpotentials = zeros(numNeurons, length); %starting matrix for all neuron inputs

%the loop below creates an mxlength matrix with entries from all neurons
for m = 1:numNeurons
    actionpotentials(m, :) = actionPot(vElec, relaxation, dt); %each row is one neuron
end

actionSignal = invActPot(actionpotentials, dt); %what the sharks brain resolves


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
correlated = (length^2)*ifft((ampV/ampM)*(filtered_V));


%% Least squares comparison
%When this was written I was having trouble in that the correlated signal
%was getting a DC offset that would have thrown out the least squares
%analysis so I manually get them to centre at the same DC point, obviously
%not ideal but will try and work it out

% middle = max(vElec) - abs(min(vElec)); %getting the centre of vElec
% offsetCorr = 0.5*(max(correlated) - abs(min(correlated))) - middle;
% correlated = correlated - offsetCorr;
% ampCor = 0.5*(max(correlated) - min(correlated));
% correlated = correlated*ampV/ampCor;
r = vElec - real(correlated);
sumRsquared = sum(r.^2)/length; %sum of square of differences



end

