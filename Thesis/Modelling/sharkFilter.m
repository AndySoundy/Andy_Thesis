function [ output_args ] = sharkFilter( actSignal, Msignal, dt )
%Here the user will input the action signal (or the Vsignal) and the
%vesitbular signal (Msignal), then we take the DFT of both signals and
%essentially use the DFT of Msignal to throw away frequencies that aren't
%integer multiples of the vestibular frequency (omega).
%   Detailed explanation goes here

fs = 1/dt; %sampling frequency
length = size(actSignal, 2); 
NFFT = 2^nextpow2(length); %next power of 2 from length of signal
F_actSignal = fft(actSignal, NFFT)/length; 
F_Msignal = fft(Msignal, NFFT)/length;

%Trying to get a filter that will pick out the frequencies of omega and
%2*omega
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
filtered_V = P_filtervec.*(F_actSignal); %filters actSignal

end

