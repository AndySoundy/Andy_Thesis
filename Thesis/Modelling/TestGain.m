%Thesis, This is just to test the gain function got from Camperi:2007 to
%reproduce the graph that they produce

voltage = linspace(-100*10^-6, 100, 100*10^-6);
gain = 1.6 + 62./(1 + 0.9*exp(voltage*(10^6)/11.5));%gives the firing rate in Hz

plot(voltage, gain)
xlabel('Voltage (V)');
ylabel('Frequency (Hz)');
title('Gain function vs frequency');
