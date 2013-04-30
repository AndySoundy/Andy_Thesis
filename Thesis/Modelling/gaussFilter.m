function [ filter ] = gaussFilter( peak, bandwidth, length, fs, sigma )
%To try and get less leakage into the time domain this filter is going to
%involve a gaussian of max height == 1 that goes to zero outside the
%bandwidth. Gaussian form of exp(-(x - nu)^2/(2*sigma^2))
%   Because the gaussian goes to zero at exp(-x^2) you could choose to end
%   it outside the bandwidth, say twice the bandwidth for example.

nu = peak; %centering the gaussian at the peak value omega
freqVec = fs*linspace(0, 1, length); %need this to get a squared value
gaussPulse = exp(-(freqVec - nu).^2/(2*sigma^2)); %creates a gaussian over all frequency
cutoff = squareFilter(peak, bandwidth, length); %cuts out frequencies outside bandwidth

filter = gaussPulse.*cutoff; %now have a gaussian that == 0 outside bandwidth



end

