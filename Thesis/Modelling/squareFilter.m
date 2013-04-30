function [ filter ] = squareFilter( peak, bandwidthHz, length, dt )
%This is the simplest filter and just consists of a vector 1's centered at
%the peak frequency and zeros at bandwidth/2 on either side
df = 1/(dt*length);
bandwidth = bandwidthHz/df; %gets the bandwidth in bins not Hz
halfBand = round(bandwidth/2);
filter = zeros(1, length); %empty vector to fill 

if peak > halfBand
    filter((peak - halfBand):(peak + halfBand)) = ones(1, bandwidth+1);
else 
    filter(1:peak + halfBand) = ones(1, peak + halfBand);
end

end

