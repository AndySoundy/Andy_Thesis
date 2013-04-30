function [ noise ] = gaussNoise( sigma, mean, time )
%gaussNoise, This function accepts some scalar variance sigma and some
%vector time and from there gives a range of values in a gaussian
%distribution around the given mean.

noise = randn(1, size(time,2))*sigma + mean;


end

