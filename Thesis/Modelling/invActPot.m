function [ actionSignal ] = invActPot( action, dt )
%Takes an actionPotential from a group of m neurons and then tries to
%decode it as some sort of coherent voltage signal.
%   By gathering all the neuron data into one matrix we can look at some
%   column n and by the number of neurons with spikes at n we can infer the
%   probability of a spike at n which is proportional to the firing rate
%   which is dependent on the voltage.

prob_n = zeros(1, size(action, 2)); %empty vector to fill

%this loop will approximate the probability of a spike at point n
for n = 1:size(action, 2)
    actionTime = action(:, n); %all action potentials at time n
    spikes = size(find(actionTime == 1), 1); %# of spikes at time n
    activeNeurons = size(action, 1) - size(find(actionTime > 0), 1); %# of active neurons at time n
    prob_n(1, n) = spikes/activeNeurons; %prob of a spike at time n
end

firing_rate = prob_n/dt; %gives the firing rate in Hz
actionSignal = (11.5/(10^6))*log((1/0.9)*(1 - (62./(firing_rate - 1.6)))); 
                            %rearranging the gain fntn gives us the
                            %actionSignal as a fntn of the firing_rate

    

end

