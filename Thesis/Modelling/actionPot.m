function [ actionPotential ] = actionPot( Vsignal, relaxation, dt )
%This will take the Vsignal input as the potential experienced at the
%individual neuron and then translate it into some series of action
%potentials which will be characterised by the gain function and the
%relaxation time of the neurons.
%The gain function used is = 1.6 + 62/(1 + 0.9*exp((10^6)*Vsignal/11.5))
%Going to try and use the gain fntn to get the spacing between spikes in s
%and then create a firing rate vector composed of 1's and 0's to represent
%the action potential of this neuron

actionPotential = zeros(1, size(Vsignal,2)); %empty vector atm
firing_rate = 1.6 + 62./(1 + 0.9*exp((10^6)*Vsignal/11.5)); %gives the firing rate in Hz
relaxation_index = round(relaxation/dt); %gives the time in dt's of the relaxation time
random = rand(1, size(Vsignal,2)); %creates a vector of #'s between 0 and 1
n = 1; %sets the index for the while loop
relax_time = linspace(0, 1, relaxation_index); 
relax_decay = exp(-relax_time); %creates an exponential decay during the relaxation time

while n <= size(Vsignal,2)
    prob_n = firing_rate(1,n)*dt; %probability of a spike at n given Vsignal
    
    if random(1,n) > prob_n
        actionPotential(1,n) = 0;
        n = n + 1; %increments n to the location of the latest spike
        
    else actionPotential(1, n) = 1; %spike w/ relaxation time
        if (n + 1 + relaxation_index) <=  size(Vsignal,2)
            actionPotential(1, (n):(n + relaxation_index - 1)) = relax_decay;
            n = (n + 1 + relaxation_index); %increments n to the location of the latest spike
        
        else n = n + 1; %just increments n by 1
        end
        
    end
end
    
                                                    

end

