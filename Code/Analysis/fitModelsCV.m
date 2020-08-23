clear; clc
load('NeuronInfo.mat');
data = NeuronInfo;

nStart = 20;
RNG_val = 1;
%%
clear os_temp os

all_fits = struct(); % initialize an empty structure 
for ind = 1:length(data) % for all neurons
    fprintf('n %i of %i\n', ind, length(data))
    
    %put data into os
    os_temp(ind).Rat = data(ind).Rat;
    os_temp(ind).Ses = data(ind).Ses;
    os_temp(ind).Typ = data(ind).Typ;
    os_temp(ind).Outcome = data(ind).Outcome;
    os_temp(ind).Choice = data(ind).Choice;
    os_temp(ind).Licks = data(ind).Licks;
    os_temp(ind).InPort = data(ind).InPort;
    os_temp(ind).Window = data(ind).Window;
    os_temp(ind).Evs = data(ind).Evs;
        
        
 
    % get spike counts
    os_temp(ind).spikeCount_RD = data(ind).Spikes(:,6); 
    
    %pick included trials -- forced, and in port
    os_temp(ind).Included = data(ind).Choice==0 & data(ind).InPort==1;
    
    % spikeCount is a temporary field
    % fit RD
    os_temp(ind).spikeCount = os_temp(ind).spikeCount_RD;
    ms = helper_AllModeled_CV(os_temp(ind), 'StartingPoints', nStart, 'RNG', RNG_val);
    os_temp(ind).mod_RD = ms;
    
    % remove spikeCount to avoid future confusion
    os(ind) = rmfield(os_temp(ind), 'spikeCount');
end


save('ModelFits_CV.mat','os')

fprintf('Finished\n')