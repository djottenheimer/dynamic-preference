function [LH, probSpike, V, mean_predictedSpikes, RPE] = AllModeled_pref(startValues, spikeCounts, rewards, included)

slope = startValues(1);
intercept = startValues(2);

%logistic curve
midpoint = startValues(3) * length(rewards); %x0, which trial is 50/50 preference
steepness = startValues(4); %k, how quickly does preference switch 

% initialize parameters
trials = length(rewards);
RPE = zeros(trials,1);
preference = zeros(trials,1);

%iterate
for t = 1:trials
    preference(t) = 1 / (1 + exp(-steepness*(t-midpoint)));
end

%flip preference for water
preference(rewards==0) = 1 - preference(rewards==0);

%normalize preference
preference = (preference - min(preference))/(max(preference) - min(preference));
V=preference;


%preference is in terms of sucrose, so need to reverese for water trials
RPE(rewards==0) = preference(rewards==0);
RPE(rewards==1) = preference(rewards==1);

rateParam = exp(slope*RPE + intercept);

    

% given a poisson process, what is the probability of each spike?
probSpike = poisspdf(spikeCounts(included), rateParam(included)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(included);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end
