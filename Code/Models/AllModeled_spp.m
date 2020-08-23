function [LH, probSpike, V, mean_predictedSpikes, RPE] = AllModeled_spp(startValues, spikeCounts, rewards, included)

slope = startValues(1);
intercept = startValues(2);

%logistic curve
midpoint = startValues(3) * length(rewards); %Xo, which trial is 50/50 preference. Inputted as fraction of session
steepness = startValues(4); %k, how quickly does preference switch 
weight = startValues(5); %0 is pure set point, 1 is pure rel

% initialize parameters
trials = length(rewards);
V = zeros(trials + 1,1);
setpoint = zeros(trials + 1,1);%set point
RPE = zeros(trials,1);
preference = zeros(trials,1);

setpoint(1) = 1;
%iterate
for t = 1:trials
    setpoint(t+1) = 1-t/trials;
    preference(t) = 1 / (1 + exp(-steepness*(t-midpoint)));
end

setpoint=setpoint(1:trials);

%flip preference for water
preference(rewards==0) = 1 - preference(rewards==0);

%normalize preference
preference = (preference - min(preference))/(max(preference) - min(preference));
V=preference;
    

%perform weighting
RPE(rewards==0) = weight*(preference(rewards==0)) + (1 - weight)*setpoint(rewards==0);
RPE(rewards==1) = weight*(preference(rewards==1)) + (1 - weight)*setpoint(rewards==1);

rateParam = exp(slope*RPE + intercept);

    

% given a poisson process, what is the probability of each spike?
probSpike = poisspdf(spikeCounts(included), rateParam(included)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(included);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end
