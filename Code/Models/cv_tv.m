function [LH, probSpike, V, mean_predictedSpikes, RPE] = cv_tv(startValues, spikeCounts, rewards, train, test)

slope = startValues(1);
intercept = startValues(2);

% initialize parameters
trials = length(rewards);
%V = zeros(trials + 1, 1);
taskV = zeros(trials + 1,1);
V = zeros(trials,1);
RPE = zeros(trials,1);

taskV(1) = 1;

%V(1) = Vinit;
% Call learning rule
for t = 1:trials
        
    taskV(t + 1) = taskV(t) - 1/trials;
      
end

RPE=taskV(1:end-1);

rateParam = exp(slope*RPE + intercept);

% given a poisson process, what is the probability of each spike?
probSpike = poisspdf(spikeCounts(test), rateParam(train)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(test);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = taskV(1:end-1);
