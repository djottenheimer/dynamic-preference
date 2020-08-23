function [LH, probSpike, V, mean_predictedSpikes, RPE] = AllModeled_mean(startValues, spikeCounts, included)


probSpike = poisspdf(spikeCounts(included), mean(spikeCounts(included))); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = mean(spikeCounts(included));

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;