function [LH, probSpike, V, mean_predictedSpikes, RPE] = cv_mean(startValues, spikeCounts, train, test)


probSpike = poisspdf(spikeCounts(test), mean(spikeCounts(train))); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = mean(spikeCounts(train));

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;