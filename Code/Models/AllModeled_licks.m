function [LH, probSpike, V, mean_predictedSpikes, RPE] = AllModeled_licks(startValues, spikeCounts, licks, included)
% firing rate only correlates with current reward; no real learning happens

slope = startValues(1);
intercept = startValues(2);

%normalize licks
licks=(licks-mean(licks))/std(licks);


rateParam = exp(slope*licks + intercept);
rateParam(rateParam < 0) = 0.1; % this is vestigial

probSpike = poisspdf(spikeCounts(included), rateParam(included)); % mask rateParam to exclude trials where the animal didn't lick fast enough

mean_predictedSpikes = rateParam(included);

if any(isinf(log(probSpike)))
    LH = 1e9;
else
    LH = -1 * sum(log(probSpike));
end

V = NaN;
RPE = NaN;