function ms = helper_AllModeled(os, varargin)
% helper_RW_RPE    Fits neural models to test RW model predictions ~ RPE
%   ms = helper_RW_RPE(os, varargin)
%   INPUTS
%       os: behavioral data structure
%           .spikeCount: number of spikes in a particular time bin
%           .rewards: 0 for maltodextrin, 1 for sucrose
%           .Included: 1 if the animal licked within 2s, 0 if the animal did not (logical)
%               Included will always be more than or equal to the number of spike trials
%       varargin
%           StartingPoints: determines how many points to optimize from
%           ParticularModel: cell array of strings of models to use
%           RNG: random number generator seed (default = 1)
%   OUTPUTS
%       ms: model structure of fits




p = inputParser;
p.addParameter('StartingPoints', 1)
p.addParameter('ParticularModel', []);
p.addParameter('RNG', []);
p.parse(varargin{:});

if ~isempty(p.Results.RNG)
    rng(p.Results.RNG)
end
              
% Initialize models
if isempty(p.Results.ParticularModel)
    modelNames = {'mean','tv','pref','spp','licks'};
    %modelNames = {'mean','abre'};
else
    modelNames = p.Results.ParticularModel;
end


% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

% set boundary conditions
alpha_range = [0 1];
weighting_range = [0 1];
slope_range = [0 10]; % reward sensitivity, only positive (increased firing for higher value)
intercept_range = [-10 10]; % baseline spiking (converted with exp)
sucDecay_range = [0 0.05]; 
watDecay_range = [0 0.05]; %full decay in 20 trials 0.05
sucStart_range = [0 1];
steepness_range = [0 3];
midpoint_range = [-0.2 1.2];

%for divisive normalization
% baseline_range = [0 100];
% Rmax_range = [0 300];
sigma_range = [0 2];



for currMod = modelNames
    currMod = currMod{:};
    
    % initialize output variables
    runs = p.Results.StartingPoints;
    LH = zeros(runs, 1);
    
    if strcmp(currMod, 'mean')
        paramNames = {''};
        numParam = 0;        
        [LH, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            AllModeled_mean([], os.spikeCount, os.Included);
        bestFit = 1;
        ms.(currMod).bestParams = [];
        hess = [];
    elseif strcmp(currMod, 'licks')
        paramNames = {'slope','intercept'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2);
           -slope_range(1); -intercept_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :)] = ...
                fmincon(@AllModeled_licks, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Licks, os.Included);
        end
        [~, bestFit] = min(LH);
        ms.(currMod).bestParams = allParams(bestFit, :);
        [~, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            AllModeled_licks(ms.(currMod).bestParams, os.spikeCount, os.Licks, os.Included);
        
        % get the hessian for confidence intervals
        [~, ~, ~, ~, ~, ~, hess] = fmincon(@AllModeled_licks, allParams(bestFit,:), A, b, [], [], [], [], [], options, os.spikeCount, os.Licks, os.Included);
    elseif strcmp(currMod, 'tv')
        paramNames = {'slope','intercept'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2); 
           -slope_range(1); -intercept_range(1);];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :)] = ...
                fmincon(@AllModeled_tv, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, os.Included);
        end
        [~, bestFit] = min(LH);
        ms.(currMod).bestParams = allParams(bestFit, :);
        [~, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            AllModeled_tv(ms.(currMod).bestParams, os.spikeCount, os.Outcome, os.Included);
        
        % get the hessian for confidence intervals
        [~, ~, ~, ~, ~, ~, hess] = fmincon(@AllModeled_tv, allParams(bestFit,:), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, os.Included);
    elseif strcmp(currMod, 'pref')
        paramNames = {'slope','intercept','midpoint','steepness'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1) ...
                       rand(runs, 1)*diff(midpoint_range) + midpoint_range(1) ...
                       rand(runs, 1)*diff(steepness_range) + steepness_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2);  midpoint_range(2);  steepness_range(2);
           -slope_range(1); -intercept_range(1); -midpoint_range(1); -steepness_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :)] = ...
                fmincon(@AllModeled_pref, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, os.Included);
        end
        [~, bestFit] = min(LH);
        ms.(currMod).bestParams = allParams(bestFit, :);
        [~, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            AllModeled_pref(ms.(currMod).bestParams, os.spikeCount, os.Outcome, os.Included);
        
        % get the hessian for confidence intervals
        [~, ~, ~, ~, ~, ~, hess] = fmincon(@AllModeled_pref, allParams(bestFit,:), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, os.Included);
    elseif strcmp(currMod, 'spp')
        paramNames = {'slope','intercept','midpoint','steepness','weighting'};
        startValues = [rand(runs, 1)*diff(slope_range) + slope_range(1) ...
                       rand(runs, 1)*diff(intercept_range) + intercept_range(1) ...
                       rand(runs, 1)*diff(midpoint_range) + midpoint_range(1) ...
                       rand(runs, 1)*diff(steepness_range) + steepness_range(1) ...
                       rand(runs, 1)*diff(weighting_range) + weighting_range(1)];
        numParam = size(startValues, 2);
        allParams = zeros(runs, numParam);
        A=[eye(numParam); -eye(numParam)];
        b=[ slope_range(2);  intercept_range(2);  midpoint_range(2);  steepness_range(2);  weighting_range(2);
           -slope_range(1); -intercept_range(1); -midpoint_range(1); -steepness_range(1); -weighting_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :)] = ...
                fmincon(@AllModeled_spp, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, os.Included);
        end
        [~, bestFit] = min(LH);
        ms.(currMod).bestParams = allParams(bestFit, :);
        [~, ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
            AllModeled_spp(ms.(currMod).bestParams, os.spikeCount, os.Outcome, os.Included);
        
        % get the hessian for confidence intervals
        [~, ~, ~, ~, ~, ~, hess] = fmincon(@AllModeled_spp, allParams(bestFit,:), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, os.Included);
    else 
        error('RW model: Model name not found')
    end
    % for the model we're on, get
    ms.(currMod).paramNames = paramNames;
    % get log likelihood
    ms.(currMod).LH = -1 * LH(bestFit, :);
    % get the BIC
    ms.(currMod).BIC = log(length(os.spikeCount))*numParam  - 2*ms.(currMod).LH;
    % get the AIC
    ms.(currMod).AIC = 2*numParam - 2*ms.(currMod).LH;
    % get the AICc
    ms.(currMod).AICc = ms.(currMod).AIC + (2*numParam^2 + 2*numParam)/(length(os.spikeCount) - numParam - 1);
    
    % get confidence intervals
    ms.(currMod).CIvals = sqrt(diag(inv(hess)))'*1.96;
end