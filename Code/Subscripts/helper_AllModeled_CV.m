function ms = helper_AllModeled_CV(os, varargin)
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
    %modelNames = {'mean','outW','outS','tv','pref','spp','licks'};
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

%cross-validation parameters
train_fraction = 0.8;
num_repeats = 50;

%for divisive normalization
% baseline_range = [0 100];
% Rmax_range = [0 300];
sigma_range = [0 2];

%make partitions so all models are modeled with same sets of partitions
trains={};trains={};
for repeat = 1:num_repeats
    %generate new train and test groups
    trials = find(os.Included);
    groups = cat(1,ones(floor(train_fraction*length(trials)),1),zeros(ceil((1-train_fraction)*length(trials)),1));
    groups = groups(randperm(length(groups)));

    %pick training set
    train = zeros(length(os.Included),1);
    train(trials(groups==1))=1;
    train=logical(train);
    trains{repeat,1}=train;

    %pick test set
    test = zeros(length(os.Included),1);
    test(trials(groups==0))=1;
    test=logical(test);
    tests{repeat,1}=test;
end


for currMod = modelNames
    currMod = currMod{:};
   
    
    LHbest = NaN(num_repeats,1);
    for repeat = 1:num_repeats
        train = trains{repeat,1};
        test = tests{repeat,1};
              
        % initialize output variables
        runs = p.Results.StartingPoints;
        LH = zeros(runs, 1);
        
        
        if strcmp(currMod, 'mean')
            paramNames = {''};
            numParam = 0;
            [LHbest(repeat,1), ms.(currMod).probSpike, ms.(currMod).V, ms.(currMod).mean_predictedSpikes, ms.(currMod).RPEs] = ...
                cv_mean([], os.spikeCount, train, test);
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
                    fmincon(@AllModeled_licks, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Licks, train);
            end
            [~, bestFit] = min(LH);
            ms.(currMod).bestParams = allParams(bestFit, :);
            
            LHbest(repeat,1) = AllModeled_licks(ms.(currMod).bestParams, os.spikeCount, os.Licks, test);
            
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
                    fmincon(@AllModeled_tv, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, train);
            end
            [~, bestFit] = min(LH);
            ms.(currMod).bestParams = allParams(bestFit, :);
            LHbest(repeat,1) = AllModeled_tv(ms.(currMod).bestParams, os.spikeCount, os.Outcome, test);
            
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
                    fmincon(@AllModeled_pref, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, train);
            end
            [~, bestFit] = min(LH);
            ms.(currMod).bestParams = allParams(bestFit, :);
            LHbest(repeat,1) = AllModeled_pref(ms.(currMod).bestParams, os.spikeCount, os.Outcome, test);
            
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
                    fmincon(@AllModeled_spp, startValues(r, :), A, b, [], [], [], [], [], options, os.spikeCount, os.Outcome, train);
            end
            [~, bestFit] = min(LH);
            ms.(currMod).bestParams = allParams(bestFit, :);
            LHbest(repeat,1) = AllModeled_spp(ms.(currMod).bestParams, os.spikeCount, os.Outcome, test);
            
        else
            error('RW model: Model name not found')
        end
    end
    % for the model we're on, get
    ms.(currMod).paramNames = paramNames;
    % get cross-validated log likelihood
    ms.(currMod).CVLH = -1 * mean(LHbest);
    
end