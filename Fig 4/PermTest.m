function out = PermTest(respTrials, varargin)
% 50%-vs-50% random split test for additive vs multiplicative effects
% Real regression on stimulus-averaged responses.
%
% respTrials : [nStim x 2 x nTrials]  (OFF/ON)
% stimVals   : [nStim x 1] stimulus regressor
%
% Optional:
%   nPerm : number of permutations (default 200)
%   alpha : significance level (default 0.05)

% -------------------------
% Parse inputs
% -------------------------
p = inputParser;
addParameter(p, 'nPerm', 200);
addParameter(p, 'alpha', 0.05);
parse(p, varargin{:});

nPerm = p.Results.nPerm;
alpha = p.Results.alpha;

PercentTrials=0.5;
rng('default')

[nStim, ~, nTrials] = size(respTrials);

% -------------------------
% Real regression on stimulus-averaged responses
% -------------------------
mean_OFF = squeeze(nanmean(respTrials(:,1,:),3));  % nStim x 1
mean_ON  = squeeze(nanmean(respTrials(:,2,:),3));  % nStim x 1

X_real = [ones(nStim,1) mean_OFF];  % intercept + slope
b_real = X_real \ mean_ON;

offset_real = b_real(1);  % additive
gain_real   = b_real(2);  % multiplicative

% -------------------------
% Null distribution via 50%-50% random trial splits
% -------------------------
gain_null   = zeros(nPerm,1);
offset_null = zeros(nPerm,1);

for pIdx = 1:nPerm
    meanA = zeros(nStim,1);
    meanB = zeros(nStim,1);

    for s = 1:nStim
        % Combine all trials for this stimulus (OFF + ON)
        trials_all = squeeze(respTrials(s,1,:));   % OFF trials
        trials_all = [trials_all; squeeze(respTrials(s,2,:))]; % ON trials
        nTotal = numel(trials_all);

        nSelect = round(PercentTrials*nTotal);

        permIdx = randperm(nTotal);

        idxA = permIdx(1:nSelect);
        idxB = permIdx(nSelect+1:2*nSelect);

        meanA(s) = nanmean(trials_all(idxA));
        meanB(s) = nanmean(trials_all(idxB));
    end

    % Regression across stimuli: groupB = intercept + slope * groupA
    Xperm = [ones(nStim,1) meanA];
    b_perm = Xperm \ meanB;

    gain_null(pIdx)   = b_perm(2);  % slope
    offset_null(pIdx) = b_perm(1);  % intercept
end

% -------------------------
% Significance testing using null distribution
% -------------------------
gain_CI   = prctile(gain_null,   [100*alpha/2 100*(1-alpha/2)]);
offset_CI = prctile(offset_null, [100*alpha/2 100*(1-alpha/2)]);

sig_mult = gain_real > gain_CI(2);
sig_div  = gain_real < gain_CI(1);

sig_add  = offset_real > offset_CI(2);  % positive additive
sig_sub  = offset_real < offset_CI(1);  % negative additive

% -------------------------
% Output
% -------------------------
out.gain_real    = gain_real;
out.offset_real  = offset_real;

out.gain_null    = gain_null;
out.offset_null  = offset_null;

out.gain_CI      = gain_CI;
out.offset_CI    = offset_CI;

out.sig_mult     = sig_mult;
out.sig_div      = sig_div;
out.sig_add      = sig_add;
out.sig_sub      = sig_sub;

end
