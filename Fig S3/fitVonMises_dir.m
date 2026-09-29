function fit = fitVonMises_dir(theta, resp)
%
% Parameters:
% p(1) = baseline
% p(2) = amplitude 1  (dominant peak, enforced A1 >= A2 post-fit)
% p(3) = amplitude 2
% p(4) = preferred direction
% p(5) = kappa (width)

theta = double(theta(:));
resp  = double(resp(:));

% -------------------------
% Initial guesses
% -------------------------
b0 = min(resp);
A0 = max(resp) - min(resp);
[~,i0] = max(resp);
mu0 = theta(i0);
k0  = 1;
p0 = [b0, A0, A0/2, mu0, k0];

% kappa constraints
minFWHM = deg2rad(45);
maxFWHM = deg2rad(180);
kappa_from_fwhm = @(fwhm) log(2) ./ (1 - cos(fwhm/2));
kappa_max = kappa_from_fwhm(minFWHM);  % narrowest allowed
kappa_min = kappa_from_fwhm(maxFWHM);  % broadest allowed

% -------------------------
% Bounds
% -------------------------
lb = [0, 0, 0, 0,        kappa_min]; 
ub = [ 1, Inf, Inf, 2*pi,   kappa_max];

opts = optimset('Display','off');

% -------------------------
% Model
% -------------------------
vmDouble = @(p,x) ...
    p(1) + ...
    p(2)*exp(p(5)*(cos(x - p(4))      - 1)) + ...
    p(3)*exp(p(5)*(cos(x - p(4) - pi) - 1));

% -------------------------
% Fit — try multiple mu0 starting points and keep best
% -------------------------
mu0_candidates = linspace(0, 2*pi, 9);
mu0_candidates(end) = [];   

bestSSR = Inf;
pHat    = p0;   % fallback

for ci = 1:numel(mu0_candidates)
    p0c    = [b0, A0, A0/2, mu0_candidates(ci), k0];
    try
        pTry   = lsqcurvefit(vmDouble, p0c, theta, resp, lb, ub, opts);
        resid  = resp - vmDouble(pTry, theta);
        ssr    = sum(resid.^2);
        if ssr < bestSSR
            bestSSR = ssr;
            pHat    = pTry;
        end
    catch
        continue
    end
end

% -------------------------
% Store results
% -------------------------
fit.b      = pHat(1);
fit.A1     = pHat(2);
fit.A2     = pHat(3);
fit.mu     = mod(pHat(4), 2*pi);   % direction in [0, 2*pi)
fit.kappa  = pHat(5);

% -------------------------
% Break A1/A2 degeneracy: enforce A1 >= A2
% (A1,A2,mu) and (A2,A1,mu+pi) produce identical curves;
% standardise so A1 is always the dominant peak
% -------------------------
if fit.A2 > fit.A1
    [fit.A1, fit.A2] = deal(fit.A2, fit.A1);
    fit.mu = mod(fit.mu + pi, 2*pi);
end

% Amplitude summaries
fit.Amean  = (fit.A1 + fit.A2) / 2;
fit.Asym   = (fit.A1 - fit.A2) / (fit.A1 + fit.A2 + eps);   % in [0,1], >=0 by construction

% -------------------------
% Fine curve for plotting
% -------------------------
thFine = linspace(0, 2*pi, 360);
fit.thetaFine = thFine;
fit.curve     = vmDouble(pHat, thFine);
fit.respMean  = resp;
fit.theta     = theta;

% Peak response amplitude above baseline
fit.peakResp  = max(fit.curve) - fit.b;

fit.oppResp = pHat(2)*exp(-2*pHat(5)) + pHat(3);

end