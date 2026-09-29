% =========================================================================
% Per-animal / per-position processing of visual + photostim trials.
%
% For each field-of-view (F), animal (ani), and imaging position (pos):
%   1. Split raw trials into visual-stim-only and photostim-only trials.
%   2. Concatenate photostim + visual trials along the spatial-frequency
%      dimension so "condition 1" = photostim, conditions 2:17 = vis stims.
%   3. Compute stim-evoked responses (mean over the StimON window).
%   4. Identify ROI x direction combos that are significantly *excited*
%      by the laser (responsiveEachDirON), using a paired t-test against
%      baseline.
%   5. Among laser-responsive ROI/dir combos, identify ROIs that are
%      significantly *inhibited* (SigInhAD_NZON) by comparing photostim
%      vs. visual-stim response distributions.
%   6. Pool results across animals/positions and apply a response-range
%      threshold (RangeThr).
% =========================================================================

baselineWindow = 1:15;   % frames used to estimate pre-stim baseline

for F = 1:2

    nAni = Params.n_ani{F};

    % Only these three are actually needed after the ani/pos loop below
    % (they get cat()'d across ani/pos at the end of each field F).
    % Everything else computed inside the loop is scratch for that single
    % (ani,pos) iteration and does not need to be stored per-(ani,pos).
    StimsSF{F}                      = cell(nAni, max(Params.n_pos{F}));
    SigInhAD_NZON{F}                = cell(nAni, max(Params.n_pos{F}));
    VisStimTrialsSFCatMeanONResp{F} = cell(nAni, max(Params.n_pos{F}));

    for ani = 1:nAni
        for pos = 1:Params.n_pos{F}(ani)

            if isempty(data{F}{ani,pos})
                continue
            end

            raw = data{F}{ani,pos}.Raw;

            % --- Split raw trials into visual-only and photostim-only ---
            visTrials   = raw(:,:,2:9,:,:,:);
            photoTrials = squeeze(raw(:,:,1,1,:,:));

            % --- Fold direction x spatial-frequency dims together, then
            %     prepend photostim as "condition 1" ---
            [nROI, nFrame, nDir, nSF, nRep] = size(visTrials);
            visSFCat = reshape(visTrials, nROI, nFrame, nDir*nSF, nRep);

            photoVisSFCat = NaN(nROI, nFrame, 1 + nDir*nSF, size(photoTrials,3), size(photoTrials,4));
            photoVisSFCat(:,:,1,:,:)      = photoTrials;
            photoVisSFCat(:,:,2:end,:,:)  = visSFCat;

            % --- Stim-evoked responses (mean over StimON window) ---
            visStimSFCat  = squeeze(mean(photoVisSFCat(:,Params.StimONResponse,:,:,:), 2, 'omitnan'));
            visStimSFCat2 = squeeze(mean(visSFCat(:,Params.StimONResponse,:,:), 2, 'omitnan'));
            visStimSFCatMean = squeeze(mean(visStimSFCat, 4, 'omitnan'));

            StimsSF{F}{ani,pos} = squeeze(mean(visTrials(:,Params.StimONResponse,:,:,:), 2, 'omitnan'));

            % --- Laser-ON trials only (2nd condition along dim 4) ---
            visSFCatON = squeeze(visSFCat(:,:,:,2,:));

            base    = squeeze(mean(visSFCatON(:,baselineWindow,:,:), 2));
            stimOn  = squeeze(mean(visSFCatON(:,Params.StimONResponse,:,:), 2));
            baseRepAvg   = squeeze(mean(base, 3));
            stimOnRepAvg = squeeze(mean(stimOn, 3));

            [nROIon, nDirOn] = size(visSFCatON, [1 3]);
            isResponsiveDir = false(nROIon, nDirOn);

            for ROI = 1:nROIon
                for Dir = 1:nDirOn
                    isSig = ttest(squeeze(base(ROI,Dir,:)), squeeze(stimOn(ROI,Dir,:)), alphaVal);
                    if isSig && stimOnRepAvg(ROI,Dir) > baseRepAvg(ROI,Dir)
                        isResponsiveDir(ROI,Dir) = true;
                    end
                end
            end
            % --- Keep only responses from ROI/Dir combos flagged above,
            %     for use in the inhibition test below ---
            respOnly = NaN(size(visStimSFCat2));
            for ROI = 1:size(respOnly,1)
                for Dir = 1:size(respOnly,2)
                    if isResponsiveDir(ROI,Dir)
                        respOnly(ROI,Dir,:,:) = visStimSFCat2(ROI,Dir,:,:);
                    end
                end
            end

            % --- Test each ROI for laser-driven inhibition: compare the
            %     photostim condition (1) vs. visual-stim conditions (2)
            %     across the pooled Dir x Rep distribution ---
            respPermuted = permute(respOnly, [1 3 2 4]);
            respFlat = reshape(respPermuted, size(respPermuted,1), size(respPermuted,2), []);

            nROIresp = size(respFlat, 1);
            isInhibited = false(nROIresp, 1);
            for ROI = 1:nROIresp
                photoResp = respFlat(ROI,1,:);
                visResp   = respFlat(ROI,2,:);
                [h, ~] = ttest2(photoResp, visResp, 'Alpha', alphaVal);
                if h && ~isnan(mean(visResp,3,'omitnan')) && ...
                        mean(visResp,3,'omitnan') < mean(photoResp,3,'omitnan')
                    isInhibited(ROI) = true;
                end
            end
            SigInhAD_NZON{F}{ani,pos} = isInhibited;

            % --- Build the "mean" version restricted to responsive
            %     ROI/Dir combos (condition 1 = photostim, kept as-is) ---
            meanResp = visStimSFCatMean;
            meanOnResp = NaN(size(meanResp));
            meanOnResp(:,1,:) = meanResp(:,1,:);
            for ROI = 1:size(meanOnResp,1)
                for Dir = 2:size(meanOnResp,2)
                    if isResponsiveDir(ROI, Dir-1)
                        meanOnResp(ROI,Dir,:) = meanResp(ROI,Dir,:);
                    end
                end
            end
            VisStimTrialsSFCatMeanONResp{F}{ani,pos} = meanOnResp;

        end
    end

    % --- Pool across animals/positions for this field ---
    VisStimTrialsSFCatONRespAll{F} = cat(1, VisStimTrialsSFCatMeanONResp{F}{:});
    SigInhAD_NZONAll{F} = logical(cat(1, SigInhAD_NZON{F}{:}));

    InhAnyDir{F} = VisStimTrialsSFCatONRespAll{F}(SigInhAD_NZONAll{F}, :, :);

    StimsSFAll{F}     = cat(1, StimsSF{F}{:});
    StimsSFInhAll{F}  = StimsSFAll{F}(SigInhAD_NZONAll{F}, :, :, :, :);

    % --- Response-range threshold: keep ROIs whose (photostim-only)
    %     response range across directions exceeds rangeSizeResp ---
    nUnits = size(VisStimTrialsSFCatONRespAll{F}, 1);
    RangeThr{F} = false(nUnits, 1);
    for i = 1:nUnits
        if range(VisStimTrialsSFCatONRespAll{F}(i,:,1)) > rangeSizeResp
            RangeThr{F}(i) = true;
        end
    end

end