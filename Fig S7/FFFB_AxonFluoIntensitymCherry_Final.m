% ================================
%  Fluorescence vs Distance to mCherry+ Cell
%  MATLAB 2017 Compatible (RM-ANOVA block requires fitrm/ranova, R2014b+)
% ================================
%
%  INPUTS:  laminar_profiles.mat  (results, cell_counts, mouse_names,
%                                  slice_names, BIN_UM)
%
%  PIPELINE PER SLICE:
%   1. Load fluo profile and cell-count profile for the slice.
%   2. Smooth fluo profile with a moving average.
%   3. Find the bin of peak fluorescence; normalise profile to that peak.
%   4. Select mCherry+ cells within MAX_DIST_FROM_PEAK um of the peak bin.
%   5. For each selected cell, shift the x-axis so that cell = 0 um,
%      and extract HALF_RANGE_UM on each side (no mirroring).
%   6. Average per-cell curves within each slice.
%   7. Repeated-measures one-way ANOVA across slices (factor = distance),
%      over the same folded range used in the mean/SEM plot.
%   8. Average across slices and plot mean +/- SEM with shaded error band.

clear; close all; clc;

% ---- Parameters ----------------------------------------------------------
SMOOTH_WIN         = 20;    % moving-average window (bins)
MAX_DIST_FROM_PEAK = 150;   % um: only include mCherry cells within this
                             %     distance of the fluorescence peak bin
HALF_RANGE_UM      = 300;   % um: extract profile this far left AND right
                             %     of each mCherry cell
PLOT_COLOR         = [0 0.8 0.6];  % RGB: change plot colour here
ROOT_DIR           = 'D:\Confocal data spectral unmixing\Analysis of axon density around looped cells';

% ---- Load data -----------------------------------------------------------
fprintf('Select the laminar_profiles.mat file...\n');
[mat_file, mat_path] = uigetfile(fullfile(ROOT_DIR, '*.mat'), 'Select laminar_profiles.mat');
if isequal(mat_file, 0), disp('Cancelled.'); return; end
load(fullfile(mat_path, mat_file), ...
    'results', 'cell_counts', 'mouse_names', 'slice_names', 'BIN_UM');

HALF_RANGE_BINS = round(HALF_RANGE_UM / BIN_UM);   % half-window in bins

% ---- Create plots output folder ------------------------------------------
plots_dir = fullfile(ROOT_DIR, 'plots');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
fprintf('Figures will be saved to: %s\n', plots_dir);
fig_idx = 0;   % global figure counter for file naming
n_exp           = numel(results);
x_axis          = (-HALF_RANGE_BINS : HALF_RANGE_BINS) * BIN_UM;  % common x (um)
n_pts           = numel(x_axis);

% ---- Process each experiment separately ----------------------------------
for e = 1:n_exp
    fprintf('\n=== Experiment %d ===\n', e);

    n_mice = size(results{1,e}, 1);
    all_slice_curves = {};   % one entry per valid slice

    for m = 1:n_mice
        n_slices = size(results{1,e}, 2);

        for s = 1:n_slices
            fluo_raw = results{1,e}{m, s};
            cell_cnt = cell_counts{1,e}{m, s};

            if isempty(fluo_raw) || isempty(cell_cnt)
                continue;
            end

            n_bins   = numel(fluo_raw);
            bin_ctrs = ((0:n_bins-1) + 0.5) * BIN_UM;   % bin centres in um

            % ---- Step 1: smooth -----------------------------------------
            fluo_sm = movingAverage(fluo_raw, SMOOTH_WIN);

            % ---- Step 2: find peak and normalise ------------------------
            [peak_val, peak_bin] = max(fluo_sm);
            fluo_norm = fluo_sm / peak_val;
            peak_um   = bin_ctrs(peak_bin);

            % ---- Intermediate plot A: raw / smoothed / normalised -------
            figure;
            subplot(2,1,1); hold on;
            plot(bin_ctrs, fluo_raw, 'Color', [0.7 0.7 0.7], 'LineWidth', 1,   'DisplayName', 'Raw');
            plot(bin_ctrs, fluo_sm,  'Color', PLOT_COLOR,     'LineWidth', 1.8, 'DisplayName', 'Smoothed');
            xline2(peak_um, 'k--', 'Peak');
            xlabel('Position along pia (um)'); ylabel('Fluorescence - bg');
            title(sprintf('Exp%d | %s | slice %d  -  Profile', e, mouse_names{1,e}{m}, s), ...
                'Interpreter', 'none');
            legend('show'); grid on;

            subplot(2,1,2); hold on;
            plot(bin_ctrs, fluo_norm, 'Color', PLOT_COLOR, 'LineWidth', 1.8);
            xline2(peak_um, 'k--', 'Peak');
            xlabel('Position along pia (um)'); ylabel('Norm. fluorescence (peak=1)');
            title('Normalised'); grid on;
            fig_idx = fig_idx + 1; savePlot(gcf, plots_dir, fig_idx, sprintf('A_profile_exp%d_m%d_s%d', e, m, s));
close
            % ---- Step 3: find eligible mCherry cells --------------------
            % Expand counts into list of bin indices (one entry per cell)
            cell_bin_list = [];
            for b = 1:min(n_bins, numel(cell_cnt))
                if cell_cnt(b) > 0
                    cell_bin_list = [cell_bin_list, repmat(b, 1, cell_cnt(b))]; %#ok<AGROW>
                end
            end

            if isempty(cell_bin_list)
                fprintf('  Exp%d | %s | slice %d: no mCherry cells, skipping\n', ...
                    e, mouse_names{1,e}{m}, s);
                continue;
            end

            dist_from_peak = abs(bin_ctrs(cell_bin_list) - peak_um);
            eligible       = cell_bin_list(dist_from_peak <= MAX_DIST_FROM_PEAK);

            fprintf('  Exp%d | %s | slice %d: %d cells total, %d within %g um of peak\n', ...
                e, mouse_names{1,e}{m}, s, numel(cell_bin_list), numel(eligible), MAX_DIST_FROM_PEAK);

            if isempty(eligible)
                continue;
            end

            % ---- Intermediate plot B: cell selection --------------------
            figure; hold on;
            plot(bin_ctrs, fluo_norm, 'Color', PLOT_COLOR, 'LineWidth', 1.8);
            xline2(peak_um, 'k--', 'Peak');
            scatter(bin_ctrs(cell_bin_list), zeros(1,numel(cell_bin_list)), ...
                40, [0.7 0.7 0.7], 'v', 'filled', 'DisplayName', 'All cells');
            scatter(bin_ctrs(eligible), zeros(1,numel(eligible)), ...
                60, [0.8 0.2 0.2], 'v', 'filled', ...
                'DisplayName', sprintf('Within %g um', MAX_DIST_FROM_PEAK));
            plot(peak_um + [-MAX_DIST_FROM_PEAK MAX_DIST_FROM_PEAK], [0 0], ...
                'r-', 'LineWidth', 2);
            xlabel('Position along pia (um)'); ylabel('Norm. fluorescence');
            title(sprintf('Cell selection | Exp%d | %s | slice %d', ...
                e, mouse_names{1,e}{m}, s), 'Interpreter', 'none');
            legend('show'); grid on;
            fig_idx = fig_idx + 1; savePlot(gcf, plots_dir, fig_idx, sprintf('B_cellselect_exp%d_m%d_s%d', e, m, s));
close
            % ---- Step 4: extract per-cell shifted curves ---------------
            % For each cell at bin cb, extract fluo_norm over cb-HALF:cb+HALF.
            % x=0 = cell position. Bins outside profile range -> NaN.
            per_cell_full = NaN(numel(eligible), n_pts);   % full symmetric window

            for ci = 1:numel(eligible)
                cb = eligible(ci);
                for k = 1:n_pts
                    src_bin = cb - HALF_RANGE_BINS + (k-1);
                    if src_bin >= 1 && src_bin <= n_bins
                        per_cell_full(ci, k) = fluo_norm(src_bin);
                    end
                end
            end

            % ---- Step 4b: fold each curve (average both sides) ----------
            % x_axis runs from -HALF to +HALF. Centre index = HALF_RANGE_BINS+1.
            % For distance d (bins 1..HALF_RANGE_BINS): average index (ctr-d) and (ctr+d).
            % Result is a one-sided curve: x_fold = 0, BIN_UM, 2*BIN_UM, ...
            ctr          = HALF_RANGE_BINS + 1;             % centre index in full window
            x_fold       = (0:HALF_RANGE_BINS) * BIN_UM;   % one-sided x axis (um)
            per_cell_fold = NaN(numel(eligible), HALF_RANGE_BINS + 1);

            for ci = 1:numel(eligible)
                row = per_cell_full(ci, :);
                for d = 0:HALF_RANGE_BINS
                    left_val  = row(ctr - d);   % d bins left of cell
                    right_val = row(ctr + d);   % d bins right of cell
                    per_cell_fold(ci, d+1) = nanmean([left_val, right_val]);
                end
            end

            % ---- Intermediate plot C: full shifted curves ---------------
            figure; hold on;
            for ci = 1:size(per_cell_full, 1)
                plot(x_axis, per_cell_full(ci,:), ...
                    'Color', [PLOT_COLOR 0.3], 'LineWidth', 1);
            end
            xline2(0, 'k--', 'Cell');
            xlabel('Distance from mCherry+ cell (um)');
            ylabel('Norm. fluorescence');
            title(sprintf('Per-cell shifted curves | Exp%d | %s | slice %d', ...
                e, mouse_names{1,e}{m}, s), 'Interpreter', 'none');
            grid on;
            fig_idx = fig_idx + 1; savePlot(gcf, plots_dir, fig_idx, sprintf('C_percell_exp%d_m%d_s%d', e, m, s));

            % ---- Intermediate plot D: folded (two-sides averaged) -------
            figure; hold on;
            for ci = 1:size(per_cell_fold, 1)
                plot(x_fold, per_cell_fold(ci,:), ...
                    'Color', [PLOT_COLOR 0.3], 'LineWidth', 1);
            end
            slice_mean = nanmean(per_cell_fold, 1);
            plot(x_fold, slice_mean, 'Color', PLOT_COLOR, 'LineWidth', 2.5, ...
                'DisplayName', 'Slice mean');
            xline2(0, 'k--', 'Cell');
            xlabel('Distance from mCherry+ cell (um)');
            ylabel('Norm. fluorescence');
            title(sprintf('Per-cell folded curves | Exp%d | %s | slice %d', ...
                e, mouse_names{1,e}{m}, s), 'Interpreter', 'none');
            grid on;
            fig_idx = fig_idx + 1; savePlot(gcf, plots_dir, fig_idx, sprintf('D_folded_exp%d_m%d_s%d', e, m, s));

            all_slice_curves{end+1} = slice_mean; %#ok<AGROW>
close
        end % slices
    end % mice

    if isempty(all_slice_curves)
        fprintf('  No valid slices for Exp%d.\n', e);
        continue;
    end

    % ---- Step 5: assemble slices x distance matrix -----------------------
    % Same matrix / same range used below for the mean+SEM plot.
    curve_mat  = NaN(numel(all_slice_curves), HALF_RANGE_BINS + 1);
    for k = 1:numel(all_slice_curves)
        L = numel(all_slice_curves{k});
        curve_mat(k, 1:L) = all_slice_curves{k};
    end
    x_fold_plot = (0:HALF_RANGE_BINS) * BIN_UM;

    % ---- Step 6: repeated measures one-way ANOVA (factor = Distance) -----
    % n=1 unit = one slice's folded mean curve. Repeated factor = distance
    % bin, over exactly the range plotted in the mean/SEM figure below.
    validRows = all(~isnan(curve_mat), 2);
    if any(~validRows)
        fprintf('  Excluding %d/%d slices from RM-ANOVA (incomplete data over full range).\n', ...
            sum(~validRows), numel(validRows));
    end
    curve_mat_anova = curve_mat(validRows, :);

    if size(curve_mat_anova, 1) < 3
        warning('Exp%d: fewer than 3 complete slices, skipping RM-ANOVA.', e);
    else
        nBinsFold = size(curve_mat_anova, 2);
        varNames  = arrayfun(@(k) sprintf('d%d', k), 1:nBinsFold, 'UniformOutput', false);
        Tanova    = array2table(curve_mat_anova, 'VariableNames', varNames);

        withinDesign = table(x_fold_plot(:), 'VariableNames', {'Distance'});
        withinDesign.Distance = categorical(withinDesign.Distance);

        rm = fitrm(Tanova, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), ...
            'WithinDesign', withinDesign);
        ranovaTbl = ranova(rm, 'WithinModel', 'Distance');

        fprintf('\n--- Exp%d: Repeated measures one-way ANOVA (factor = Distance, n = %d slices) ---\n', ...
            e, size(curve_mat_anova, 1));
        disp(ranovaTbl);
effectRow = '(Intercept):Distance';   % change to '(Intercept):Position' etc. if that's your factor name
errorRow  = 'Error(Distance)';

F_stat = ranovaTbl.F(effectRow);
df1    = ranovaTbl.DF(effectRow);
df2    = ranovaTbl.DF(errorRow);
p_val  = ranovaTbl.pValueGG(effectRow);   % Greenhouse-Geisser corrected, recommended with many levels

disp(sprintf('RM one-way ANOVA: F(%.0f,%.0f) = %.2f, p = %.3g', df1, df2, F_stat, p_val));

        try
            sphericityTbl = mauchly(rm);
            fprintf('Mauchly sphericity test (Exp%d):\n', e);
            disp(sphericityTbl);
        catch ME
            warning('Exp%d: could not compute sphericity test: %s', e, ME.message);
        end
    end

    % ---- Step 7: grand mean +/- SEM across slices -------------------------
    grand_mean = nanmean(curve_mat, 1);
    grand_sem  = nanstd(curve_mat, 0, 1) ./ sqrt(sum(~isnan(curve_mat), 1));

    % ---- Final plot: mean +/- SEM ---------------------------------------
    figure; hold on;
    fill([x_fold_plot, fliplr(x_fold_plot)], ...
         [grand_mean + grand_sem, fliplr(grand_mean - grand_sem)], ...
         PLOT_COLOR, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(x_fold_plot, grand_mean, 'Color', PLOT_COLOR, 'LineWidth', 2.5);
    xline2(0, 'k--', 'mCherry cell');
    xlabel('Distance from mCherry+ cell (um)');
    ylabel('Normalised fluorescence (peak = 1)');
    title(sprintf('Experiment %d  -  Mean +/- SEM  (n = %d slices)', ...
        e, numel(all_slice_curves)));
    grid on;
    fig_idx = fig_idx + 1; savePlot(gcf, plots_dir, fig_idx, sprintf('E_grandmean_exp%d', e));
close
    %% plot for fig
    figure; hold on;
    fill([x_fold_plot, fliplr(x_fold_plot)], ...
         [grand_mean + grand_sem, fliplr(grand_mean - grand_sem)], ...
         PLOT_COLOR, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(x_fold_plot, grand_mean, 'Color', PLOT_COLOR, 'LineWidth', 1);
    ylim([0 1])
xlabel('distance from V1_{PM} neuron (um)'); ylabel('tdTomato Fluorescence Intensity (norm)');
xlabel('distance from PM_{V1} neuron (um)'); ylabel('tdTomato Fluorescence Intensity (norm)');

%     title(sprintf('Experiment %d  -  Mean +/- SEM  (n = %d slices)', ...
%         e, numel(all_slice_curves)));
    ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
% despline(gca);
set(gcf,'Units','centimeters','Position',[10 10 4.5 3.5]);
set(gca,'FontSize',7,'LineWidth',0.8,'TickDir','out','TickLength',[0.01, 0.01]);

    fig_idx = fig_idx + 1; savePlot(gcf, plots_dir, fig_idx, sprintf('E_grandmean_exp%d', e));
    
saveas(gcf,[plots_dir filesep num2str(e) ' grandmean FF' ],'svg');
    
end % experiments

fprintf('\nDone.\n');


%% ---- Local Functions ----------------------------------------------------

function y = movingAverage(x, win)
% MOVINGAVERAGE  Centred moving average; edges use available samples only.
    y = x;
    n = numel(x);
    hw = floor(win / 2);
    for i = 1:n
        i1 = max(1, i - hw);
        i2 = min(n, i + hw);
        y(i) = mean(x(i1:i2));
    end
end

function xline2(x_val, style, lbl)
% XLINE2  MATLAB-2017-compatible vertical reference line.
    yl = ylim;
    plot([x_val x_val], yl, style, 'LineWidth', 1.2);
    text(x_val, yl(2), lbl, 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

function savePlot(fig, plots_dir, idx, name)
% SAVEPLOT  Save figure as PNG into plots_dir.
    fname = fullfile(plots_dir, sprintf('%03d_%s.png', idx, name));
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, fname, '-dpng', '-r150');
    fprintf('    Saved: %s\n', fname);
end