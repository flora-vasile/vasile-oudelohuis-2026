saveDir_FigS5 = [Params.MainDir filesep 'plot\Fig S5'];
if ~exist(saveDir_FigS5,'dir')
    mkdir(saveDir_FigS5)
end

%FFFB_FitExc_Final

% --- OPTIONS ---
nBins       = 3;
showN       = false;   % write number of values per bin on plot
% ---------------

for F=1:2
    
    for ani=1:Params.n_ani{F}
        for pos=1:Params.n_pos{F}(ani)
            if ~isempty(data{F}{ani,pos})
                CellStat{F}{ani,pos}.SigExcAD=logical(CellStat{F}{ani,pos}.SigExcAD);
                VisRespExcSes{F}{ani,pos}=squeeze(nanmean(nanmean(data{F}{ani,pos}.StimAv(CellStat{F}{ani,pos}.SigExcAD,Params.StimONResponse,1,:),2),4)); 
                
                AllDirRespAllRepSes{F}{ani,pos}=squeeze(nanmean(data{F}{ani,pos}.StimAv(:,Params.StimONResponse,:,:),2));
                AllDirRespSes{F}{ani,pos}=squeeze(nanmean(nanmean(data{F}{ani,pos}.StimAv(:,Params.StimONResponse,:,:),2),4));
                
                dprimeSes{F}{ani,pos}=(AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.SigExcAD,2)...
                    -AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.SigExcAD,1))./...
                    sqrt(1/2*(nanvar(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.SigExcAD,2,:),0,3)...
                    +nanvar(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.SigExcAD,1,:),0,3)));
                                excitedAD{F}{ani,pos}=logical(CellStat{F}{ani,pos}.SigExcAD);

            end
        end
    end

    excited{F} = logical(cat(1,excitedAD{F}{:}));
    VisRespExc{F}=cat(1,VisRespExcSes{F}{:});
    dprimeExc{F}=cat(1,dprimeSes{F}{:});
    temp=squeeze(VisStimSFCatTrialMeanAll{F}(excited{F},:,1));
    RangeRespExc{F}=range(temp,2);

   
end

% Variables to plot on y-axis
yVars   = {'slopeExc', 'interceptExc'};
yLabels = {'Slope', 'Offset'};
% X-axis variables
xVars   = {'RangeRespExc', 'VisRespExc', 'dprimeExc'};
xLabels = {'Range of Visual Response (Ev./F0)', 'Mean Visual Response (Ev./F0)', 'dprime'};

%% --- ANCOVA: FF vs FB, controlling for covariate ---
ancovaResults = struct();

for xV = 1:3
    for yV = 1:2
        xData_all = [];
        yData_all = [];
        group_all = [];

        for F = 1:2
            xData = eval([xVars{xV} '{' num2str(F) '}']);
            yData = eval([yVars{yV} '{' num2str(F) '}']);
            valid = ~isnan(xData) & ~isnan(yData);

            xData_all = [xData_all; xData(valid)];
            yData_all = [yData_all; yData(valid)];
            group_all = [group_all; repmat(F, sum(valid), 1)]; % 1=FB, 2=FF
        end

        tbl = table(xData_all, yData_all, categorical(group_all), ...
                     'VariableNames', {'X','Y','Group'});

        % 1) Slope-homogeneity check: does the X-Y relationship differ by group?
        mdl_int   = fitlm(tbl, 'Y ~ X*Group');
        coefNames = mdl_int.CoefficientNames;
        intIdx    = find(contains(coefNames, ':'));            % interaction term
        t_int     = mdl_int.Coefficients.tStat(intIdx);
        p_int     = mdl_int.Coefficients.pValue(intIdx);

        % 2) Main ANCOVA: group difference controlling for X
        mdl_main  = fitlm(tbl, 'Y ~ X + Group');
        coefNames2 = mdl_main.CoefficientNames;
        grpIdx    = find(contains(coefNames2, 'Group') & ~contains(coefNames2, ':'));
        t_grp     = mdl_main.Coefficients.tStat(grpIdx);
        p_grp     = mdl_main.Coefficients.pValue(grpIdx);
        df        = mdl_main.DFE;

        ancovaResults(xV,yV).xVar          = xVars{xV};
        ancovaResults(xV,yV).yVar          = yVars{yV};
        ancovaResults(xV,yV).t_group       = t_grp;
        ancovaResults(xV,yV).p_group       = p_grp;
        ancovaResults(xV,yV).df            = df;
        ancovaResults(xV,yV).t_interaction = t_int;
        ancovaResults(xV,yV).p_interaction = p_int;
        ancovaResults(xV,yV).n             = height(tbl);

        fprintf(['%-12s vs %-16s | ANCOVA group effect: t(%d) = %6.3f, p = %.4g', ...
                 '  ||  slope-heterogeneity (X*Group): t = %6.3f, p = %.4g\n'], ...
                 yVars{yV}, xVars{xV}, df, t_grp, p_grp, t_int, p_int);
    end
end

figure('Position', [100, 100, 2000, 300]);



plotIdx = 1;
for xV = 1:3
    for yV = 1:2

        ax = subplot(3, 2, plotIdx);
        hold on;

        % --- Build shared bin edges and centers ---
        % Pool valid x data across both F for shared edges
        allX = [];
        for F = 1:2
            xData = eval([xVars{xV} '{' num2str(F) '}']);
            yData = eval([yVars{yV} '{' num2str(F) '}']);
            valid = ~isnan(xData) & ~isnan(yData);
            allX  = [allX; xData(valid)];
        end


                xEdges = prctile(allX, linspace(0, 100, nBins + 1));
                xEdges = unique(xEdges);   % remove duplicates if data has ties
                xCenters = xEdges(1:end-1);


        nBinsAct = numel(xEdges) - 1;   % single source of truth for both methods
        % ------------------------------------------

        for F = 1:2
            xData = eval([xVars{xV} '{' num2str(F) '}']);
            yData = eval([yVars{yV} '{' num2str(F) '}']);

            valid = ~isnan(xData) & ~isnan(yData);
            xData = xData(valid);
            yData = yData(valid);

            yMean = NaN(1, nBinsAct);
            ySEM  = NaN(1, nBinsAct);
            nPerBin = zeros(1, nBinsAct);

            for b = 1:nBinsAct
                inBin = xData > xEdges(b) & xData <= xEdges(b+1);
                vals  = yData(inBin);
                nPerBin(b) = numel(vals);
                if numel(vals) >= 2
                    yMean(b) = nanmean(vals);
                    ySEM(b)  = nanstd(vals) / sqrt(numel(vals));
                elseif numel(vals) == 1
                    yMean(b) = vals;
                    ySEM(b)  = NaN;
                end
            end

            validBins = ~isnan(yMean);
            xPlot = double(xCenters(validBins));
            yPlot = double(yMean(validBins));
            ePlot = ySEM(validBins);
            nPlot = nPerBin(validBins);

            % Shaded error region
            fill([xPlot, fliplr(xPlot)], ...
                 [yPlot + ePlot, fliplr(yPlot - ePlot)], ...
                 Params.colors{F}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

            % Main line
            plot(xPlot, yPlot, ...
                 'Color', Params.colors{F}, ...
                 'LineWidth', 1);

            % Optionally write N per bin
            if showN
                for b = 1:numel(xPlot)
                    text(xPlot(b), yPlot(b), [num2str(nPlot(b))], ...
                         'Color', Params.colors{F}, ...
                         'FontSize', 7, ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'bottom');
                end
            end
        end

        xlabel(xLabels{xV}, 'FontSize', 11);
        ylabel(yLabels{yV}, 'FontSize', 11);

        box off;
        hold off;

        plotIdx = plotIdx + 1;
        
        set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.03, 0.03]);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
    end
end

set(gcf,'Units','centimeters','Position',[10 10 5.5 6.8]);

PlotTitle='_Fit controls Resp, range and dprime ';
basePath = fullfile(saveDir_FigS5,[PlotTitle]);
saveFigAll(gcf, basePath);