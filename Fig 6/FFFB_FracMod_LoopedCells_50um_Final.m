

% This script compares functional properties of looped cells to those of
% nearby non-looped neighbours within an increasingly large distance radius.
%
% Analysis steps:
% • Collect all responsive looped and non-looped cells across animals and sessions
% • Store cell coordinates, session identity, response class (Exc/Inh/Mod),
%   direction index (DI), and absolute DI
% • For each looped cell, identify non-looped neighbours from the same session
%   that fall within a cumulative distance threshold (D ? DistThr)
% • Remove looped cells that have fewer than two eligible neighbours
% • For the remaining looped cells:
%     – Compute the real fraction (or mean) of each metric in looped cells
%     – Generate a null distribution by bootstrap resampling:
%         · For each looped cell, randomly select one neighbour from its local pool
%         · Repeat across all looped cells and iterations
% • Compute bootstrap mean, 95% confidence intervals, and two-sided p-values
%   comparing looped-cell values to the neighbour null distribution
% • Repeat the entire procedure for a range of increasing distance thresholds
%   (e.g. 10 ? 200 µm), such that each threshold includes all closer neighbours
% • Plot looped-cell values as a function of distance threshold, overlaid with
%   neighbour confidence intervals, p-values, and the number of contributing
%   looped cells at each distance

saveDir_Fig6 = [Params.MainDir filesep 'plot\Fig 6'];
if ~exist(saveDir_Fig6,'dir')
    mkdir(saveDir_Fig6)
end

%% PARAMETERS
DistVec = 50;   % µm
nDist   = numel(DistVec);
ThNbneigh=2;
Niter    = 10000;      % bootstrap iterations (neighbours)

for F = 1:2
    for d = 1:nDist
        Results(F).Dist(d) = DistVec(d);
    end
end

for F = 1:2
    fprintf('Running F=%d\n',F);
    
    %% -------------------------
    % 1) BUILD LOOP + NEIGHBOR LISTS
    %% -------------------------
    allLoop = [];             % all loop cells
    allNeighCandidates = [];  % all non-loop cells
    loopSessionIdx = [];      % to quickly find same-session pseudo-neighbors
    
    for ani = 1:Params.n_ani{F}
        for pos = 1:Params.n_pos{F}(ani)
             if isempty(data{F}{ani,pos}), continue; end
           AllDirRespAllRepSes{F}{ani,pos}=squeeze(mean(data{F}{ani,pos}.StimAv(:,Params.StimONResponse,:,:),2));
                AllDirRespSes{F}{ani,pos}=squeeze(mean(mean(data{F}{ani,pos}.StimAv(:,Params.StimONResponse,:,:),2),4));
                
                dprime{F}{ani,pos}=(AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,2)...
                    -AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,1))...
                    ./sqrt(1/2*(var(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,2,:),0,3)...
                    +var(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,1,:),0,3)));
                
            Absdprime{F}{ani,pos}=abs(dprime{F}{ani,pos});
            VisResp{F}{ani,pos}=AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,1);            
            
            coords = CellStat{F}{ani,pos}.Coordyx_um(CellStat{F}{ani,pos}.responsiveAD,:);
             isLoop = CellStat{F}{ani,pos}.Loop(CellStat{F}{ani,pos}.responsiveAD);
            isExc  = CellStat{F}{ani,pos}.SigExcAD(CellStat{F}{ani,pos}.responsiveAD);
            isInh  = CellStat{F}{ani,pos}.SigInhAD(CellStat{F}{ani,pos}.responsiveAD);
            isMod  = isExc | isInh;
            DI = dprime{F}{ani,pos};
            AbsDI = Absdprime{F}{ani,pos};
            
            loopIdx = find(isLoop);
            nonloopIdx = find(~isLoop);
            
            origROIidx = find(CellStat{F}{ani,pos}.responsiveAD);

            
            % store all loop cells (no filtering yet)
            for k = loopIdx(:)'
                allLoop(end+1).coords = coords(k,:);
                allLoop(end).session = [ani,pos];
                allLoop(end).roi        = origROIidx(k);
                allLoop(end).isExc = isExc(k);
                allLoop(end).isInh = isInh(k);
                allLoop(end).isMod = isMod(k);
                allLoop(end).DI = DI(k);
                allLoop(end).AbsDI = AbsDI(k);
                allLoop(end).VisResp = VisResp{F}{ani,pos}(k);
                
            end
            
            % store all non-loop candidates
            for k = nonloopIdx(:)'
                allNeighCandidates(end+1).coords = coords(k,:);
                allNeighCandidates(end).session = [ani,pos];
                allNeighCandidates(end).roi      = origROIidx(k);
                allNeighCandidates(end).isExc = isExc(k);
                allNeighCandidates(end).isInh = isInh(k);
                allNeighCandidates(end).isMod = isMod(k);
                allNeighCandidates(end).DI = DI(k);
                allNeighCandidates(end).AbsDI = AbsDI(k);
                allNeighCandidates(end).VisResp = VisResp{F}{ani,pos}(k);

            end
        end
    end
    
    nLoopTotal  = numel(allLoop);
    nNeighTotal = numel(allNeighCandidates);
    Ncells = nLoopTotal + nNeighTotal;
    
    %% -------------------------
    % 2) BUILD COMPACT ARRAYS
    %% -------------------------
    coordsAll  = zeros(Ncells,2);
    sessionAll = zeros(Ncells,2);
    roiAll     = nan(Ncells,1);   
    isExcAll   = false(Ncells,1);
    isInhAll   = false(Ncells,1);
    isModAll   = false(Ncells,1);
    DIAll      = zeros(Ncells,1);
    AbsDIAll   = zeros(Ncells,1);
    VisRespAll = zeros(Ncells,1);
    
    % loops
    for i = 1:nLoopTotal
        coordsAll(i,:)  = allLoop(i).coords;
        sessionAll(i,:) = allLoop(i).session;
        roiAll(i)       = allLoop(i).roi;
        isExcAll(i)     = allLoop(i).isExc;
        isInhAll(i)     = allLoop(i).isInh;
        isModAll(i)     = allLoop(i).isMod;
        DIAll(i)        = allLoop(i).DI;
        AbsDIAll(i)     = allLoop(i).AbsDI;
        VisRespAll(i) = allLoop(i).VisResp;

    end
    
    % neighbors
    for j = 1:nNeighTotal
        idx = nLoopTotal + j;
        coordsAll(idx,:)  = allNeighCandidates(j).coords;
        sessionAll(idx,:) = allNeighCandidates(j).session;
        roiAll(idx)       = allNeighCandidates(j).roi;
        isExcAll(idx)     = allNeighCandidates(j).isExc;
        isInhAll(idx)     = allNeighCandidates(j).isInh;
        isModAll(idx)     = allNeighCandidates(j).isMod;
        DIAll(idx)        = allNeighCandidates(j).DI;
        AbsDIAll(idx)     = allNeighCandidates(j).AbsDI;
        VisRespAll(idx) = allNeighCandidates(j).VisResp;

    end

    
    %% -------------------------
% 3) PRECOMPUTE NEIGHBORS PER LOOP
%% -------------------------
for d = 1:nDist
    
    DistThr = DistVec(d);
    fprintf('F=%d | Distance threshold = %d µm\n',F,DistThr);

    %% -------------------------
    % 3) PRECOMPUTE NEIGHBORS PER LOOP 
    %% -------------------------

    groupPerLoop = cell(nLoopTotal,1);
    neighborsPerLoop = cell(nLoopTotal,1);
    invalidLoop = false(nLoopTotal,1);
    neighGlobalIdx = (nLoopTotal+1):(nLoopTotal+nNeighTotal);

    for i = 1:nLoopTotal
        loopCoord = coordsAll(i,:);
        loopSess  = sessionAll(i,:);

        sameSessMask = (sessionAll(neighGlobalIdx,1)==loopSess(1)) & ...
                       (sessionAll(neighGlobalIdx,2)==loopSess(2));
        sameNeighIdx = neighGlobalIdx(sameSessMask);

        neighList = [];
        if ~isempty(sameNeighIdx)
            neighCoords = coordsAll(sameNeighIdx,:);
            D = sqrt(sum((neighCoords - loopCoord).^2,2));
            neighList = sameNeighIdx(D <= DistThr);
        end

        if isempty(neighList) || numel(neighList)<ThNbneigh
            invalidLoop(i) = true;
            continue
        end

        neighborsPerLoop{i} = neighList(:);
        groupPerLoop{i} = unique([i; neighList(:)]);
    end

% --- REMOVE INVALID LOOPS (UNCHANGED LOGIC) ---
validIdx = find(~invalidLoop);
nLoop = numel(validIdx);

% --- Unique neighbours actually used across the valid loop cells ---
allNeighUsed = unique(vertcat(neighborsPerLoop{validIdx}));
nNeighUsed = numel(allNeighUsed);

animalsUsed = unique(sessionAll(validIdx,1));
nAnimalsUsed = numel(animalsUsed);

Results(F).nLoop(d) = nLoop;
Results(F).nNeighUsed(d) = nNeighUsed;

fprintf('F=%d | DistThr=%d um | nLoopCells=%d | nUniqueNeighbours=%d | nAnimals=%d\n', ...
    F, DistThr, nLoop, nNeighUsed, nAnimalsUsed);


    %% -------------------------
    % 4) BOOTSTRAP 
    %% -------------------------

    realExc = mean(isExcAll(validIdx));
    realInh = mean(isInhAll(validIdx));
    realMod = mean(isModAll(validIdx));
    realDI  = mean(DIAll(validIdx));
    realAbsDI = mean(AbsDIAll(validIdx));
    realVisResp = mean(VisRespAll(validIdx));

    bootExc = nan(Niter,1);
    bootInh = nan(Niter,1);
    bootMod = nan(Niter,1);
    bootDI  = nan(Niter,1);
    bootAbsDI = nan(Niter,1);
    bootVisResp = nan(Niter,1);


    for it = 1:Niter
        sampledIdx = nan(nLoop,1);
        for k = 1:nLoop
            g = groupPerLoop{validIdx(k)};
            g(g==validIdx(k)) = [];
            sampledIdx(k) = g(randi(numel(g)));
        end

        bootExc(it) = mean(isExcAll(sampledIdx));
        bootInh(it) = mean(isInhAll(sampledIdx));
        bootMod(it) = mean(isModAll(sampledIdx));
        bootDI(it)  = mean(DIAll(sampledIdx));
        bootAbsDI(it) = mean(AbsDIAll(sampledIdx));
        bootVisResp(it) = mean(VisRespAll(sampledIdx));

    end

metrics  = {'Exc','Inh','Mod','DI','AbsDI','VisResp'};
realVals = [realExc realInh realMod realDI realAbsDI realVisResp];
bootVals = {bootExc,bootInh,bootMod,bootDI,bootAbsDI, bootVisResp};


    for m = 1:6
        dist = bootVals{m};
        dist = dist(~isnan(dist));

        Results(F).Metric(m).Real(d) = realVals(m);
        Results(F).Metric(m).CI(:,d) = prctile(dist,[2.5 97.5]);
        Results(F).Metric(m).Mean(d) = mean(dist);
        Results(F).Metric(m).dist(:,d) = dist;
        Results(F).Metric(m).p(d) = ...
            2*min(mean(dist>realVals(m)),mean(dist<realVals(m)));
    end

end

%nKept{F}=size(keepLoopIdx,1);

end 
  

%% Plot 
metricNames = {'Pos.','Neg.','Mod.'};
nlpname{1}= 'V1_{ND}';
lpname{1}= 'V1_{PM}';
nlpname{2}= 'PM_{ND}';
lpname{2}= 'PM_{V1}';

figure;


for F = 1:2
    
    for m = 1:length(metricNames)
        subplot(2,length(metricNames),(3-F-1)*length(metricNames) + m); hold on;

        realVals = Results(F).Metric(m).Real;
        CI = Results(F).Metric(m).CI;
        dist = Results(F).Metric(m).dist;
        pvals = Results(F).Metric(m).p;

        
        
          % histogram / density of bootstrap distribution
        histogram(dist,40,'Normalization','probability', ...
            'FaceAlpha',0.6,'EdgeColor','none','FaceColor',[0.4 0.4 0.4],'binwidth',0.01);

        % plot vertical lines
        yl = ylim;
        plot([realVals realVals], yl, 'LineWidth', 1.5, 'Color', 'r');
        plot([CI(1) CI(1)], yl, ':', 'Color', [0.5 0.5 0.5]);
        plot([CI(2) CI(2)], yl, ':', 'Color', [0.4 0.5 0.5]);

if m==1
        text(0.05,(yl(1)+yl(2))./1.3 , nlpname{F}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',5, 'Color', [0.4 0.4 0.4]);
        text(realVals(d)+0.03,(yl(1)+yl(2))./1.3, lpname{F}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',5, 'Color', 'r');
        
        ylabel('Probability');
end
ylim(yl)
        
        
if F==2
        
            title(sprintf('%s %s', 'Frac.', metricNames{m}),'Fontweight','normal');
end
        
if pvals(d) < 0.001
        text(realVals(d)+0.03,(yl(1)+yl(2))./2 , '***', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',10,'Fontweight','bold', 'Color', 'r');
elseif pvals(d) < 0.01
        text(realVals(d)+0.03,(yl(1)+yl(2))./2 , '**', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',10,'Fontweight','bold', 'Color', 'r');    
    elseif pvals(d) < 0.05
        text(realVals(d)+0.03,(yl(1)+yl(2))./2 , '*', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',10,'Fontweight','bold', 'Color', 'r');
    end
        box off; 
set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.02, 0.01]);
% offsetAxes(gca);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';    
    
    
    end
end

 set(gcf,'Units','centimeters','Position',[10 10 6.7 4.3]);

PlotTitle=[' Looped vs bootstrapped neighbours ' num2str(DistVec) 'um ' num2str(Niter) ...
        ' bootstraps Dist Threshold - Th Nb neighbours ' num2str(ThNbneigh)];
basePath = fullfile(saveDir_Fig6,[PlotTitle]);
saveFigAll(gcf, basePath);
