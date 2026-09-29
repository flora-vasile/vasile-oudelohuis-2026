

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

saveDir_FigS7 = [Params.MainDir filesep 'plot\Fig S7'];
if ~exist(saveDir_FigS7,'dir')
    mkdir(saveDir_FigS7)
end



%% PARAMETERS
DistVec = 40:10:100;   % µm
nDist   = numel(DistVec);
ThNbneigh=2;
Niter    = 10000;      % bootstrap iterations (neighbours)
clear Results

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
LaserAlone{F}{ani,pos}=squeeze(nanmean(nanmean(data{F}{ani,pos}.Raw(:,Params.StimONResponse,1,1,2,:),2),4));
             
             
            coords = CellStat{F}{ani,pos}.Coordyx_um;
             isLoop = CellStat{F}{ani,pos}.Loop;
            isExc  = CellStat{F}{ani,pos}.SigExcLaserAlone;
            isInh  = CellStat{F}{ani,pos}.SigInhLaserAlone;
            isMod  = isExc | isInh;
            resp   = LaserAlone{F}{ani,pos};
            
            loopIdx = find(isLoop);
            nonloopIdx = find(~isLoop);
            
            % store all loop cells (no filtering yet)
            % store all loop cells (no filtering yet)
            for k = loopIdx(:)'
                allLoop(end+1).coords = coords(k,:);
                allLoop(end).session = [ani,pos];
                allLoop(end).isExc = isExc(k);
                allLoop(end).isInh = isInh(k);
                allLoop(end).isMod = isMod(k);
                allLoop(end).resp = resp(k);
            end
            
            % store all non-loop candidates
            for k = nonloopIdx(:)'
                allNeighCandidates(end+1).coords = coords(k,:);
                allNeighCandidates(end).session = [ani,pos];
                allNeighCandidates(end).isExc = isExc(k);
                allNeighCandidates(end).isInh = isInh(k);
                allNeighCandidates(end).isMod = isMod(k);
                allNeighCandidates(end).resp = resp(k);
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
    isExcAll   = false(Ncells,1);
    isInhAll   = false(Ncells,1);
    isModAll   = false(Ncells,1);
    respAll      = zeros(Ncells,1);
    
    % loops
    for i = 1:nLoopTotal
        coordsAll(i,:)  = allLoop(i).coords;
        sessionAll(i,:) = allLoop(i).session;
        isExcAll(i)     = allLoop(i).isExc;
        isInhAll(i)     = allLoop(i).isInh;
        isModAll(i)     = allLoop(i).isMod;
        respAll(i)        = allLoop(i).resp;
    end
    
    % neighbors
    for j = 1:nNeighTotal
        idx = nLoopTotal + j;
        coordsAll(idx,:)  = allNeighCandidates(j).coords;
        sessionAll(idx,:) = allNeighCandidates(j).session;
        isExcAll(idx)     = allNeighCandidates(j).isExc;
        isInhAll(idx)     = allNeighCandidates(j).isInh;
        isModAll(idx)     = allNeighCandidates(j).isMod;
        respAll(idx)        = allNeighCandidates(j).resp;
    end


    %% -------------------------
% 3) PRECOMPUTE NEIGHBORS PER LOOP
%% -------------------------
for d = 1:nDist
    
    DistThr = DistVec(d);
    fprintf('F=%d | Distance threshold = %d µm\n',F,DistThr);

    %% -------------------------
    % 3) PRECOMPUTE NEIGHBORS PER LOOP (IDENTICAL LOGIC)
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

Results(F).nLoop(d) = nLoop;

    %% -------------------------
    % 4) BOOTSTRAP 
    %% -------------------------

    realExc = mean(isExcAll(validIdx));
    realInh = mean(isInhAll(validIdx));
    realMod = mean(isModAll(validIdx));
    realresp  = mean(respAll(validIdx));

    bootExc = nan(Niter,1);
    bootInh = nan(Niter,1);
    bootMod = nan(Niter,1);
    bootresp  = nan(Niter,1);

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
        bootresp(it)  = mean(respAll(sampledIdx));
    end

    metrics = {'Exc','Inh','Mod','PhotoResp'};
    realVals = [realExc realInh realMod realresp];
    bootVals = {bootExc,bootInh,bootMod,bootresp};

    for m = 1:4
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

end 
  

%% Plot for figure
metricNames = {'Exc'};
nlpname{1}= 'V1_{ND}';
lpname{1}= 'V1_{PM}';
nlpname{2}= 'PM_{ND}';
lpname{2}= 'PM_{V1}';
figure;

for F = 1:2
    for m = 1
        subplot(2,1,(3-F-1)*1 + m); hold on;

        realVals = Results(F).Metric(m).Real;
        CI = Results(F).Metric(m).CI;
        pvals = Results(F).Metric(m).p;

        % CI shading (neighbors)
        f1=fill([DistVec fliplr(DistVec)], ...
             [CI(1,:) fliplr(CI(2,:))], ...
             [0.5 0.5 0.5], ...
             'FaceAlpha',0.25, ...
             'EdgeColor','none');

         % looped cells
        p1=plot(DistVec, realVals, '-o','MarkerSize',2,...
            'Color','r', ...
            'MarkerFaceColor', 'r', ...
            'LineWidth', 1);
legend({[nlpname{F} ' CI'] lpname{F} }); legend boxoff

for d = 1:numel(DistVec)
    % significance marker
    if pvals(d) < 0.05
        text(DistVec(d), realVals(d), '*', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'FontSize',10);
    end
end

        title(sprintf('%s – %s', Params.ExpName{F}, metricNames{m}));
        if F==1
        xlabel('Distance threshold (µm)');
            ylabel('Fraction');
        end
        box off
        set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.02, 0.01]);
% offsetAxes(gca);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';  
    end
end

 set(gcf,'Units','centimeters','Position',[10 10 1.5 6]);

 PlotTitle=['Looped vs bootstrapped neighbours (distances)' num2str(Niter) ' bootstraps Dist Threshold - Laser only Exc '];
basePath = fullfile(saveDir_FigS7,[PlotTitle]);
saveFigAll(gcf, basePath);
