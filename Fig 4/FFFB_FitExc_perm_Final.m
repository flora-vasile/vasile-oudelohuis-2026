saveDir_Fig4 = fullfile(Params.MainDir, 'plot', 'Fig 4');
if ~exist(saveDir_Fig4,'dir')
    mkdir(saveDir_Fig4)
end

minRespRange = 0.04;

%% Linear regression

for F=1:2
    
    for ani=1:Params.n_ani{F}
        for pos=1:Params.n_pos{F}(ani)
            if isempty(data{F}{ani,pos})
                continue
            end
            
            raw = data{F}{ani,pos}.Raw;              % [cell x time x stim x dir x laser x rep]
            
            % Visual-only stimuli (conditions 2-9) and photostim-only condition (1,1)
            visTrials{F}{ani,pos}   = raw(:,:,2:9,:,:,:);
            photoTrials{F}{ani,pos} = squeeze(raw(:,:,1,1,:,:));
            
            % Collapse stim-type x direction into a single "condition" dimension
            [nCells, nTime, nStim, nDir, nLaser, nRep] = size(visTrials{F}{ani,pos});
            visTrialsCat{F}{ani,pos} = reshape(visTrials{F}{ani,pos}, ...
                nCells, nTime, nStim*nDir, nLaser, nRep);
            
            % Prepend the photostim-only condition as condition 1.
            photoTrialsPad = reshape(photoTrials{F}{ani,pos}, nCells, nTime, 1, nLaser, nRep);
            allTrialsCat{F}{ani,pos} = cat(3, photoTrialsPad, visTrialsCat{F}{ani,pos});
            
            % Average over the stimulus-on window, then over repeats
            condResp{F}{ani,pos}=squeeze(nanmean(allTrialsCat{F}{ani,pos}(:,Params.StimONResponse,:,:,:),2));
            
            stimResp{F}{ani,pos}=squeeze(nanmean(visTrials{F}{ani,pos}(:,Params.StimONResponse,:,:,:,:),2));
            
            excitedCell{F}{ani,pos}=logical(CellStat{F}{ani,pos}.SigExcAD);
        end
    end
    
    %% Pool sessions/animals together
    
    excited{F} = logical(cat(1,excitedCell{F}{:}));
    respAll{F}=cat(1,condResp{F}{:});
    
    respExc{F}=respAll{F}(excited{F},:,:,:);
    
    %% Linear regression Excited cells
    
    nExcited = size(respExc{F},1);
    stimVals   = (1:size(respExc{F},2))';   % or actual stimulus values
    
    for n = 1:nExcited
        valid{F}(n) = false;
        
        respN = squeeze(respExc{F}(n,:,:,:));
        
        
        offMean = squeeze(mean(respN(:,1,:),3));
        if sum(~isnan(offMean)) < 3 || range(offMean) < minRespRange || max(offMean) < minRespRange
            continue
        end
        
        valid{F}(n) = true;
        %_labelShuffle
        res{F}(n) = PermTest( ...
            respN, ...
            'nPerm', 5000, ...
            'alpha', 0.05 );
    end
    
end


%% Bar plot
%% Define categories
categories = {'Mult','Div', 'Add','Sub'};
nCat = numel(categories);
nConditions = 2; % FB, FF

counts = zeros(nConditions, nCat);

for F = 1:2
    
    nNeurons = numel(res{F});
    
    for n = 1:nNeurons
        
        % Only consider valid neurons
        if ~valid{F}(n)
            continue
        end
        
        % Extract significance
        isMult = res{F}(n).sig_mult;
        isAdd  = res{F}(n).sig_add;
        isSub = res{F}(n).sig_sub;
        isDiv = res{F}(n).sig_div;
        
        %         Categorise
        %Option 1
        if isMult
            counts(F,1) = counts(F,1) + 1; % Multiplicative
        end
        if isDiv
            counts(F,2) = counts(F,2) + 1; % Divisive
        end
        if isAdd
            counts(F,3) = counts(F,3) + 1; % Additive
        end
        if isSub
            counts(F,4) = counts(F,4) + 1; % Subtractive
        end
    end
    
    totals(F) = sum(valid{F});
end

%% Fractions
fractions = counts ./ totals(:);

%% Plot grouped bar chart
% plot FF first and FB second

figure; hold on;
b = bar([fractions(2,:);fractions(1,:)]', 'grouped');

% b = bar(fractions', 'grouped');

% Colors (use whatever you already use for FB / FF)
b(1).FaceColor = Params.colors{2};  % FB
b(2).FaceColor = Params.colors{1};  % FF

xticks(1:nCat);
xticklabels(categories);
ylabel('Fraction of excited neurons');
ylim([0 1]);
box off;

legend(Params.ExpName{2},Params.ExpName{1}, ...
    'Location','northwest'); legend boxoff


%% Fisher exact test per category
fprintf('\nFisher''s exact test (FB vs FF):\n');
fprintf('-------------------------------------------\n');

for c = 1:nCat
    tab = [ counts(1,c), totals(1)-counts(1,c);
        counts(2,c), totals(2)-counts(2,c) ];
    
    try
        [~, p, stats] = fishertest(tab);
        
        orText = sprintf(' OR=%.2f', stats.OddsRatio);
        
    catch
        p = NaN;
        orText = '';
    end
    
    %     fprintf('%s: FB=%d / FF=%d    p=%.4g%s\n', ...
    %         categories{c}, counts(1,c), counts(2,c), p);
    
    fprintf('%s: FB=%d / FF=%d    p=%.4g%s\n', ...
        categories{c}, counts(1,c), counts(2,c), p, orText);
    
    
    xmid = mean([b(1).XData(c)+b(1).XOffset, ...
        b(2).XData(c)+b(2).XOffset]);
    ymax = max(fractions(:,c));
    text(xmid, ymax + 0.08, num2str(p), ...
        'HorizontalAlignment','center', ...
        'FontSize',4);
end

ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
%despline(gca);
set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
set(gcf,'Units','centimeters','Position',[10 10 5 3.2]);

PlotTitle='_Gain categorisation - perm Excited';
basePath = fullfile(saveDir_Fig4,[PlotTitle]);
saveFigAll(gcf, basePath);