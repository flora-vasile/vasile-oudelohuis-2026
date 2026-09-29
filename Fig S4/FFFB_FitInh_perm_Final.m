saveDir_FigS4 = fullfile(Params.MainDir, 'plot', 'Fig S4');
if ~exist(saveDir_FigS4,'dir')
    mkdir(saveDir_FigS4)
end

clear res valid

minRespRange = 0.04;

%% Linear regression

for F=1:2      
    
    for ani=1:Params.n_ani{F}
        for pos=1:Params.n_pos{F}(ani)
            if ~isempty(data{F}{ani,pos})
                visTrials{F}{ani,pos}=data{F}{ani,pos}.Raw(:,:,2:9,:,:,:);
                photoTrials{F}{ani,pos}=squeeze(data{F}{ani,pos}.Raw(:,:,1,1,:,:));
                
                visTrialsCat{F}{ani,pos}=reshape(visTrials{F}{ani,pos},size(visTrials{F}{ani,pos},1),size(visTrials{F}{ani,pos},2),size(visTrials{F}{ani,pos},3)*size(visTrials{F}{ani,pos},4),size(visTrials{F}{ani,pos},5),size(visTrials{F}{ani,pos},6));
                allTrialsCat{F}{ani,pos}(:,:,1,:,:)=photoTrials{F}{ani,pos};
                allTrialsCat{F}{ani,pos}(:,:,2:17,:,:)=visTrialsCat{F}{ani,pos};
                
                condResp{F}{ani,pos}=squeeze(nanmean(allTrialsCat{F}{ani,pos}(:,Params.StimONResponse,:,:,:),2));
                condMean{F}{ani,pos}=squeeze(nanmean(condResp{F}{ani,pos},4));
                               
               condRespOn{F}{ani,pos}=NaN(size(condResp{F}{ani,pos}));
                condRespOn{F}{ani,pos}(:,1,:,:)=condResp{F}{ani,pos}(:,1,:,:);
                for ROI=1:size(condRespOn{F}{ani,pos},1)
                    
                    for Dir=2:size(condRespOn{F}{ani,pos},2)
                        if dirResponsive{F}{ani,pos}(ROI,Dir-1)
                            condRespOn{F}{ani,pos}(ROI,Dir,:,:)=condResp{F}{ani,pos}(ROI,Dir,:,:);
                        end
                    end
                end
        end
    end
end
    
respAll{F}=cat(1,condRespOn{F}{:});

respInh{F}=respAll{F}(inhibited{F},:,:,:);

nInhibited = size(respInh{F},1);
stimVals   = (1:size(respInh{F},2))';   % or actual stimulus values

% Rebuild fresh each run so leftover data from a previous run (with a
% different number of neurons) can't leave valid{F}/res{F} mismatched.
valid{F} = false(nInhibited,1);

for n = 1:nInhibited
    respN = squeeze(respInh{F}(n,:,:,:));  

    
    offMean = squeeze(nanmean(respN(:,1,:),3));
    if sum(~isnan(offMean)) < 3 || range(offMean) < minRespRange || max(offMean) < minRespRange
        continue
    end 

    valid{F}(n) = true;
        respN(any(any(isnan(respN),3),2),:,:)=[];

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
ylabel('Fraction of inhibited neurons');
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

PlotTitle='_Gain categorisation - perm Inhibited';
basePath = fullfile(saveDir_FigS4,[PlotTitle]);
saveFigAll(gcf, basePath);