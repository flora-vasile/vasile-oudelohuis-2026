saveDir_Fig5 = [MainDir filesep 'plot\Fig 5'];
if ~exist(saveDir_Fig5,'dir')
    mkdir(saveDir_Fig5)
end

%% ================= PARAMETERS =================
CellSubSelecName = 'SFTF pooled responsive';

nNeuronsList = [2 5 10 20 50 100];  % neuron-count bins
nRepeats = 1;  %20                  % subsampling repetitions

%% ================= MAIN LOOP =================
for F = 1:2
    % Initialize output containers
    for nIdx = 1:numel(nNeuronsList)
        OFF_perf{F}{nIdx} = cell(n_ani{F}, max(n_pos{F}));
    end
    
    for ani = 1:n_ani{F}
        for pos = 1:n_pos{F}(ani)
            
            if isempty(data{F}{ani,pos})
                continue
            end
            
            % Cell selection
            cellMask = logical(CellStat{F}{ani,pos}.responsiveAD);
            nNeuronsTotal = sum(cellMask);
            
            if nNeuronsTotal < min(nNeuronsList)
                continue
            end
            
            disp(['F=' num2str(F) ' Mouse=' num2str(ani) ' Pos=' num2str(pos) ...
                ' nCells=' num2str(nNeuronsTotal)])
            
            %% ----------- DATA PREPARATION -----------
            Vis = squeeze(nanmean(data{F}{ani,pos}.Raw(cellMask, Params.StimONResponse, 2:9, :, :, :),2));
            
            % concatenate trials from 2 SFTF
            Vis2 = permute(Vis, [1 2 4 5 3]);
            MeanTime = reshape(Vis2, size(Vis2,1), size(Vis2,2), size(Vis2,3), []);
            
            % Dimensions
            [nNeurons, nDir, nLaser, nReps] = size(MeanTime);
            
            % Split by laser
            data_off = squeeze(MeanTime(:,:,1,:)); % neurons × dir × reps
            
            % Reshape
            data_off = reshape(permute(data_off, [1 3 2]), nNeurons, []);
            
            dirIdentity = repelem(1:nDir, nReps);
            
            %% ----------- DECODER FOR EACH NEURON COUNT -----------
            for nIdx = 1:numel(nNeuronsList)
                
                nNeuronsSub = nNeuronsList(nIdx);
                
                % Skip if not enough neurons
                if nNeurons < nNeuronsSub
                    OFF_perf{F}{nIdx}{ani,pos} = NaN;
                    continue
                end
                
                OFF_rep = nan(nRepeats,1);
                
                for rep = 1:nRepeats
                    
                    idxNeurons = randperm(nNeurons, nNeuronsSub);
                    
                    % ----- Laser OFF -----
                    X = data_off(idxNeurons,:)';
                    y = dirIdentity(:);
                    
                    t = templateSVM('KernelFunction','linear','Standardize',true);
                    Mdl = fitcecoc(X, y, 'Learners', t, 'KFold', 5);
                    OFF_rep(rep) = (1 - kfoldLoss(Mdl)) * 100;
                end
                
                % Store mean across repetitions
                OFF_perf{F}{nIdx}{ani,pos} = mean(OFF_rep);
            end % neuron count
        end % pos
    end % ani
end % F


%% Plot
sem = @(x) nanstd(x) ./ sqrt(sum(~isnan(x)));

chanceLevel = 12.5; % adjust if different

nBins = numel(nNeuronsList);

mean_OFF = nan(2, nBins);
sem_OFF  = nan(2, nBins);

for F = 1:2
    for nIdx = 1:nBins
        
        vals_OFF = [];
        
        for ani = 1:n_ani{F}
            for pos = 1:n_pos{F}(ani)
                
                if ~isempty(OFF_perf{F}{nIdx}{ani,pos})
                    vals_OFF(end+1) = OFF_perf{F}{nIdx}{ani,pos}; %#ok<SAGROW>
                end
                
            end
        end
        
        mean_OFF(F,nIdx) = nanmean(vals_OFF);
        sem_OFF(F,nIdx)  = sem(vals_OFF);
    end
end

%% ---- Count total number of sessions used ----

nSessionsUsed = zeros(1,2);

for F = 1:2
    nAnimalsUsed = 0;

    for ani = 1:n_ani{F}
        animalUsed = false;

        for pos = 1:n_pos{F}(ani)

            sessionUsed = false;

            for nIdx = 1:nBins
                val = OFF_perf{F}{nIdx}{ani,pos};
                if ~isempty(val) && ~isnan(val)
                    sessionUsed = true;
                    break
                end
            end

            if sessionUsed
                nSessionsUsed(F) = nSessionsUsed(F) + 1;
                animalUsed = true;
            end

        end

        if animalUsed
            nAnimalsUsed = nAnimalsUsed + 1;
        end
    end

    fprintf('F=%d: %d sessions used in the plot, %d animals used\n', F, nSessionsUsed(F), nAnimalsUsed);
end



%% Plot

figure; hold on;

groupOrder = [2 1]; % FF first, then FB
groupNames = {'PM','V1'};
colors = [Params.colors{2};   % FF
    Params.colors{1}];  % FB

for sp = 1:2
    
    F = groupOrder(sp);
    col = colors(sp,:);
    
    
    %% ---- OFF mean + SEM (mean across all animals & positions) ----
    x = nNeuronsList;
    y = mean_OFF(F,:);  % already computed
    e = sem_OFF(F,:);   % already computed
    
    fill([x fliplr(x)], [y-e fliplr(y+e)], col, ...
        'FaceAlpha',0.25,'EdgeColor','none');
    
    plot(x, y, '-', 'Color', col, 'LineWidth',1,'Marker','o','MarkerSize',1.8,'MarkerFaceColor',col);
    text(20,80,groupNames{sp},'Fontweight','normal', 'Color', col, 'FontSize',5);
end

%% ---- Formatting ----
plot(xlim, [chanceLevel chanceLevel], '--','Color',[0.5 0.5 0.5], 'LineWidth',1);
text(80,20,'Chance', 'FontSize',5);
xlabel('Population size');
ylabel('Decoding performance (%)');
xlim([1.5 110]); ylim([0 100]);
box off;
offsetAxes(gca);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
set(gca,'FontSize',5,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);

set(gcf,'Units','centimeters','Position',[10 10 3.5 2.75]);

PlotTitle='_Decoder performance vs population size - Sessions';
basePath = fullfile(saveDir_Fig5,[tv PlotTitle]);
saveFigAll(gcf, basePath);
