saveDir_FigS3 = [Params.MainDir filesep 'plot\Fig S3'];
if ~exist(saveDir_FigS3,'dir')
    mkdir(saveDir_FigS3)
end

color_OFF{1}=[ .6 1 .8];
color_OFF{2}=[.8 .6 1];
rangeSizeResp=0.04;
nShuffle = 1; %100;

for F=1:2

        for ani=1:Params.n_ani{F}
        for pos=1:Params.n_pos{F}(ani)
            if ~isempty(data{F}{ani,pos})
    VisTrialsBestSFMeanTime{F}{ani,pos}=squeeze(nanmean(data{F}{ani,pos}.VisBestSF(:,Params.StimONResponse,:,:,:),2));
    VisTrialsBestSFMeanTimeTrials{F}{ani,pos}=squeeze(nanmean(VisTrialsBestSFMeanTime{F}{ani,pos},4));
    responsive{F}{ani,pos} = CellStat{F}{ani,pos}.responsiveADBestSF;

            end
        end
        end
    BestSFTFAll{F}=cat(1,VisTrialsBestSFMeanTimeTrials{F}{:});
    BestSFTFAllTrials{F}= cat(1,VisTrialsBestSFMeanTime{F}{:});
    resp= cat(1,responsive{F}{:});
    
    % Range of responses logical
    PassRangeTh=false(size(BestSFTFAll{F},1),1);
    for n=1:size(BestSFTFAll{F},1)
        if range(BestSFTFAll{F}(n,:,1)) > rangeSizeResp | range(BestSFTFAll{F}(n,:,2)) > rangeSizeResp
            PassRangeTh(n)=true;
        end
    end
    
    BestSFTFAllTrialsResp{F}=BestSFTFAllTrials{F}(resp&PassRangeTh,:,:,:);
    
    nDir = size(BestSFTFAllTrialsResp{F},2);
    dirDeg = linspace(0,315,nDir)';
    fwhm_from_kappa = @(k) 2 * acos(1 - log(2)/k);
    dPeak{F}    = [];   % delta amplitude at preferred direction (peakResp ON - OFF)
    dOpp{F}     = [];   % delta amplitude at opposite direction (oppResp ON - OFF)
    dMu{F}      = [];
    dFWHM{F}       = [];
    dB{F}       = [];
    
    dPeak_shuff{F}   = [];
    dOpp_shuff{F}    = [];
    dMu_shuff{F}     = [];
    dFWHM_shuff{F} = [];
    dB_shuff{F}      = [];
    
    
    for n=1:size(BestSFTFAllTrialsResp{F},1)
        disp(['F = ' num2str(F) ' | N = ' num2str(n)])
        theta = linspace(0,2*pi,9);
        theta(end) = [];
        
        % --- Trial-level data: [nDir x 2 x nTrials]
        respTrials = squeeze(BestSFTFAllTrialsResp{F}(n,:,:,:));
        
        % Skip neurons with no trial dimension
        if ndims(respTrials) < 3
            continue
        end
        
        nTrials = size(respTrials,3);
        
        % ===== REAL FIT =====
        respMean_OFF = squeeze(mean(respTrials(:,1,:),3));
        respMean_ON  = squeeze(mean(respTrials(:,2,:),3));
        
        fOff_real{F}(n) = fitVonMises_dir(theta, respMean_OFF);
        fOn_real{F}(n)  = fitVonMises_dir(theta, respMean_ON);
        
        dPeak{F}(end+1) = fOn_real{F}(n).peakResp - fOff_real{F}(n).peakResp;
        dOpp{F}(end+1)  = fOn_real{F}(n).oppResp  - fOff_real{F}(n).oppResp;
        dFWHM{F}(end+1) = rad2deg(fwhm_from_kappa(fOn_real{F}(n).kappa)) - ...
            rad2deg(fwhm_from_kappa(fOff_real{F}(n).kappa));
        dB{F}(end+1)    = fOn_real{F}(n).b        - fOff_real{F}(n).b;
        
        dmu_real = circ_dist(fOn_real{F}(n).mu, fOff_real{F}(n).mu);
        dMu{F}(end+1)   = rad2deg(dmu_real);
        
        % ===== SHUFFLE =====
        for s = 1:nShuffle
            
            % One permutation shared across all directions (preserves trial structure)
            shuffIdx = randperm(2*nTrials);
            
            respOFF = nan(nDir, nTrials);
            respON  = nan(nDir, nTrials);
            
            for d = 1:nDir
                trialsOFF  = squeeze(respTrials(d,1,:));   % [nTrials x 1]
                trialsON   = squeeze(respTrials(d,2,:));   % [nTrials x 1]
                allTrials  = [trialsOFF; trialsON];        % [2*nTrials x 1]
                allTrials  = allTrials(shuffIdx);          % same permutation for every direction
                respOFF(d,:) = allTrials(1:nTrials);
                respON(d,:)  = allTrials(nTrials+1:end);
            end
            
            respMean_OFF_sh = mean(respOFF, 2);
            respMean_ON_sh  = mean(respON,  2);
            
            try
                fOff_sh = fitVonMises_dir(theta, respMean_OFF_sh);
                fOn_sh  = fitVonMises_dir(theta, respMean_ON_sh);
            catch
                continue
            end
            
            dPeak_shuff{F}(end+1) = fOn_sh.peakResp - fOff_sh.peakResp;
            dOpp_shuff{F}(end+1)  = fOn_sh.oppResp  - fOff_sh.oppResp;
            dFWHM_shuff{F}(end+1) = rad2deg(fwhm_from_kappa(fOn_sh.kappa)) - ...
                rad2deg(fwhm_from_kappa(fOff_sh.kappa));
            dB_shuff{F}(end+1)    = fOn_sh.b        - fOff_sh.b;
            
            % Orientation-correct circular distance for shuffle
            dmu_sh                = circ_dist(2*fOn_sh.mu, 2*fOff_sh.mu) / 2;
            dMu_shuff{F}(end+1)   = rad2deg(dmu_sh);
        end
        
    end
end


%% Example: plot some neurons
for F = 1:2
    nNeurons = numel(fOff_real{F});
    exampleNeurons = randperm(nNeurons, min(10, nNeurons));
    
    figure; hold on
    
    for k = 1:numel(exampleNeurons)
        n = exampleNeurons(k);
        subplot(2,5,k)
        
        fitOff = fOff_real{F}(n);
        fitOn  = fOn_real{F}(n);
        
        plot(rad2deg(theta), fitOff.respMean, 'o','MarkerFaceColor',color_OFF{F},'MarkerEdgeColor','None','MarkerSize',3); hold on
        plot(rad2deg(theta), fitOn.respMean,  'o','MarkerFaceColor',Params.colors{F},'MarkerEdgeColor','None','MarkerSize',3);
        
        plot(rad2deg(fitOff.thetaFine), fitOff.curve,'Color', color_OFF{F},'LineWidth',1);
        plot(rad2deg(fitOn.thetaFine),  fitOn.curve, 'Color', Params.colors{F},'LineWidth',1);
        
        yl = ylim;
        plot([rad2deg(fitOff.mu) rad2deg(fitOff.mu)], yl, '--','Color', color_OFF{F})
        plot([rad2deg(fitOn.mu)  rad2deg(fitOn.mu)],  yl, '--', 'Color', Params.colors{F})
        
        % Convert kappa to FWHM in degrees
        fwhm_from_kappa = @(k) 2 * acos(1 - log(2)/k);
        fwhmOff = rad2deg(fwhm_from_kappa(fitOff.kappa));
        fwhmOn  = rad2deg(fwhm_from_kappa(fitOn.kappa));
        
        % Add kappa and FWHM as text
        text(0.05, 0.95, sprintf('OFF: \\kappa=%.2f (%.0f°)', fitOff.kappa, fwhmOff), ...
            'Units','normalized','FontSize',12, ...
            'VerticalAlignment','top');
        text(0.05, 0.82, sprintf('ON:  \\kappa=%.2f (%.0f°)', fitOn.kappa, fwhmOn), ...
            'Units','normalized','FontSize',12, ...
            'VerticalAlignment','top');
        
        xlabel('Direction (deg)');
        ylabel('Response');
        xticks(0:90:360);
        xlim([0 360]);
        title(sprintf('Neuron %d', n));
        box off
        pause(0.1)
        
        set(gca, 'Fontsize', 12,'linewidth',1)
        ax = gca;
        ax.XColor = 'k';
        ax.YColor = 'k';
        set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]); %4.6
    end
    
    legend({'OFF mean','ON mean','OFF fit','ON fit','OFF pref','ON pref'}, ...
        'Location','northoutside','Orientation','horizontal');
    set(gcf,'Units','centimeters','Position',[10 10 17 4]);
    
    PlotTitle='_ von Mises Example fits ';
basePath = fullfile(saveDir_FigS3,[ PlotTitle]);
saveFigAll(gcf, basePath);
    
end


%% Plot
color_shuff_face = [0.85 0.85 0.85];
color_shuff_edge = [0.3  0.3  0.3 ];

figure;

edgesMu   = linspace(-200,   200,   50);   % orientation: (-pi/2, pi/2] in deg
edgesFWHM = linspace(-100, 100, 30);


for F = 1:2
    row = 3-F;
    % --- delta Preferred direction ---
    subplot(2,2,(row-1)*2+1); hold on
    histogram(dMu_shuff{F}, edgesMu, 'Normalization','probability', ...
        'FaceColor', color_shuff_face, 'EdgeColor', color_shuff_edge);
    histogram(dMu{F},       edgesMu, 'Normalization','probability', ...
        'FaceColor', color_OFF{F},    'EdgeColor', Params.colors{F});
    yl = ylim;
    plot([0 0], yl, 'k:', 'HandleVisibility','off');
    title('\DeltaPref. direction (deg)'); box off
    %xlim([-90 90])
    set(gca, 'Fontsize', 12,'linewidth',1)
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
    
    % --- delta tuning width (deg) ---
    subplot(2,2,(row-1)*2+2); hold on
    histogram(dFWHM_shuff{F}, edgesFWHM, 'Normalization','probability', ...
        'FaceColor', color_shuff_face, 'EdgeColor', color_shuff_edge);
    histogram(dFWHM{F},       edgesFWHM, 'Normalization','probability', ...
        'FaceColor', color_OFF{F},    'EdgeColor', Params.colors{F});
    yl = ylim;
    plot([0 0], yl, 'k:', 'HandleVisibility','off');
    title('\DeltaTuning width (deg)'); box off
    xlabel('Degrees')
    set(gca, 'Fontsize', 12,'linewidth',1)
    ax = gca; ax.XColor = 'k'; ax.YColor = 'k';
    set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
   
end

set(gcf,'Units','centimeters','Position',[10 10 6 4]);

PlotTitle='_Fit von Mises pref dir and tuning width ';
basePath = fullfile(saveDir_FigS3,[ PlotTitle]);
saveFigAll(gcf, basePath);
