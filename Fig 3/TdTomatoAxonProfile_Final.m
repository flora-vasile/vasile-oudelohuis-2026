%% Input data from ImageJ
FF=1;
FB=0;

if FF
    colorTrace=2;
elseif FB
    colorTrace=1;
end


%% Normalise to cortical depth
EndnValues=200;
DepthVect=[1:EndnValues];


for ani= 1:size(FtdTprofileImgJ,2)
for n=1:size(FtdTprofileImgJ{ani},2)
    temp=downsample(movmean(FtdTprofileImgJ{ani}{n}(:,2),10),floor(length(FtdTprofileImgJ{ani}{n})./EndnValues));
FtdTprofileds{ani}(:,n)=temp(1:EndnValues);
end
end

%% Save Matrix with tdTomato profile
SaveDir = ['D:\Confocal data spectral unmixing\ tdTomato axon profile'];

if ~exist(SaveDir,'dir')
    mkdir(SaveDir)
end

if FF
save([SaveDir filesep tv 'TdTomatoFluoProfileFF']);
elseif FB
save([SaveDir filesep tv 'TdTomatoFluoProfileFB']);    
end

%% Normalise to max
for ani= 1:size(FtdTprofileds,2)
for n=1:size(FtdTprofileds{ani},2)
    FtdTprofileNorm{ani}(:,n)=FtdTprofileds{ani}(:,n)./max(FtdTprofileds{ani}(:,n));
end
end

%% Mean for each animal

for ani= 1:size(FtdTprofileds,2)
meanSliceFtdTprofileNorm(:,ani)=mean(FtdTprofileNorm{ani},2);
end

%% sem
semSliceFtdTprofileNorm=std(meanSliceFtdTprofileNorm,0,2)./sqrt(size(FtdTprofileds,2));

%% Plot Profile

figure;
for ani= 1:size(FtdTprofileds,2)
plot( meanSliceFtdTprofileNorm(:,ani),DepthVect, 'Linewidth', 0.4, 'Color',[0.6 1 1]); hold on
end
plot( mean(meanSliceFtdTprofileNorm,2),DepthVect, 'Linewidth', 1, 'Color',[0.4 1 1]);
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', fliplr(yt))
xlabel('tdTomato Fluorescence (Norm)'); ylabel('Cortical Depth (a.u.)'); 
box off
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
offsetAxes(gca);
yticks([])
set(gca, 'Fontsize', 6,'LineWidth', 0.8,'TickDir','out','TickLength',[0.02, 0.01])
set(gcf,'Units','centimeters','Position',[10 10 0.8 4])

if FF
    PlotTitle = 'FF Axon Layer distribution';
else
    PlotTitle = 'FB Axon Layer distribution';
end
h=suptitle(titlefig);
set(h,'FontSize',7)

basePath = fullfile(saveDir_FigS5,[PlotTitle]);
saveFigAll(gcf, basePath);
