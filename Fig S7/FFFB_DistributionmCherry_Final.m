
saveDir_FigS7 = fullfile('E:\2P\Pop Analysis\Final_datasets\LoopV1PM_deconv_final\plot', 'plot', 'Fig S7');
if ~exist(saveDir_FigS7, 'dir')
    mkdir(saveDir_FigS7)
end

FF=1;
%% Distribution in Fraction

for ani= 1:size(CoordmCherryCellsNorm,2)
    for n=1:size(CoordmCherryCellsNorm{ani},2)
NumberDepth{ani}(:,n)=histcounts(CoordmCherryCellsNorm{ani}{n},nbins);
FractionDepth{ani}(:,n)= NumberDepth{ani}(:,n)./size(CoordmCherryCellsNorm{ani}{n},1);  
    end
end

%% Mean for each animal

for ani= 1:size(FractionDepth,2)
meanSliceFractionDepth(:,ani)=mean(FractionDepth{ani},2);
end

%% sem
semSliceFractionDepth=std(meanSliceFractionDepth,0,2)./sqrt(size(FractionDepth,2));

%% Plot Profile

figure;
for ani= 1:size(FractionDepth,2)
plot( meanSliceFractionDepth(:,ani),DepthBin, 'Linewidth', 0.4, 'Color',[1 0.4 0.4]); hold on
end
plot( mean(meanSliceFractionDepth,2),DepthBin, 'Linewidth', 1, 'Color','r');
yt = get(gca, 'YTick');
set(gca, 'YTickLabel', fliplr(yt))
xlabel('mCherry Cells Fraction'); ylabel('Cortical Depth (a.u.)'); 
box off
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
offsetAxes(gca);
yticks([])
set(gca, 'Fontsize', 6,'LineWidth', 0.8,'TickDir','out','TickLength',[0.02, 0.01])
set(gcf,'Units','centimeters','Position',[10 10 0.8 4])

if FF
PlotTitle='Distribution of mCherry Cells in PM' ;
else
PlotTitle='Distribution of mCherry Cells in V1' ;
end
basePath = fullfile(saveDir_FigS7,[PlotTitle]);
saveFigAll(gcf, basePath);


