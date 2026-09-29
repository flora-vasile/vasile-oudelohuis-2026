saveDir_FigS4 = fullfile(Params.MainDir, 'plot', 'Fig 4');
if ~exist(saveDir_FigS4,'dir')
    mkdir(saveDir_FigS4)
end
for F=1:2      
for C=1:2      
 
AllDirsResp{F}{C}=squeeze(mean(mean(data{F}{C}.StimAv(:,Params.StimONResponse,:,:),2),4));
AllDirsRespReps{F}{C}=squeeze(mean(data{F}{C}.StimAv(:,Params.StimONResponse,:,:),2));

dprimeSes{F}{C}=(AllDirsResp{F}{C}(CellStat{F}{C}.responsiveAD,2)-AllDirsResp{F}{C}(CellStat{F}{C}.responsiveAD,1))...
    ./sqrt(1/2*(var(AllDirsRespReps{F}{C}(CellStat{F}{C}.responsiveAD,2,:),0,3)+var(AllDirsRespReps{F}{C}(CellStat{F}{C}.responsiveAD,1,:),0,3)));
end
end

figure; hold on
for F=1:2
subplot(2,1,3-F); hold on
IndexHist(ModAllDirs{F}{1}, ModAllDirs{F}{2}, Params.colors{F}, [0.5 0.5 0.5], ' ','ChrimsonR','Control','d prime',0.02)
xlim([-0.6 0.6])

 set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
end
    box off

set(gcf,'Units', 'centimeters','Position',[10 10 3 6.5])

PlotTitle='_dprime Chrimson vs Control';
basePath = fullfile(saveDir_FigS4,[PlotTitle]);
saveFigAll(gcf, basePath);
