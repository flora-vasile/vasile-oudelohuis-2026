function HeatMapPlot(Input,ylab,xlab,titl,rangecolorbar)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
imagesc(Input); hold on
set(gca,'YDir','reverse') 
ylabel(ylab, 'Fontsize', 15)
xlabel(xlab, 'Fontsize', 15)
title(titl, 'Fontsize', 15)

a=colorbar;
caxis(rangecolorbar);
colormap(jet);
a.Label.String = 'Evenst/F0';
%ylabel(a,'Laser ON - Laser OFF','FontSize',16,'Rotation',270);
set(gca, 'Fontsize', 12);
if size(Input,2)<25
xticks([0:6:24]);
plot([6.5 6.5],[0 size(Input,1)],'k--', 'lineWidth', 1);
plot([12.5 12.5],[0 size(Input,1)],'k--', 'lineWidth', 1);
else
xticks([0:15:60]); 
plot([15.5 15.5],[0 size(Input,1)],'k--', 'lineWidth', 1);
plot([30.5 30.5],[0 size(Input,1)],'k--', 'lineWidth', 1);
end
xticklabels([0:1:4]);

xlim([0.5 size(Input,2)+0.5]); ylim([0.5 size(Input,1)+0.5]);

end

