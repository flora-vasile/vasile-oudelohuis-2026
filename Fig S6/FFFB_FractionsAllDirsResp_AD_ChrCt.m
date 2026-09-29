saveDir_FractionAllDirs{1} = [MainDir filesep 'plot\FB\Fraction All Dirs'];
saveDir_FractionAllDirs{2} = [MainDir filesep 'plot\FF\Fraction All Dirs'];

colors = [
    0.8 0.8 0.8;    % grey
    1 0.6 0.4;
    0.6 0.8 1]; % blue

%1 0.6 0.6;    % red

FFFBName{1}='Feedback'; FFFBName{2}='Feedforward';

for F=1:2      
    if ~exist(saveDir_FractionAllDirs{F},'dir')
        mkdir(saveDir_FractionAllDirs{F})
    end
    
       
for C=1:2
nROIs{F}{C}=sum(ResponsiveAD{F}{C});
FracExc{F}{C}=(sum(ExcAD{F}{C})./nROIs{F}{C})*100;
FracInh{F}{C}=(sum(InhAD{F}{C})./nROIs{F}{C})*100;
FracNoChange{F}{C}=((nROIs{F}{C}-sum(ExcAD{F}{C})-sum(InhAD{F}{C}))./nROIs{F}{C})*100;
end

labels={'No change','Excited','Inhibited'};

% z test Exc
zstat = myztest(sum(ExcAD{F}{1}),sum(ExcAD{F}{2}),nROIs{F}{1},nROIs{F}{2});
pzExc{F} =  pvaluefromz(zstat);
% z test Inh
zstat = myztest(sum(InhAD{F}{1}),sum(InhAD{F}{2}),nROIs{F}{1},nROIs{F}{2});
pzInh{F} =  pvaluefromz(zstat);
end

% bars 
    figure;
    for F=1:2
    subplot(2,1,3-F); hold on
    b1 = bar([FracExc{F}{1} FracExc{F}{2}; FracInh{F}{1} FracInh{F}{2}],'FaceColor','flat','LineWidth',0.7,'barWidth',0.8);hold on
    
    b1(1).CData = [1 0.6 0.4;[0.6 0.8 1]]; %b1(2).CData = [0.6 0.8 1]; %b1(3).CData =[0.6 0.8 1];  b1(4).CData = [0.8 0.8 0.8];
    
%     legend({b1(1),b1(2)},'ChrimsonR+light','Light'); legend boxoff
    
    xticks([0.9 1.1 1.9 2.1]); xticklabels({'ChrimsonR','Control','ChrimsonR','Control'});
    xtickangle(20)
    
        text(1, max([FracExc{F}{1} FracInh{F}{1}]+1) ,'Excited','FontSize',6,...
            'HorizontalAlignment','Center','VerticalAlignment','bottom')
        text(2, max([FracExc{F}{1} FracInh{F}{1}]+1) ,'Inhibited','FontSize',6,...
            'HorizontalAlignment','Center','VerticalAlignment','bottom')
        
%         text(F, 92 ,[num2str(sum(ResponsiveAD{F}))], 'fontweight', 'bold' ,'FontSize',18,'color','w',...
%             'HorizontalAlignment','Center','VerticalAlignment','bottom')

    
    ylabel('Percentage');
         xlim([0.5 2.5]);
         
         set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
% despline(gca);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
end
    box off

set(gcf,'Units', 'centimeters','Position',[10 10 3 7])
for F=1:2      
saveas(gcf,[saveDir_FractionAllDirs{F} filesep tv '_Fraction modulated cells Chrimson vs Control '],'fig');
saveas(gcf,[saveDir_FractionAllDirs{F} filesep tv '_Fraction modulated cells Chrimson vs Control '],'svg');
screen2png([saveDir_FractionAllDirs{F} filesep tv '_Fraction modulated cells Chrimson vs Control ']);

end



