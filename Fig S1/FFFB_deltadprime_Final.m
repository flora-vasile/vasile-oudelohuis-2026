saveDir_FigS1 = [Params.MainDir filesep 'plot\Fig S1'];
if ~exist(saveDir_FigS1,'dir')
    mkdir(saveDir_FigS1)
end


for F=1:2
    
    for ani=1:Params.n_ani{F}
        for pos=1:Params.n_pos{F}(ani)
            if ~isempty(data{F}{ani,pos})
                ROICoordyxResp{F}{ani,pos}=CellStat{F}{ani,pos}.Coordyx_um(CellStat{F}{ani,pos}.responsiveAD,:);
                %Modulation for responsive cells
                AllDirRespAllRepSes{F}{ani,pos}=squeeze(mean(data{F}{ani,pos}.StimAv(:,Params.StimONResponse,:,:),2));
                AllDirRespSes{F}{ani,pos}=squeeze(mean(mean(data{F}{ani,pos}.StimAv(:,Params.StimONResponse,:,:),2),4));
                
                dprime{F}{ani,pos}=(AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,2)...
                    -AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,1))...
                    ./sqrt(1/2*(var(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,2,:),0,3)...
                    +var(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,1,:),0,3)));
                
                
                % Number of points
                num_points = size(ROICoordyxResp{F}{ani,pos}, 1);
                
                % Initialize an empty matrix to store the Euclidean distances
                distances{F}{ani,pos} = nan(num_points, num_points);
                mod_diff{F}{ani,pos} = nan(num_points, num_points);
                % Calculate the Euclidean distance between each pair of points
                for i = 1:num_points
                    for j = 1:num_points
                        if i ~= j
                            % Compute the Euclidean distance: sqrt((x2-x1)^2 + (y2-y1)^2)
                            dist = sqrt((ROICoordyxResp{F}{ani,pos}(i,1) - ROICoordyxResp{F}{ani,pos}(j,1))^2 + (ROICoordyxResp{F}{ani,pos}(i,2) - ROICoordyxResp{F}{ani,pos}(j,2))^2);
                            distances{F}{ani,pos}(i,j) = dist;
                            % Compute the absolute difference in Modulation values
                            mod_diff{F}{ani,pos}(i,j) = abs(dprime{F}{ani,pos}(i) - dprime{F}{ani,pos}(j));
                            
                        else
                            distances{F}{ani,pos}(i,j) = 0;  % Distance from a point to itself is 0
                            mod_diff{F}{ani,pos}(i,j) = 0;    % No modulation difference for the same point
                        end
                    end
                end
                
                % Define bin size for distances
                bin_size = 10; %um
                maxbin=300;
                binsin=[1:maxbin/bin_size];
                xv=[0:bin_size:length(binsin)*bin_size];

                % Create bins for distances
                bin_edges = 0:bin_size:maxbin;
                
                % Bin the distances and calculate the average modulation difference
                for i = 1:num_points
                    mod_diff_all_bins{F}{ani,pos}{i} = cell(length(bin_edges)-1, 1);
                    bin_counts{i} = zeros(length(bin_edges)-1, 1);
                    for j = 1:num_points
                        if i ~= j & distances{F}{ani,pos}(i,j)<=maxbin
                            % Find the bin index for the distance
                            dist = distances{F}{ani,pos}(i,j);
                            bin_index = find(dist >= bin_edges, 1, 'last');
                            % Add the modulation difference to the corresponding bin
                            mod_diff_all_bins{F}{ani,pos}{i}{bin_index} = [mod_diff_all_bins{F}{ani,pos}{i}{bin_index}, mod_diff{F}{ani,pos}(i,j)];
                            bin_counts{i}(bin_index) = bin_counts{i}(bin_index) + 1;
                        end
                    end
                end
                
                % Calculate the average modulation difference for each bin and each neuron
                for i = 1:length(mod_diff_all_bins{F}{ani,pos})
                    mod_diff_binned{F}{ani,pos}{i} = zeros(length(bin_edges)-1, 1);
                    for k = 1:length(mod_diff_all_bins{F}{ani,pos}{i})
                        if bin_counts{i}(k) > 0
                            % Calculate the mean modulation difference for the bin
                            mod_diff_binned{F}{ani,pos}{i}(k) = mean(mod_diff_all_bins{F}{ani,pos}{i}{k});
                        else
                            mod_diff_binned{F}{ani,pos}{i}(k) = NaN;
                        end
                    end
                end
                
                Mat2{F}{ani,pos}=cat(2,mod_diff_binned{F}{ani,pos}{:});  %dim1: neurons, dim2: bins
                Ses{F}{ani,pos}=nanmean(Mat2{F}{ani,pos},2);
            end
        end
    end

    meanSes{F}=cat(2,Ses{F}{:})';
    meanAllSes{F}=nanmean(meanSes{F},1);
    
end


%% One way rm anova
% Find rows (sessions) without NaNs in each group
validRowsG1 = all(~isnan(meanSes{1}(:,binsin)), 2);
validRowsG2 = all(~isnan(meanSes{2}(:,binsin)), 2);

% Filter out animals with NaNs
FB_clean = meanSes{1}(validRowsG1, binsin);
FF_clean = meanSes{2}(validRowsG2, binsin);

semAllSes{1}=nanstd(FB_clean,1)./sqrt(size(FB_clean,1));
semAllSes{2}=nanstd(FF_clean,1)./sqrt(size(FF_clean,1));

% Normality test
for i=1:size(FB_clean,2)   
Normtest{1}(i)=kstest(FB_clean(:,i));  
end
for i=1:size(FF_clean,2)   
Normtest{2}(i)=kstest(FF_clean(:,i));   
end


[pFB,tbl] = friedman(FB_clean,1,'off');
fprintf('p = %.16e\n', pFB)
disp(tbl)

[pFF,tbl] = friedman(FF_clean,1,'off');
fprintf('p = %.16e\n', pFF)
disp(tbl)

%% Figure
figure;
F=1;
p1=plot(xv(binsin), mean(FB_clean,1), 'LineWidth', 1, ...
    'MarkerSize', 4, 'MarkerEdgeColor', Params.colors{F}, 'MarkerFaceColor', Params.colors{F},'Color', Params.colors{F}); hold on

h1 = fill([xv(binsin) flip(xv(binsin))],[mean(FB_clean,1)-semAllSes{F}(binsin) flip(mean(FB_clean,1)+semAllSes{F}(binsin))],Params.colors{F},'linestyle','none');
set(h1,'facealpha',.3)
F=2;
p2=plot(xv(binsin),mean(FF_clean,1), 'LineWidth', 1, ...
    'MarkerSize', 4, 'MarkerEdgeColor', Params.colors{F}, 'MarkerFaceColor', Params.colors{F},'Color', Params.colors{F}); hold on

h1 = fill([xv(binsin) flip(xv(binsin))],[mean(FF_clean,1)-semAllSes{F}(binsin) flip(mean(FF_clean,1)+semAllSes{F}(binsin))],Params.colors{F},'linestyle','none');
set(h1,'facealpha',.3)
xlabel('Distance between cell pairs (\mum)');
xl=xlim;
plot(xl,[1 1],'k:');
xlim([0 length(binsin)*bin_size]);
ylim([0.1 0.35])
ylabel('|delta d prime| ');
h=[p1;p2];
lgd=legend(h,'FF','FB', 'Location','southeast');  legend boxoff
box off

set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
despline(gca);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';


set(gcf,'Units','centimeters','Position',[10 10 4.5 3.5]);

PlotTitle='_delta dprime';
basePath = fullfile(saveDir_FigS1,[ PlotTitle]);
saveFigAll(gcf, basePath);

