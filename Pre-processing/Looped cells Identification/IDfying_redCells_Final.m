%% Input: x and y vector with coordinates of red cells
% data.mat file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   Preambule   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if paste measurments from ImageJ in Matlab as IJ
for j = 1:data.Info.nPlanes
    if ~isempty(IJ{j})
        x = IJ{j}(:,3);
        y = IJ{j}(:,4);
    else
        x=[];
        y=[];
    end
    saveFig = true;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% get infos
    if isfield(data.s2p{j},'rescale_factor')
        rescale_factor = data.s2p{j}.rescale_factor;
    else
        answer = inputdlg('structural stack-to-functional imaging rescaling factor');
        rescale_factor = answer{1,1}; %  ~=1 if z-stack and functional recording have been acquired with different resolution
    end
    if saveFig
        mainDir =  data.Info.mainDir;
        if exist(mainDir,'dir')
            disp(mainDir)
        else
            mainDir = uigetdir('X:\camille.mazo\2P_processed','Remap main dir');
            fprintf('Main directory has changed to %s \n',mainDir)
            data.Info.mainDir = mainDir;
        end
        saveDir = [data.Info.mainDir, '\Figures\Bead_Pos'];
        if ~ exist(saveDir)
            mkdir(saveDir)
        end
    end
    
    %% fill the s2p {iscell=1} masks in a matrix
    for i = 1:size(data.s2p{j}.stat,2)
        if data.s2p{j}.iscell(i)
            size_masks(i) = size(data.s2p{j}.stat{1, i}.xpix,2);
        end
    end
    size_masks_max = max(size_masks);
    cellMask_x = NaN(size(data.s2p{j}.stat,2),size_masks_max);
    cellMask_y = NaN(size(data.s2p{j}.stat,2),size_masks_max);
    for i = 1:size(data.s2p{j}.stat,2)
        if data.s2p{j}.iscell(i)
            cellMask_x(i,1:size_masks(i)) =  double(data.s2p{j}.stat{1, i}.xpix);
            cellMask_y(i,1:size_masks(i)) =  double(data.s2p{j}.stat{1, i}.ypix);
        end
    end
    
    %% plot the masks on the target image
    target = data.s2p{j}.meanImg;
    if rescale_factor ~=1
        cellMask_x = round(cellMask_x/2);
        cellMask_y = round(cellMask_y/2);
        target = imresize(target,rescale_factor);
    end
    
    BeadCell_Masks{j} = figure;
    hold on
    imagesc(target)
    axis equal
    xlim([0 size(target,1)]); ylim([0 size(target,2)]);
    axis ij
    colormap('gray')
    
    %% plot the bead+ cells
    x= round(x); y=round(y);
    plot(x,y,'r+')
    
    test = [];
    counter = 1;
    for i = 1:length(x)
        cellID = any(x(i) == cellMask_x,2) & any(y(i) == cellMask_y,2);
        if sum(cellID)
            try
                test(counter,1) = find(cellID==1);
            catch
                ROIlist = find(cellID);
                plot(x(i),y(i),'+', 'Color',[1 1 1])
                f1 = figure; hold on;
                imagesc(target); axis equal
                xlim([0 size(target,1)]); ylim([0 size(target,2)]);
                axis ij; colormap('gray')
                for ii = 1:size(ROIlist,1)
                    plot(cellMask_x(ROIlist(ii),:),cellMask_y(ROIlist(ii),:),'k');
                    text(cellMask_x(ROIlist(ii),1),cellMask_y(ROIlist(ii),1),num2str(ROIlist(ii)),...
                        'Color',[1 1 1],'HorizontalAlignment','Center')
                end
                plot(x(i),y(i),'+', 'Color',[1 1 1])
                disp(num2str(ROIlist))
                test(counter,1) = input('Which ROI is bead+?\n')
                close(f1)
            end
            test(counter,2) = true;
            plot(cellMask_x(test(counter,1),:),cellMask_y(test(counter,1),:),'r');
            counter = counter+1;
        else
            disp(['beads+ #',num2str(i),' is not an active cell'])
        end
        
    end
    plot(x,y,'+','Color',[1 1 1])
    title({['Plane' num2str(j)],'Red mask: active & bead+ cell', 'White cross: structurally identified bead+ cell'})
    
    %% logical indexing with ROI masks matching the number of "iscell" ROIs
    if ~isempty(test)
        for i = 1:size(test,1)
            ROInum = test(i,1);
            data.beads_pos{j}(i,1) = sum(data.s2p{j}.iscell(1:ROInum));
        end
    else
        data.beads_pos{j}=[];
    end
end % end loop through planes

selec = cell(1,data.Info.nPlanes);
for j = 1:data.Info.nPlanes
    selec{j} = false(1,sum(data.s2p{j}.iscell));
    selec{j}(data.beads_pos{j}) = true;
end
selec = cat(2,selec{:});
%% save...
%  - the fig
for j = 1:data.Info.nPlanes
    figure(BeadCell_Masks{j})
    screen2png([saveDir '\MeanImg_Mask_p' num2str(j)]);
end

%  - the data
data.beadID.bead_pos = selec;
data.beadID.x = x;
data.beadID.y = y;
data.beadID.IJmeasurments = IJ; %[area mean min max x y ch]

if exist(data.Info.mainDir,'dir') == 7 && isfield(data.Info,'saveName')
    mainDir = data.Info.mainDir;
    saveName = data.Info.saveName;
elseif isfield(data.Info,'saveName')
    mainDir = uigetdir('E:\TempData\','select dir to save data');
else
    mainDir = uigetdir('E:\TempData\','select dir to save data');
    answer = inputdlg('save name (animal_position)');
    saveName = answer{1,1};
end
try
save([mainDir, '\Deconv\analysis\Spks_', saveName],'data','-append')
catch
save([mainDir, '\analysis\', saveName],'data','-append')
end

disp('data saved');