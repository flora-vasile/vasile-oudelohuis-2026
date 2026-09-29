%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_ROI = 1;
foldername = 'Z:\Users\Devanshi\Flora Muscimol Experiments\FF196\FF196_FF_310726_Muscimol2'; % folder where all the files are located

DataID = 'FF196_pos1_muscimol2';
skipNofFrames = 1; % skip the extraction of number of frames from tiff files, coz it's an input

isFlyBack = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parametrize
s2p_folder = dir([foldername filesep 'suite2p']);
match = {}; counter = 1;
for i = 1:length(s2p_folder)
    if ~isempty(regexp(s2p_folder(i).name,'^plane'))
        match{end+1} = s2p_folder(i).name;
        index(counter) = i;
        counter = counter +1;
    end
end
if length(index)>1
    index = index(1:end-isFlyBack);
    plane = 'multiples';
else
    plane = 'single';
end
disp(['number of planes: ', num2str(length(index))])

dt = datestr(now,'yyyymmdd');

filetype = 'tif'; % type of files to be processed
saveFolder = [foldername '\Deconv\analysis'];
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

oldfolder = cd(saveFolder);
diaryName = ['FtracesLog_' dt];
diary off; diary(diaryName)
cd(oldfolder)
saveName = ['SpksTrace_' DataID '_' dt];
disp(saveFolder); disp(saveName);

%% Chop all the Fall file of each plane
for j = 1:length(index)
    path = [s2p_folder(index(j)).folder filesep s2p_folder(index(j)).name];
    disp(s2p_folder(index(j)).name)
    load([path '\Fall.mat'])
    
    % Y and X regstration offset
    yOffset = ops.yoff1;
    xOffset = ops.xoff1;
    % list of files
    fileList = ops.filelist;
    
    if skipNofFrames==0
        files = subdir(fullfile(foldername,[GenericFile_name,filetype]));   % list of filenames (will NOT search all subdirectories)
        [H, Data] = opentif(files(1).name); %opentif_SI2019b(files(1).name);
        n2pFiles = length(files);
        disp(['Found ',num2str(n2pFiles), ' tif files'])
        numberOfFrames = zeros(n2pFiles,1);
        f =  waitbar(0);
        for i = 1:n2pFiles
            info = imfinfo(files(i).name);
            numberOfFrames(i) = numel(info)/(size(Data,3)*(size(Data,5)));
            waitbar(i/n2pFiles,f,[num2str(i/n2pFiles*100,2),'%'])
        end
        close(f)
    else
        try
    numberOfFrames = ops.frames_per_file;
        catch
            numberOfFrames = double(ops.nframes/length(ops.filelist))*ones(1,length(ops.filelist));
        end
    end
    
    disp(['number of frames per 2P files btwn ' num2str(max(numberOfFrames)) ' and ' num2str(min(numberOfFrames))])
    
    
    if all_ROI == 1
        Deconvtraces_all{j}.spks = spks;
        Deconvtraces_all{j}.RawF = F;
    else
        Deconvtraces_all{j}.spks = spks(is_cell(1:lastROI),:);
    end
    
    Deconvtraces_all{j}.s2p.iscell   = iscell(:,1);
    Deconvtraces_all{j}.s2p.meanImg  = ops.meanImg;
    Deconvtraces_all{j}.s2p.meanImgE = ops.meanImgE;
    Deconvtraces_all{j}.s2p.stat     = stat;
    Deconvtraces_all{j}.s2p.yOffset  = yOffset;
    Deconvtraces_all{j}.s2p.xOffset  = xOffset;
    
    
    try
        Deconvtraces_all{j}.header.frameFrequency=H.SI.hRoiManager.scanFrameRate;
    catch
        if j == 1
        indx = listdlg('PromptString','Sampling Frequency, per plane','SelectionMode','single',...
            'ListString',{'30','15','10','7.5','6'});
        frameFrequency =  30.4812/indx; %frameFrequency =  60.9624/indx;
        disp(['Sampling Frequency = ' num2str(frameFrequency)]);
        end
    Deconvtraces_all{j}.header.frameFrequency = frameFrequency;
        
        %    15.2406;
    end
    Deconvtraces_all{j}.header.numberOfFrames = numberOfFrames;
    
    Deconvtraces_all{j}.Info.saveName = DataID;
    Deconvtraces_all{j}.Info.mainDir = foldername;
    Deconvtraces_all{j}.Info.fileList = fileList;
    Deconvtraces_all{j}.Info.plane = plane;
end


%% save
fileName = [saveFolder,'\',saveName,'.mat'];
save(fileName,'Deconvtraces_all')
diary off