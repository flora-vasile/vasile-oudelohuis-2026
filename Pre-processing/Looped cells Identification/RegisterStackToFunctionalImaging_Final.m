function data = RegisterStackToFunctionalImaging_v2(data, crop,regStack, varargin)
global j
for j = 1:data.Info.nPlanes
    
%% load the images (/ stacks)
try
    target = data.s2p{j}.meanImg;
catch
    try
        target = data.s2p.meanImg;
    catch
        [target_filename, path_target] = uigetfile(['X:\camille.mazo\2P_processed\','*.tif'],...
            'Select target for reg');
        target = imread([path_target target_filename]);
    end
end

if ~isempty(varargin)
    path = varargin{1};
    [source_filename, path] = uigetfile([path 'Processed' filesep '*.tif'],...
        ['select image (stack) to be registered, channel 1 plane' num2str(j)]);
else
    [source_filename, path] = uigetfile(['X:\camille.mazo\2P_processed\','*.tif'],...
        'select image to be registered, channel 1');
end
info = imfinfo([path source_filename]);
num_images = numel(info);

% rescale the target if functional imaging and structural stack have
% different resolution
Width = info.Width; Height = info.Width;
if size(target) ~= [Width Height]
    if Width == Height && size(target,1) == size(target,2)
        rescale_factor = Width/size(target,1);
        target = imresize(target,rescale_factor);
    else
        disp('not a square?1?')
    end
else
    rescale_factor = 1;
end

% crop the target
target_crop = target(crop(1):512-2*crop(1),crop(2):512-2*crop(2));
fft_target = fft2(target_crop);


%% find best plane by looking at the error in each
% load the planes in structural z-stack and register to the target
for i = 1:num_images
    source = imread([path source_filename], i ,'Info', info);
    source_crop = source(crop(1):512-2*crop(1),crop(2):512-2*crop(2));
    fft_source = fft2(source_crop);
    [output, ~] = dftregistration(fft_target,fft_source,1);
    outputs(i,:) = output;
end
% find best plane = lower error.
[~,b] = min(outputs,[],1);
source = imread([path source_filename], b(1) ,'Info', info);

fprintf('Best plane is %g \nShift in x and y are: %g, %g\n' , b(1), outputs(b(1),3), outputs(b(1),4))
fileName = [path,'\reg_metrics_p', num2str(j), '.mat'];
save(fileName,'outputs')


%% load corresponding image in channel 2
try % Using "NewAlignment"
    source_filename_c2 = [source_filename(1:end-5), '2.tif'];
    path_2 = path;
    info_2 = imfinfo([path_2 source_filename_c2]);
catch
    try
        source_filename_c2 = [source_filename(1:end-8), '2',source_filename(end-6:end)];
        path_2 = path;
        info_2 = imfinfo([path_2 source_filename_c2]);
    catch % manual
        [source_filename_c2, path_2] = uigetfile({[path '\*.*']},'select image to be registered, channel 2');
        info_2 = imfinfo([path_2 source_filename_c2]);
    end
end
source_c2 = imread([path_2 source_filename_c2], b(1) ,'Info', info_2);


%% apply x and y shifts to the z-projection

max_shift = max(abs(outputs(b(1),3)), abs(outputs(b(1),4)));
source_reg = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
source_reg(max_shift+1+outputs(b(1),3):size(source,1)+max_shift+1+outputs(b(1),3)-1,max_shift+1+outputs(b(1),4):size(source,1)+max_shift+1+outputs(b(1),4)-1) = source;

if max_shift==0 %% bug fix 230201
    max_shift=1;
end

source_reg_crop = source_reg(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1);

source_c2_reg = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
source_c2_reg(max_shift+1+outputs(b(1),3):size(source,1)+max_shift+1+outputs(b(1),3)-1,max_shift+1+outputs(b(1),4):size(source,1)+max_shift+1+outputs(b(1),4)-1) = source_c2;
source_c2_reg_crop = source_c2_reg(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1);

% for display purposes
a = min(min(source_reg_crop)) ;
isnan(source_reg_crop);
source_reg_crop(ans) = a;%-10;

a = min(min(source_c2_reg_crop)) ;
isnan(source_c2_reg_crop);
source_c2_reg_crop(ans) = a;%-10;

% save the images
source_reg_crop = int16(source_reg_crop);
source_c2_reg_crop = int16(source_c2_reg_crop);
target = int16(target);
clear source_c1c2_reg_crop
source_c1c2_reg_crop(:,:) = source_reg_crop;
source_c1c2_reg_crop(:,:,2) = source_c2_reg_crop;

filename_c1 = [path source_filename(1:end-4) '_1plAVGreg.tif'];
filename_c2 = [path_2 source_filename_c2(1:end-4) '_1plAVGreg.tif'];
filename_target = [path 'target_p' num2str(j) '.tif'];
filename_c1c2 = [path source_filename(1:end-4) 'c2' ,source_filename(end-6:end-4) '_1plAVGreg.tif'];

tt=Tiff(filename_c1,'w');
tt2=Tiff(filename_c2,'w');
tt3=Tiff(filename_target,'w');
tt4=Tiff(filename_c1c2,'w');

tagstruct.ImageLength = size(source_reg_crop,1);
tagstruct.ImageWidth = size(source_reg_crop,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat=Tiff.SampleFormat.Int;

tt.setTag(tagstruct);
tt.write(source_reg_crop);
tt.writeDirectory()
close(tt)

tt2.setTag(tagstruct);
tt2.write(source_c2_reg_crop);
tt2.writeDirectory()
close(tt2)

tt3.setTag(tagstruct);
tt3.write(target);
tt3.writeDirectory()
close(tt3)

for ii = 1:size(source_c1c2_reg_crop,3)
    tt4.setTag(tagstruct);
    tt4.write(source_c1c2_reg_crop(:,:,ii));
    tt4.writeDirectory()
end
close(tt4)

%% Average 3 z-planes
if num_images >1
    nPlaneToAvg = 3;
    if  b(1)~= 1 &&  b(1)~=num_images
        for i = 1:nPlaneToAvg
            source(:,:,i) = imread([path source_filename], b(1)-2+i ,'Info', info);
            source_c2(:,:,i) = imread([path_2 source_filename_c2], b(1)-2+i ,'Info', info_2);
        end
    elseif b(1)== 1
        for i = 1:nPlaneToAvg
            source(:,:,i) = imread([path source_filename], b(1)+i-1 ,'Info', info);
            source_c2(:,:,i) = imread([path_2 source_filename_c2], b(1)+i-1 ,'Info', info_2);
        end
    elseif b(1)== num_images
        for i = 1:nPlaneToAvg
            source(:,:,i) = imread([path source_filename], b(1)-3+i ,'Info', info);
            source_c2(:,:,i) = imread([path_2 source_filename_c2], b(1)-3+i ,'Info', info_2);
        end
    end
    source = mean(source,3); % max(source,[],3);
    source_c2 = mean(source_c2,3); % max(source_c2,[],3);
    
    % max_shift = max(abs(outputs(b(1),3)), abs(outputs(b(1),4)));
    source_reg = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
    source_reg(max_shift+1+outputs(b(1),3):size(source,1)+max_shift+1+outputs(b(1),3)-1,max_shift+1+outputs(b(1),4):size(source,1)+max_shift+1+outputs(b(1),4)-1) = source;
    source_reg_crop = source_reg(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1);
    
    source_c2_reg = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
    source_c2_reg(max_shift+1+outputs(b(1),3):size(source,1)+max_shift+1+outputs(b(1),3)-1,max_shift+1+outputs(b(1),4):size(source,1)+max_shift+1+outputs(b(1),4)-1) = source_c2;
    source_c2_reg_crop = source_c2_reg(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1);
    
    
    % for display purposes
    a = min(min(source_reg_crop)) ;
    isnan(source_reg_crop);
    source_reg_crop(ans) = a;%-10;
    
    a = min(min(source_c2_reg_crop)) ;
    isnan(source_c2_reg_crop);
    source_c2_reg_crop(ans) = a;%-10;
    
    % save the images
    source_reg_crop = int16(source_reg_crop);
    source_c2_reg_crop = int16(source_c2_reg_crop);
    
    filename_c1 = [path source_filename(1:end-4) '_3plAVGreg.tif'];
    filename_c2 = [path source_filename_c2(1:end-4) '_3plAVGreg.tif'];
    filename_c1c2 = [path source_filename(1:end-4) 'c2' ,source_filename(end-6:end-4) '_3plAVGreg.tif'];
    
    tt=Tiff(filename_c1,'w');
    tt2=Tiff(filename_c2,'w');
    tt4=Tiff(filename_c1c2,'w');
    
    tagstruct.ImageLength = size(source_reg_crop,1);
    tagstruct.ImageWidth = size(source_reg_crop,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.SampleFormat=Tiff.SampleFormat.Int;
    
    tt.setTag(tagstruct);
    tt.write(source_reg_crop);
    tt.writeDirectory()
    close(tt)
    
    tt2.setTag(tagstruct);
    tt2.write(source_c2_reg_crop);
    tt2.writeDirectory()
    close(tt2)
    
    clear source_c1c2_reg_crop
    source_c1c2_reg_crop(:,:) = source_reg_crop;
    source_c1c2_reg_crop(:,:,2) = source_c2_reg_crop;
    source_c1c2_reg_crop = int16(source_c1c2_reg_crop);
    for ii = 1:size(source_c1c2_reg_crop,3)
        tt4.setTag(tagstruct);
        tt4.write(source_c1c2_reg_crop(:,:,ii));
        tt4.writeDirectory()
    end
    close(tt4)
    
    
    %% apply x and y shifts to the stacks
    if regStack
        for i = 1:num_images
            source(:,:,i) = imread([path source_filename], i ,'Info', info);
            source_c2(:,:,i) = imread([path_2 source_filename_c2], i ,'Info', info_2);
        end
        
        source_reg_stack(:,:,num_images) = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
        source_c2_reg_stack(:,:,num_images) = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift);
        for i = 1:num_images
            source_reg_stack(max_shift+1+outputs(b(1),3):size(source,1)+max_shift+1+outputs(b(1),3)-1,...
                max_shift+1+outputs(b(1),4):size(source,1)+max_shift+1+outputs(b(1),4)-1,...
                i) = source(:,:,i);
            source_c2_reg_stack(max_shift+1+outputs(b(1),3):size(source,1)+max_shift+1+outputs(b(1),3)-1,...
                max_shift+1+outputs(b(1),4):size(source,1)+max_shift+1+outputs(b(1),4)-1,...
                i) = source_c2(:,:,i);
        end
        source_reg_crop = source_reg_stack(max_shift:max_shift+size(source,1)-1,...
            max_shift:max_shift+size(source,2)-1,...
            :);
        source_c2_reg_crop = source_c2_reg_stack(max_shift:max_shift+size(source,1)-1,...
            max_shift:max_shift+size(source,2)-1,...
            :);
        
        % for display purposes
        a = min(min(min(source_reg_crop))) ;
        isnan(source_reg_crop);
        source_reg_crop(ans) = a-10;
        
        a = min(min(min(source_c2_reg_crop))) ;
        isnan(source_c2_reg_crop);
        source_c2_reg_crop(ans) = a-10;
        
        % save the images
        source_reg_crop = int16(source_reg_crop);
        source_c2_reg_crop = int16(source_c2_reg_crop);
        
        filename_c1 = [path source_filename(1:end-4) '_StackReg.tif'];
        filename_c2 = [path_2 source_filename_c2(1:end-4) '_StackReg.tif'];
        
        tt=Tiff(filename_c1,'w');
        tagstruct.ImageLength = size(source_reg_crop,1);
        tagstruct.ImageWidth = size(source_reg_crop,2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 16;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SampleFormat=Tiff.SampleFormat.Int;
        
        tt2=Tiff(filename_c2,'w');
        
        for ii = 1:size(source_reg_crop,3)
            tt.setTag(tagstruct);
            tt.write(source_reg_crop(:,:,ii));
            tt.writeDirectory()
        end
        close(tt)
        
        for ii = 1:size(source_c2_reg_crop,3)
            tt2.setTag(tagstruct);
            tt2.write(source_c2_reg_crop(:,:,ii));
            tt2.writeDirectory()
        end
        close(tt2)
    end
end
try
data.s2p{j}.rescale_factor = rescale_factor;
catch
    data.s2p.rescale_factor = rescale_factor; % for "SortData_PhotoStim_CM_v3.m" version older than 16 June 2020
end
%% Plot ROI Masks generated by s2p
PlotROIsMasks_2(data,rescale_factor,path)
clear outputs
end % end the loop through planes

%% save rescale factor to data
if exist(data.Info.mainDir,'dir') == 7 && isfield(data.Info,'saveName')
    mainDir = data.Info.mainDir;
    saveName = data.Info.saveName;
elseif isfield(data.Info,'saveName')
    mainDir = uigetdir('E:\TempData\','select dir to save data');
    saveName = data.Info.saveName;
else
    mainDir = uigetdir('E:\TempData\','select dir to save data');
    answer = inputdlg('save name (animal_position)');
    saveName = answer{1,1};
end

save([mainDir, '\Deconv\analysis\Spks_', saveName],'data','-append')    
end