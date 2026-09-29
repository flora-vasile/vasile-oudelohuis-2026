function path = Alignment_Final(channel2)
threshold = 60;
[stack_fileName, path] = uigetfile({'E:\TempData\*.*'},'Load a z-stack');
[~,zStack]=opentif_SI21([path stack_fileName]);

if length(size(zStack))==3 % if the zStack has 3 dimensions (bug in scanimage), reshape it
    nFrames = input('number of frames per slice: ');
    nSlices=size(zStack,3)/(2*nFrames);
    zStack = reshape(zStack,size(zStack,1),size(zStack,2),2,nFrames,nSlices);
end

data_c2 = squeeze(zStack(:,:,2,:,:,:));
data_c1 = squeeze(zStack(:,:,1,:,:,:));

%% detect motion on plane average
data_c1_avg = squeeze(mean(data_c1,3));

for i = 1:size(data_c1_avg,3)
    c1_fft = fft(data_c1_avg(:,:,i));
    crispy(:,i,1) = abs(mean(c1_fft,1));
    crispy(:,i,2) = abs(mean(c1_fft,2));
end
crispy_avg = mean(crispy,3);

for  i =1:size(crispy_avg,2)
    AUC(i) = nansum(crispy_avg(350:end,i));
end
figure; hold on
histogram(AUC)
plot([threshold threshold],ylim,'r')
xlabel('Power of high frequencies')
title('Cut off for planes selection')

%% Max projection excluding bad planes
k = false(size(AUC,2),1);
k(AUC>threshold) = true;
figure;plot(k)
xlabel('Plane #')
title('Planes selected for Max Proj')

grandMean = mean(data_c1_avg(:,:,k),3);
for i = 1:size(data_c1_avg,2)
    for ii = 1:size(data_c1_avg,2)
        maxProj(i,ii) = max(data_c1_avg(i,ii,k));
    end
end

%% Register single frames to that Max Projection
fft_target = fft2(maxProj);
for i = 1:size(data_c1,3)
    for ii = 1:size(data_c1,4)
        source = data_c1(:,:,i,ii);
        fft_source = fft2(source);
        [output, ~] = dftregistration(fft_target,fft_source,1);
        outputs(i,ii,:) = output;
    end
end

max_shift = max(max(max(abs(outputs(:,:,3)))), max(max(abs(outputs(:,:,4)))));
if max_shift == 0
    source_reg_crop = data_c1;
    source_reg_crop_c2 = data_c2;
else
    source_reg = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift,size(source,3),11);
    source_reg_c2 = NaN(size(source,1)+2*max_shift,size(source,2)+2*max_shift,size(source,3),11);
    for i = 1:size(data_c1,4)
        for ii = 1:size(data_c1,3)
            source_reg(max_shift+1+outputs(ii,i,3):size(source,1)+max_shift+1+outputs(ii,i,3)-1,max_shift+1+outputs(ii,i,4):size(source,1)+max_shift+1+outputs(ii,i,4)-1,ii,i) = data_c1(:,:,ii,i);
            source_reg_crop(:,:,ii,i) = source_reg(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1,ii,i);
            if channel2
                source_reg_c2(max_shift+1+outputs(ii,i,3):size(source,1)+max_shift+1+outputs(ii,i,3)-1,max_shift+1+outputs(ii,i,4):size(source,1)+max_shift+1+outputs(ii,i,4)-1,ii,i) = data_c2(:,:,ii,i);
                source_reg_crop_c2(:,:,ii,i) = source_reg_c2(max_shift:max_shift+size(source,1)-1,max_shift:max_shift+size(source,2)-1,ii,i);
            end
        end
    end
end
%% save tiffs
data_c1_avg = int16(squeeze(nanmean(source_reg_crop,3)));
maxProj = int16(maxProj);

saveDir = [path 'Processed' filesep];
if ~ exist(saveDir)
    mkdir(saveDir)
end
filename_c1 = [saveDir stack_fileName(1:end-4) '_c1.tif'];
filename_maxProj = [saveDir stack_fileName(1:end-4) '_c1_MaxProj.tif'];

tt0 = Tiff(filename_maxProj,'w');
tt=Tiff(filename_c1,'w');
tagstruct.ImageLength = size(data_c1_avg,1);
tagstruct.ImageWidth = size(data_c1_avg,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.SampleFormat=Tiff.SampleFormat.Int;

tt0.setTag(tagstruct);
tt0.write(maxProj);
tt0.writeDirectory()
close(tt0)

for ii = 1:size(data_c1_avg,3)
    tt.setTag(tagstruct);
    tt.write(data_c1_avg(:,:,ii));
    tt.writeDirectory()
end
close(tt)

if channel2
    data_c2_avg = int16(squeeze(nanmean(source_reg_crop_c2,3)));
    filename_c2 = [saveDir stack_fileName(1:end-4) '_c2.tif'];
    tt2=Tiff(filename_c2,'w');
    for ii = 1:size(data_c2_avg,3)
        tt2.setTag(tagstruct);
        tt2.write(data_c2_avg(:,:,ii));
        tt2.writeDirectory()
    end
    close(tt2)
    
    data_combined(:,:,1,:) = data_c1_avg;
    data_combined(:,:,2,:) = data_c2_avg;
    filename_combined = [saveDir stack_fileName(1:end-4) '_c1c2.tif'];
    tt3=Tiff(filename_combined,'w');
    for ii = 1:size(data_combined,4)
        for i = 1:2
            tt3.setTag(tagstruct);
            tt3.write(data_combined(:,:,i,ii));
            tt3.writeDirectory()
        end
    end
    close(tt3)
end


end