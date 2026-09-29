%%%%%%%%   For selection (selecting based on iscell if selection = 1) %%%%%%%%%%%%%%%%%%%
% 1) Pre-process the data
%      - Subtract neuropil or not
%      - Do dF/F. F0 is trial baseline
%      - Sort data acording to trial type
%      - Only gonna work for 2 sets of {spatial fq,temporal fq} (because
%      it's used as an logical to sort data, see trialsID structure)
% 2) Plot a bunch of graph to get an idea of the quality of the data

% trialsID structure, logical in columns are:
% 1: no visual stim
% 2-i: direction 1-(i-1) (usually 8)
% nDir+1, nDir+2: {sf,tf} 1,2 (sorted based on sf. to keep in mind if tf and sf gets decoupled)
% nDir+3, nDir+4: odd, even trials (i.e, no photostim, photostim)

% Modified to support multiple planes, 16 June 2020

function data = SortData_RandPhotoStim_CM_v3_Deconv(Deconvtraces_all, session_data, selection)
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neuropil_sub = 1;
neuropil_factor = 0.7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
% - session data info
baseline = session_data.TrialSettings.GUI1_LED_CM.baseline;
trial_length = session_data.TrialSettings.GUI1_LED_CM.trialLength;
visStimLength = session_data.TrialSettings.GUI1_LED_CM.stimLength;
xstim_vis = [0 visStimLength];
xstim_opto = [0 1];
x = [-baseline trial_length-baseline-0.2];
nRep = session_data.TrialSettings.GUI1_LED_CM.nReps;
nTrials = session_data.TrialSettings.GUI2_LED_CM.nTrials;
nTrials_per2pFile = session_data.TrialSettings.GUI1_LED_CM.nStim2PTrig;

directions = session_data.TrialSettings.GUI1_LED_CMa.stimDir;
sfValues = session_data.TrialSettings.GUI1_LED_CMa.spatFreq;
tfValues = session_data.TrialSettings.GUI1_LED_CMa.tempFreq;

if nTrials == size(session_data.TrialSettings.GUI2_LED_CMa.dirTot,2)
    disp('session_data has the expected size')
    dirTot = session_data.TrialSettings.GUI2_LED_CMa.dirTot;
    sfTot = session_data.TrialSettings.GUI2_LED_CMa.sfTot;
    LaserOnTrials = session_data.TrialSettings.GUI2_LED_CMa.LaserOnTrials;
    
else
    disp(['session_data size does not match, cutting it to ',num2str(nTrials)])
    dirTot = session_data.TrialSettings.GUI2_LED_CMa.dirTot(1,1:nTrials);
    sfTot = session_data.TrialSettings.GUI2_LED_CMa.sfTot(1,1:nTrials);
end

% - Other general session parameters
n2PTrials = length(Deconvtraces_all{1}.header.numberOfFrames);
nFrames_min = min(double(Deconvtraces_all{1}.header.numberOfFrames));
disp(['min number of frames per 2P file: ' num2str(nFrames_min)])
nFrames_min = nFrames_min - mod(nFrames_min,4);
% - Converts frames into time
timeVectTrial = 1:nFrames_min;
try
    dtCa = 1/(Deconvtraces_all{1}.header.frameFrequency/Deconvtraces_all{1}.H.SI.hFastZ.numFramesPerVolume); % legacy (before June 2020)
    disp('multiple planes')
catch
    dtCa = 1/(Deconvtraces_all{1}.header.frameFrequency);
end
timeVect = timeVectTrial*dtCa;
timeVect = timeVect - baseline;
timeBaseInd = timeVect<xstim_vis(1);
timeStimInd = timeVect>xstim_vis(1) & timeVect<xstim_vis(2);
xvect = timeVect>x(1)& timeVect<x(2);
xAxis = find(xvect);
disp(['Baseline, frames: ', num2str(find(timeBaseInd,1)), ' - ', num2str(find(timeBaseInd,1,'last'))])
disp(['Stim, frames: ', num2str(find(timeStimInd,1)), ' - ', num2str(find(timeStimInd,1,'last'))])

% - for saving
try
    mainDir = Deconvtraces_all{1}.Info.mainDir;
    if ~exist(mainDir,'Dir')
        mainDir = uigetdir('X:\camille.mazo\2P_processed','Update main dir');
    end
catch
    mainDir = uigetdir('E:\TempData','Choose main dir');
end
saveDir = [mainDir, '\Deconv\Figures\QualityCheck'];
saveDirROI = [mainDir, '\Deconv\Figures\ROIs'];
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
try
    saveName = Deconvtraces_all{1}.Info.saveName;
catch
    answer = inputdlg('savename');
    saveName = answer{1,1};
end

disp(mainDir);disp(saveName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process the data for each of the imaging plane
dFF_pl = cell(1,length(Deconvtraces_all)); ROIs_pl = cell(1,length(Deconvtraces_all));
for j = 1:length(Deconvtraces_all)
    disp(['plane ', num2str(j)])
    deconvtraces = Deconvtraces_all{j};
    if selection
        ROIs_pl{j} = logical(deconvtraces.s2p.iscell);
    else
        ROIs_pl{j} = true(size(deconvtraces.s2p.iscell,1));
    end
    nROIs(j) = size(deconvtraces.spks,1);
    
    %% Check for bleaching
    data_GrandMean = mean(deconvtraces.spks,1);
    figure; plot(data_GrandMean); title({'Session Grand Mean'})
    xlabel('Frame #'); ylabel('spikes');
    saveas(gcf,[saveDir,'\Bleaching.tif'])
    
    
    %% Chop time-serie into actual trials. Do Normalisation
    counter1 = 1;
    counter2 = 1;
    for i = 1:n2PTrials
        data_stacked_frameNorm(:,counter2:counter2+nFrames_min-1) = deconvtraces.spks(:,counter1:counter1+nFrames_min-1);
        data_stacked_frameNormRawF(:,counter2:counter2+nFrames_min-1) = deconvtraces.RawF(:,counter1:counter1+nFrames_min-1);
        counter1 = counter1 + deconvtraces.header.numberOfFrames(i);
        counter2 = counter2 + nFrames_min;
    end
    
    %% Remove negative values
    min_RawF=min(min(data_stacked_frameNormRawF));
    data_stacked_frameNormRawF=data_stacked_frameNormRawF- min_RawF;
    
    %% Normalisation (F0 Calculation on Raw F traces)
    MovmeanFrameWindow=60; % number of baselines for moving average - if frame rate
    
    data_RawFLaserTrials=reshape(data_stacked_frameNormRawF,nROIs(j),(nFrames_min/nTrials_per2pFile),nTrials);
    F0_temp=data_RawFLaserTrials(:,timeBaseInd,:);
    
    F0_temp2=mean(F0_temp,2);
    F0_temp3=movmean(F0_temp2,MovmeanFrameWindow,3);
    F0=repmat(F0_temp3,1,nFrames_min/nTrials_per2pFile,1);
    
    F0NotMovmean=repmat(F0_temp2,1,nFrames_min/nTrials_per2pFile,1); % just for plotting
    
    data_stacked=reshape(data_stacked_frameNorm,nROIs(j),nFrames_min/nTrials_per2pFile,nTrials);
    data_spks_pl{j} = data_stacked./F0;
    
    clear data_stacked_frameNorm
    
    data.s2p{j} = deconvtraces.s2p;
    
end %end the loop through the different imaging planes

%% Concatenate data from all the planes
data_spks = cat(1,data_spks_pl{:});
ROIs = cat(1,ROIs_pl{:});

%% segregate photostim trials from the controls - this is just for plotting basic graph
for i = 1:size(data_spks,3)
    if LaserOnTrials(i)==2
        data_photostim(:,:,i) = data_spks(ROIs,:,i);
    else
        data_ctrl(:,:,i) = data_spks(ROIs,:,i);
    end
end

Popdata_photostim = squeeze(nanmean(data_photostim,1));
Grandmean_photostim = mean(Popdata_photostim,2);
Grandmean_photostim_sem = std(Popdata_photostim,0,2)/sqrt(size(Popdata_photostim,2));

Popdata_ctrl = squeeze(nanmean(data_ctrl,1));
Grandmean_ctrl = mean(Popdata_ctrl,2);
Grandmean_ctrl_sem = std(Popdata_ctrl,0,2)/sqrt(size(Popdata_ctrl,2));

Grand_diff = Grandmean_photostim - Grandmean_ctrl;

figure;
subplot(1,2,1); hold on
plot(timeVect(xvect),Grandmean_photostim(xvect),'r');
plot(timeVect(xvect),Grandmean_ctrl(xvect),'b')
plot(timeVect(xvect),Grand_diff(xvect),'k')
plot(timeVect(xvect),Grandmean_photostim(xvect)+Grandmean_photostim_sem(xvect),'r:');
plot(timeVect(xvect),Grandmean_photostim(xvect)-Grandmean_photostim_sem(xvect),'r:');
plot(timeVect(xvect),Grandmean_ctrl(xvect)+Grandmean_ctrl_sem(xvect),'b:');
plot(timeVect(xvect),Grandmean_ctrl(xvect)-Grandmean_ctrl_sem(xvect),'b:');

plot([xstim_opto(1) xstim_opto(1)],ylim,'k--')
plot([xstim_opto(2) xstim_opto(2)],ylim,'k--')
title('ROIs Grand Average, all trials mixed')
xlabel('Time (s)'); ylabel('Spikes');
legend({'photostim','ctrl','photostim-ctrl'},'Location','Northwest')

%% Prepare vectors to sort out the traces based on trial identity
DirTypes = unique(dirTot);
nDir = size(unique(dirTot),2);
nFQ = size(unique(sfTot),2)-1;
nVisStimTypes = (size(unique(dirTot),2)-1)*(size(unique(sfTot),2)-1); % -1 to remove noVis trial
DirStimTypes_sorting = sort(unique(dirTot),2);

trialsID = NaN(nTrials,nDir+2);
for ii = 1:nDir
    trialsID(:,ii) = (dirTot(1:nTrials) == DirStimTypes_sorting(ii));
end
trialsID(:,nDir+1) = logical(sfTot==1); % {sf.tf} 1
trialsID(:,nDir+2) = logical(sfTot==2); % {sf.tf} 2
trialsID(:,nDir+3) = logical(LaserOnTrials==1); % no photostim trials
trialsID(:,nDir+4) = logical(LaserOnTrials==2); % photostim trials


%  -- check if noVis trials exist and if the visStim trials are balanced
if DirTypes(1) == 0
    NoVisStim_exist = 1;
    disp('There are trials with no vis stim')
else
    NoVisStim_exist = 0;
    disp('No trials with no vis stim!!! Maybe use v2 of the code')
end
nRep = session_data.TrialSettings.GUI1_LED_CM.nReps;
nRep_NoVis(1) = sum(trialsID(:,1) & trialsID(:,nDir+4));
nRep_NoVis(2) = sum(trialsID(:,1) & ~trialsID(:,nDir+4));
for i = 2:nDir
    nRep_Vis(1,i-1) = sum(trialsID(:,i) & trialsID(:,nDir+1) & trialsID(:,nDir+4));
    nRep_Vis(2,i-1) = sum(trialsID(:,i) & ~trialsID(:,nDir+1) & trialsID(:,nDir+4));
    nRep_Vis(3,i-1) = sum(trialsID(:,i) & trialsID(:,nDir+1) & ~trialsID(:,nDir+4));
    nRep_Vis(4,i-1) = sum(trialsID(:,i) & ~trialsID(:,nDir+1) & ~trialsID(:,nDir+4));
end
if any(nRep_NoVis~=nRep)
    disp('Non-Visual trials are balanced')
end
if any(nRep_Vis~=nRep)
    disp('Visual trials are balanced')
end
if nRep_NoVis(1) == nRep_Vis(1,1)
    disp('same number of visual and non-visual trials')
end

%% now actually sort the data
% -- dir 1 is not a direction, it's the no vis stim trials. therefore it does
% not come with 2 SF/TFs. --
data_sorted(:,:,1,1,1,:) = data_spks(ROIs,:,trialsID(:,1) & trialsID(:,nDir+3));
data_sorted(:,:,1,1,2,:) = data_spks(ROIs,:,trialsID(:,1) & trialsID(:,nDir+4));
if nFQ == 2
    data_sorted(:,:,1,2,1,:) = NaN(size(data_sorted,1),size(data_spks,2),nRep);
    data_sorted(:,:,1,2,2,:) = NaN(size(data_sorted,1),size(data_spks,2),nRep);
end
for i = 2:nDir
    for ii = 1:nFQ
        for iii = 1:2
            data_sorted(:,:,i,ii,iii,:) = data_spks(ROIs,:,trialsID(:,i) & trialsID(:,nDir+ii) & trialsID(:,nDir+2+iii));
        end
    end
end
sem_laserOFF = std(nanmean(data_sorted(:,xvect,1,1,1,:),6),[],1)./sqrt(size(data_sorted,1));
sem_laserON = std(nanmean(data_sorted(:,xvect,1,1,2,:),6),[],1)./sqrt(size(data_sorted,1));
subplot(1,2,2); hold on
p1 = plot(timeVect(xvect),mean(mean(data_sorted(:,xvect,1,1,1,:),6),1),'k');
p2 = plot(timeVect(xvect),mean(mean(data_sorted(:,xvect,1,1,2,:),1),6),'r');
plot(timeVect(xvect),mean(mean(data_sorted(:,xvect,1,1,1,:),6),1)-sem_laserOFF,'k:')
plot(timeVect(xvect),mean(mean(data_sorted(:,xvect,1,1,1,:),6),1)+sem_laserOFF,'k:')
plot(timeVect(xvect),mean(mean(data_sorted(:,xvect,1,1,2,:),6),1)-sem_laserON,'r:')
plot(timeVect(xvect),mean(mean(data_sorted(:,xvect,1,1,2,:),6),1)+sem_laserON,'r:')
plot([xstim_opto(1) xstim_opto(1)],ylim,'k--')
plot([xstim_opto(2) xstim_opto(2)],ylim,'k--')
xlabel('Time (s)'); ylabel('dF/F');
title('Opto Only, ROIs Grand Average')
legend([p1,p2],{'ctrl','photostim'},'Location','Northwest')
set(gcf,'Units','Normalized','Position',[0.1 0.5 0.8 0.4])
saveas(gcf,[saveDir,'\GrandMeanTraces.tif'])

%% save
data.data_sorted = data_sorted;

data.Info.baseline = baseline;
data.Info.trial_length = trial_length;
data.Info.visStimLength = visStimLength;
data.Info.timeVect = timeVect;
data.Info.xvect = xvect;
data.Info.timeStimInd = timeStimInd;
data.Info.timeBaseInd = timeBaseInd;
data.Info.xstim_vis = xstim_vis;
data.Info.xstim_opto = xstim_opto;
data.Info.x = x;
data.Info.nROIs = nROIs;
data.Info.nRep = nRep;
data.Info.nTrials = nTrials;
data.Info.nVisStimTypes = nVisStimTypes;
data.Info.directions = directions;
data.Info.sfValues = sfValues;
data.Info.tfValues = tfValues;
data.Info.LaserOnTrials = LaserOnTrials;
data.Info.NoVisStim_exist = NoVisStim_exist;
data.Info.saveName = saveName;
data.Info.mainDir = mainDir;
data.Info.neuropil_sub    = neuropil_sub;
data.Info.neuropil_factor = neuropil_factor;
data.Info.trialsID = trialsID;
data.Info.nPlanes = length(Deconvtraces_all);

save([mainDir,'\Deconv\analysis\Spks_', saveName],'data')
disp('data saved');


%% ... and do some more basic plots
trial_mean = mean(data_sorted(:,:,:,:,:,:),6);
concat = reshape(trial_mean,size(trial_mean,1)*size(trial_mean,2)*nDir*nFQ*2,1);
ymax = max(concat);
ymin = min(concat);
figure;
for i = 1:nDir-1
    for ii = 1:nFQ
        for iii = 1:2
            subtightplot(2,(nDir-1)*nFQ,i+(nDir-1)*(ii-1)+(nDir-1)*nFQ*(iii-1)); hold on
            imagesc(trial_mean(:,:,i+1,ii,iii),[ymin ymax]);
            colormap('gray')
            if i == 1 && ii == 1 && iii == 1
                ylabel('Vis Only')
            elseif i == 1 && ii == 1 && iii == 2
                ylabel('Vis + Opto')
                xlabel('Frame #')
            else
                axis off
            end
        end
    end
end
c = colorbar;
c.Position = [0.96 0.055 0.01 0.43];
c.Label.String = 'dF/F';
sgtitle({'Impact of laser on drifting gratings'})
set(gcf,'units','normalized','Position',[0.01 0.3 0.98 0.6])
saveas(gcf,[saveDir,'\PhotoStim_VisStim.tif'])

end