
global MainDir

%% Preambule

alphaVal     = 0.01;
useLoop = 1;
PowerDatasetInc=1;
Baseline=[1:15];
StimONResponse=[16:30];
ThModulation = alphaVal*5;
PMcoord=1;
DistInj=1;

if PowerDatasetInc
    if PMcoord
        load('E:\2P\Pop Analysis\Final_datasets\LoopV1PM_powers_deconv_final\PMcoordinates_2026_02_25.mat')
    end
    if DistInj
        load('E:\2P\Pop Analysis\Final_datasets\LoopV1PM_powers_deconv_final\DistInjSiteLP_2026_01_16.mat')
    end
else
    if PMcoord
        load('E:\2P\Pop Analysis\Final_datasets\LoopV1PM_deconv_final\PMcoordinates_2025_05_14.mat')
    end
    if DistInj
        load('E:\2P\Pop Analysis\Final_datasets\LoopV1PM_deconv_final\DistInjSite_2025_06_02.mat')
    end
end


tv = datestr(now, 'yyyy_mm_dd');
fprintf('\n')

MainDir = uigetdir('E:\2P\Pop Analysis','Master directory');

saveDir_plots = [MainDir filesep 'plot'];

if ~exist(saveDir_plots,'dir')
    mkdir(saveDir_plots)
end

oldfolder = cd([MainDir]);
cd(oldfolder)

%% Load data
FFFB_2_CellSelec_LoadData

%% Select responsive cells and order Directions

%% Preambule

n_animals_FB= size(data{1},1);n_animals_FF= size(data{2},1);
for i=1:n_animals_FB
    n_pos_FB(i)=sum(~cellfun('isempty',data{1}(i,:)));
end
for i=1:n_animals_FF
    n_pos_FF(i)=sum(~cellfun('isempty',data{2}(i,:)));
end


%Looped cells
if useLoop
    for i=1:n_animals_FB
        for ii=1:n_pos_FB(i)
            LoopSes{1}{i,ii}=info{1}{i,ii}.looped;
        end
    end
    
    for i=1:n_animals_FF
        for ii=1:n_pos_FF(i)
            LoopSes{2}{i,ii}=info{2}{i,ii}.looped;
        end
    end
end

% N positions and animals
n_pos{1}= n_pos_FB; n_pos{2}= n_pos_FF;
n_ani{1}= n_animals_FB; n_ani{2}= n_animals_FF;

%%
for F=1:2
    
    %% for PM, select the cells that are within PM
    if PMcoord
        FFFB_SelPM_Final
    end
    %% for V1, select the cells that are far from Injection Site
    if DistInj
        FFFB_mapROIdistToInjSite
    end
    %% select only the PM cells for FF and the cells far from PM for FB in data and Looped cells
    if PMcoord & DistInj
        for ani=1:n_ani{F}
            for pos=1:n_pos{F}(ani)
                dataSel{F}{ani,pos}=data{F}{ani,pos}(InPM{F}{ani,pos}&SelDistInjSite{F}{ani,pos},:,:,:,:,:);
                ROICoordyx{F}{ani,pos}=ROICoordyx{F}{ani,pos}(InPM{F}{ani,pos}&SelDistInjSite{F}{ani,pos},:);
                if useLoop
                    LoopSes{F}{ani,pos}=LoopSes{F}{ani,pos}(InPM{F}{ani,pos}&SelDistInjSite{F}{ani,pos});
                end
            end
        end
    end

    %% Baseline substraction
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            BaselineAv{F}{ani,pos}=mean(dataSel{F}{ani,pos}(:,Baseline,:,:,:,:),2);
            dataSel{F}{ani,pos}=dataSel{F}{ani,pos}-BaselineAv{F}{ani,pos};
        end
    end
    
    
%% Subselect data that matches the loop dataset 
% (power dataset has more laser power condtions)
if PowerDatasetInc
for ani=1:n_ani{F}
    for pos=1:n_pos{F}(ani)
        if size(dataSel{F}{ani,pos},5)>2
data_temp{F}{ani,pos}=dataSel{F}{ani,pos};
dataSel{F}{ani,pos}=NaN(size(data_temp{F}{ani,pos},1),size(dataSel{F}{ani,pos},2),9,2,2,20);
dataSel{F}{ani,pos}(:,:,[1 2 4 6 8],:,1,1:15)=data_temp{F}{ani,pos}(:,:,:,:,1,:);
dataSel{F}{ani,pos}(:,:,[1 2 4 6 8],:,2,1:15)=data_temp{F}{ani,pos}(:,:,:,:,5,:);
        end
    end
end
end
    
end
clear data
%% Structure data
for F=1:2
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)

            % All dirs cat
            dataAllSesStimCat_temp{F}{ani,pos}(:,:,:,:,1:size(dataSel{F}{ani,pos},6))= dataSel{F}{ani,pos}(:,:,2:9,1,:,:);
            dataAllSesStimCat_temp{F}{ani,pos}(:,:,:,:,size(dataSel{F}{ani,pos},6)+1:size(dataSel{F}{ani,pos},6)*2)= dataSel{F}{ani,pos}(:,:,2:9,2,:,:); % cat the 2sftf
            dataAllSesStimCat_temp2{F}{ani,pos}=permute(dataAllSesStimCat_temp{F}{ani,pos},[1 2 4 3 5]);
            dataAllSesStimCat{F}{ani,pos}=reshape(dataAllSesStimCat_temp2{F}{ani,pos}, size(dataAllSesStimCat_temp2{F}{ani,pos},1),size(dataAllSesStimCat_temp2{F}{ani,pos},2),size(dataAllSesStimCat_temp2{F}{ani,pos},3),size(dataAllSesStimCat_temp2{F}{ani,pos},4)*size(dataAllSesStimCat_temp2{F}{ani,pos},5)); clear dataAllSesStimCat_temp
            dataAllSesStimCatAll{F} = cat(1,dataAllSesStimCat{F}{:});
                        
        end
    end
    
    %% Select responsive cells (All Stimuli)
    
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            VisTrials{F}{ani,pos}=dataSel{F}{ani,pos}(:,:,2:9,:,:,:);
            OptoCatAS{F}{ani,pos}= reshape(VisTrials{F}{ani,pos},size(VisTrials{F}{ani,pos},1),size(VisTrials{F}{ani,pos},2),size(VisTrials{F}{ani,pos},6)*2*2*8);
            
            Base=squeeze(nanmean(OptoCatAS{F}{ani,pos}(:,Baseline,:),2));
            StimON=squeeze(nanmean(OptoCatAS{F}{ani,pos}(:,StimONResponse,:),2));
            BaseRepAv=squeeze(nanmean(Base,2));
            StimONRepAv=squeeze(nanmean(StimON,2));
            responsiveStimONAllDir{F}{ani,pos}=false(size(VisTrials{F}{ani,pos},1),1);
            hStimON=0;
            
            for ROI=1:size(VisTrials{F}{ani,pos},1)
                hStimON=ttest(squeeze(Base(ROI,:)),squeeze(StimON(ROI,:)),alphaVal);
                if hStimON & StimONRepAv(ROI)>BaseRepAv(ROI)
                    responsiveStimONAllDir{F}{ani,pos}(ROI)=hStimON;
                end
            end
            
            
        end
    end
    
end

%% Select Best SFTF
for F=1:2
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if ~isempty(dataSel{F}{ani,pos})
                VisTrialsMeanTimeTrials{F}{ani,pos}=squeeze(nanmean(nanmean(VisTrials{F}{ani,pos}(:,StimONResponse,:,:,:,:),2),6));
                
                for ROI=1:size(VisTrials{F}{ani,pos},1)
                    for SF=1:2
                        [M(SF),~]=max(VisTrialsMeanTimeTrials{F}{ani,pos}(ROI,:,SF,1));
                    end
                    if M(1)>M(2)
                        VisBestSF{F}{ani,pos}(ROI,:,:,:,:)=VisTrials{F}{ani,pos}(ROI,:,:,1,:,:);
                    elseif M(2)>M(1)
                        VisBestSF{F}{ani,pos}(ROI,:,:,:,:)=VisTrials{F}{ani,pos}(ROI,:,:,2,:,:);
                    end
                end
                
            end
        end
    end
end

%% Responsive Best SF                
for F=1:2
    
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if ~isempty(dataSel{F}{ani,pos})
temp=permute(VisBestSF{F}{ani,pos},[1 2 4 5 3]);
                VisTrialsBestSFAD{F}{ani,pos}=reshape(temp,size(temp,1),size(temp,2),size(temp,3),[]);

        MatOptoCat=[];
        MatOptoCat= reshape(VisTrialsBestSFAD{F}{ani,pos},size(VisTrialsBestSFAD{F}{ani,pos},1),size(VisTrialsBestSFAD{F}{ani,pos},2),[]);
        
        Base=squeeze(nanmean(MatOptoCat(:,Baseline,:),2));
        StimON=squeeze(nanmean(MatOptoCat(:,StimONResponse,:),2));
        BaseRepAv=squeeze(nanmean(Base,2));
        StimONRepAv=squeeze(nanmean(StimON,2));
        responsiveADBestSF{F}{ani,pos}=false(size(MatOptoCat,1),1);
        hStimON=0;
        
        for ROI=1:size( MatOptoCat,1)
            hStimON=ttest(squeeze(Base(ROI,:)),squeeze(StimON(ROI,:)),alphaVal);
            if hStimON & StimONRepAv(ROI)>BaseRepAv(ROI)
                responsiveADBestSF{F}{ani,pos}(ROI)=hStimON;                
            end
        end
            end
        end
    end
end
%% Excited and inhibited neurons from responsive neurons
for F=1:2
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            disp(['Animal ' num2str(ani) ' | Pos ' num2str(pos)])
            
            AllDirRespAllRep{F}{ani,pos}=squeeze(nanmean(dataAllSesStimCat{F}{ani,pos}(:,StimONResponse,:,:),2));
            AllDirResp{F}{ani,pos}=squeeze(nanmean(nanmean(dataAllSesStimCat{F}{ani,pos}(:,StimONResponse,:,:),2),4));
            
            nResponsive{F}{ani,pos} = sum(responsiveStimONAllDir{F}{ani,pos});
            
            posi{F}{ani,pos} = zeros(1,size(AllDirRespAllRep{F}{ani,pos},1)); neg{F}{ani,pos} = zeros(1,size(AllDirRespAllRep{F}{ani,pos},1)); nochange{F}{ani,pos} = zeros(1,size(AllDirRespAllRep{F}{ani,pos},1));
            
            for ROI = 1 : size(dataAllSesStimCat{F}{ani,pos},1)
                if responsiveStimONAllDir{F}{ani,pos}(ROI)
                    [h2, p2] = ttest2(squeeze(AllDirRespAllRep{F}{ani,pos}(ROI,1,:)),squeeze(AllDirRespAllRep{F}{ani,pos}(ROI,2,:)),'Alpha',alphaVal);
                    
                    if h2 && AllDirResp{F}{ani,pos}(ROI,2) > AllDirResp{F}{ani,pos}(ROI,1)
                        posi{F}{ani,pos}(ROI)= true;
                    elseif  h2 && AllDirResp{F}{ani,pos}(ROI,2) < AllDirResp{F}{ani,pos}(ROI,1)
                        neg{F}{ani,pos}(ROI) = true;
                    elseif h2 == 0
                        nochange{F}{ani,pos}(ROI) = true;
                    else
                        debug(ROI) = true;
                        disp('bug. this cannot be')
                    end
                end
            end
            
            Frac{F}{ani,pos} = [sum(posi{F}{ani,pos},'omitnan') sum(neg{F}{ani,pos},'omitnan') sum(nochange{F}{ani,pos},'omitnan')]./nResponsive{F}{ani,pos};
            
            SigExcAllDirs{F}{ani,pos}=posi{F}{ani,pos}';
            SigInhAllDirs{F}{ani,pos}=neg{F}{ani,pos}';
            
            
            
        end
    end
    
    
    %% Selection of sessions that are modulated
    
    SelModulation{F}=false(size(SigExcAllDirs{F}));
    
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if Frac{F}{ani,pos}(1)+Frac{F}{ani,pos}(2)>ThModulation
                SelModulation{F}(ani,pos)=true;
            end
        end
    end
end

%% Modulated cells photostim alone
for F=1:2
    
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if ~isempty(dataSel{F}{ani,pos})
                LaserAloneTrials{F}{ani,pos}=squeeze(dataSel{F}{ani,pos}(:,:,1,1,:,:));
                LaserAloneAllRep{F}{ani,pos}=squeeze(nanmean(LaserAloneTrials{F}{ani,pos}(:,StimONResponse,2,:),2));
                LaserAloneAllRepBase{F}{ani,pos}=squeeze(nanmean(LaserAloneTrials{F}{ani,pos}(:,Baseline,2,:),2));
                LaserAlone{F}{ani,pos}=squeeze(nanmean(nanmean(LaserAloneTrials{F}{ani,pos}(:,StimONResponse,2,:),2),4));
                LaserAloneBase{F}{ani,pos}=squeeze(nanmean(nanmean(LaserAloneTrials{F}{ani,pos}(:,Baseline,2,:),2),4));
                nROIs{F}{ani,pos} = size(LaserAloneTrials{F}{ani,pos},1);
                
                posi{F}{ani,pos} = zeros(1,size(LaserAloneAllRep{F}{ani,pos},1)); neg{F}{ani,pos} = zeros(1,size(LaserAloneAllRep{F}{ani,pos},1)); nochange{F}{ani,pos} = zeros(1,size(LaserAloneAllRep{F}{ani,pos},1));
                
                for ROI = 1 : size(LaserAloneAllRep{F}{ani,pos},1)
                    [h2, p2] = ttest(LaserAloneAllRep{F}{ani,pos}(ROI,:),LaserAloneAllRepBase{F}{ani,pos}(ROI,:),'Alpha',alphaVal);
                    if h2 && LaserAlone{F}{ani,pos}(ROI) > LaserAloneBase{F}{ani,pos}(ROI)
                        posi{F}{ani,pos}(ROI)= true;
                    elseif  h2 && LaserAlone{F}{ani,pos}(ROI) < LaserAloneBase{F}{ani,pos}(ROI)
                        neg{F}{ani,pos}(ROI) = true;
                    elseif h2 == 0
                        nochange{F}{ani,pos}(ROI) = true;
                    else
                        disp('bug. this cannot be')
                    end
                    
                end
                SigExcLaserAlone{F}{ani,pos}=posi{F}{ani,pos}';
                SigInhLaserAlone{F}{ani,pos}=neg{F}{ani,pos}';
            end
        end
    end
    
    
    
end
%% Conversion of xy from pixel to um
for F=1:2
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            ROICoordyx_um{F}{ani,pos}=0.5862 * ROICoordyx{F}{ani,pos};
        end
    end
    
    %% Structure data
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if SelModulation{F}(ani,pos)
                data{F}{ani,pos}.Raw=dataSel{F}{ani,pos};
                data{F}{ani,pos}.StimAv=dataAllSesStimCat{F}{ani,pos};
                data{F}{ani,pos}.VisBestSF=VisBestSF{F}{ani,pos};
                CellStat{F}{ani,pos}.responsiveAD=responsiveStimONAllDir{F}{ani,pos};
                CellStat{F}{ani,pos}.responsiveADBestSF=responsiveADBestSF{F}{ani,pos};
                CellStat{F}{ani,pos}.SigExcAD=SigExcAllDirs{F}{ani,pos};
                CellStat{F}{ani,pos}.SigInhAD=SigInhAllDirs{F}{ani,pos};
                CellStat{F}{ani,pos}.Loop=LoopSes{F}{ani,pos};
                CellStat{F}{ani,pos}.Coordyx_um=ROICoordyx_um{F}{ani,pos};
                CellStat{F}{ani,pos}.SigExcLaserAlone=SigExcLaserAlone{F}{ani,pos};
                CellStat{F}{ani,pos}.SigInhLaserAlone=SigInhLaserAlone{F}{ani,pos};
                %                 CellStat{F}{ani,pos}.DistInjSite=DistInjSite{F}{ani,pos}.distanceToInjSite;
            end
        end
        CellStat{F}{ani}.MouseID{F}=MouseID{F}{ani};
    end
    
    
end

disp(['FB Number of sessions passing modulation threshold: ' num2str(sum(sum(~cellfun(@isempty,data{1}))))])
disp(['FF Number of sessions passing modulation threshold: ' num2str(sum(sum(~cellfun(@isempty,data{2}))))])

Params.colors{1}=[0 0.8 0.6];
Params.colors{2}=[0.6 0.2 1];
Params.MainDir=MainDir;
Params.ExpName{1}='FB'; Params.ExpName{2}='FF';
Params.n_ani=n_ani;
Params.n_pos=n_pos;
Params.StimONResponse=StimONResponse;



VarDir=[MainDir filesep 'Variables'];
if ~exist(VarDir,'dir')
    mkdir(VarDir)
end

disp('Saving selection...')

save([VarDir filesep 'Selec_' tv],'data','CellStat','Params' );

