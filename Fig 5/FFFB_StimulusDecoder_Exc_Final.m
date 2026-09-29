saveDir_Fig5 = [MainDir filesep 'plot\Fig 5'];
if ~exist(saveDir_Fig5,'dir')
    mkdir(saveDir_Fig5)
end

%% Cell selection
CellSubSelecName= ' SFTF pooled excited';
ThNbCells=10;
% Parameters for subsampling
nNeuronsSub = 10;    % number of neurons to subsample each iteration
nRepeats = 10;       % number of subsampling repetitions

%%
for F=1:2
    
    OFF_perf{F}=cell(n_ani{F},max(n_pos{F}));
    ON_perf{F}=cell(n_ani{F},max(n_pos{F}));
    
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if ~isempty(data{F}{ani,pos})
                CellSubSelec=logical(CellStat{F}{ani,pos}.SigExcAD);
                if sum(CellSubSelec)> ThNbCells
                    disp(['Mouse ' num2str(ani) ' pos ' num2str(pos)])
                    Vis=squeeze(nanmean(data{F}{ani,pos}.Raw(CellSubSelec,Params.StimONResponse,2:9,:,:,:),2));
                    
                    % concatenate trials from 2 SFTF
                    Vis2=permute(Vis,[1 2 4 5 3]);
                    MeanTime=reshape(Vis2, size(Vis2,1),size(Vis2,2),size(Vis2,3),[]);
                    
                    
                    % Assume: data is nNeurons × 8 × 2 × nReps
                    [nNeurons, nDir, nLaser, nReps] = size(MeanTime);
                    
                    % --- Split by laser condition ---
                    data_off{F}{ani,pos} = squeeze(MeanTime(:,:,1,:));   % [nNeurons × 8 × nReps]
                    data_on{F}{ani,pos}  = squeeze(MeanTime(:,:,2,:));   % [nNeurons × 8 × nReps]
                    
                    % --- Reshape: collapse directions and repetitions into a single dimension ---
                    data_off{F}{ani,pos} = reshape(permute(data_off{F}{ani,pos}, [1 3 2]), nNeurons, nDir*nReps);
                    data_on{F}{ani,pos}  = reshape(permute(data_on{F}{ani,pos}, [1 3 2]), nNeurons, nDir*nReps);
                    
                    % --- Create direction identity vector ---
                    dirIdentity = repelem(1:nDir, nReps);  % length = nDir*nReps
                    
                    %% Initialize arrays to store results
                    OFF_perf_rep = nan(nRepeats,1);
                    ON_perf_rep  = nan(nRepeats,1);
                    
                    for rep = 1:nRepeats
                        % Randomly sample neurons (without replacement)
                        if nNeuronsSub > nNeurons
                            warning('Requested more neurons than available — using all neurons.');
                            idxNeurons = 1:nNeurons;
                        else
                            idxNeurons = randperm(nNeurons, nNeuronsSub);
                        end
                        
                        %% --- Decoder Laser OFF ---
                        X = data_off{F}{ani,pos}(idxNeurons, :)';   % Trials x Neurons
                        y = dirIdentity(:);
                        
                        t = templateSVM('KernelFunction', 'linear', 'Standardize', true);
                        Mdl = fitcecoc(X, y, 'Learners', t, 'KFold', 5);
                        cvLoss = kfoldLoss(Mdl);
                        OFF_perf_rep(rep) = (1 - cvLoss) * 100;
                        
                        %% --- Decoder Laser ON ---
                        X = data_on{F}{ani,pos}(idxNeurons, :)';   % Trials x Neurons
                        y = dirIdentity(:);
                        
                        t = templateSVM('KernelFunction', 'linear', 'Standardize', true);
                        Mdl = fitcecoc(X, y, 'Learners', t, 'KFold', 5);
                        cvLoss = kfoldLoss(Mdl);
                        ON_perf_rep(rep) = (1 - cvLoss) * 100;
                    end
                    
                    % --- Store mean decoding accuracy across repetitions ---
                    OFF_perf{F}{ani,pos} = mean(OFF_perf_rep);
                    ON_perf{F}{ani,pos} = mean(ON_perf_rep);
                    
                end
            end
        end
    end
end

%% plot sessions
for F=1:2
    for ani=1:n_ani{F}
        tempOFF{F}{ani}=cat(1,cat(1,OFF_perf{F}{ani,:}));
        tempON{F}{ani}=cat(1,cat(1,ON_perf{F}{ani,:}));
    end
    OFF_perf_ses{F}=cat(1,tempOFF{F}{:});
    ON_perf_ses{F}=cat(1,tempON{F}{:});
end

% Plot normalised and with Anova for fig
plotDecoderComparison_mixedANOVA_fig(OFF_perf_ses, ON_perf_ses,'plotColors', [Params.colors{2}; Params.colors{1}],...
    'ExpName', [CellSubSelecName ' cells - ' num2str(nNeuronsSub) ' cells subsamp - Ses'],...
    'saveDir',saveDir_Fig5);

