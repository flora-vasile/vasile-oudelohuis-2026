
    
    FBdir = dir([MainDir filesep 'FB']);
    str = sortrows({FBdir.name}');
    [s,~] = listdlg('PromptString','Select animals to analyze:', 'OKString', 'OK',...
        'SelectionMode','multiple',...
        'ListString', str, 'Name', 'Select a File');
    FBanimals = str(s);
    
    FFdir = dir([MainDir filesep 'FF']);
    str = sortrows({FFdir.name}');
    [s,~] = listdlg('PromptString','Select animals to analyze:', 'OKString', 'OK',...
        'SelectionMode','multiple',...
        'ListString', str, 'Name', 'Select a File');
    FFanimals = str(s);
    
    n_animals_FB = length(FBanimals);
    n_animals_FF = length(FFanimals);
    
    clear s v str
    
    disp('FORCED to rebuild data structure')
    % - FB
    fprintf(2,'FB animals\n')
    for i = 1:n_animals_FB
        disp(FBanimals{i})
        pos_directories = dir([FBdir(1).folder filesep FBanimals{i}]);
        match = {}; n_pos_FB(i) = 0;
        for ii = 1:length(pos_directories)
            if ~isempty(regexp(pos_directories(ii).name,'pos','ONCE'))
                match{end+1} = pos_directories(ii).name;
                n_pos_FB(i) = n_pos_FB(i) +1;
                index(i,n_pos_FB(i)) = ii;
            end
        end
        for ii = 1:n_pos_FB(i)
            path = [pos_directories(index(i,ii)).folder filesep pos_directories(index(i,ii)).name];
            disp(pos_directories(index(i,ii)).name)
            file = dir(fullfile([path filesep], '*.mat'));
            temp=load([file.folder filesep file.name]);
            
            
            if size(temp.data.data_sorted,4) == 1 % Only 1 SF TF couple
                temp.data.data_sorted(:,:,:,2,:,:) = temp.data.data_sorted(:,:,:,1,:,:);
            end
            if size(temp.data.data_sorted,6) >20 % more than 20 Rep
                temp2 = temp.data.data_sorted(:,:,:,:,:,1:20);
                temp.data.data_sorted = temp2;
            end
            
            data_FB{i,ii} = temp.data.data_sorted;
            info_FB{i,ii} = temp.data.Info;
            info_FB{i,ii}.pos = file;
            info_FB{i,ii}.animal=i;
            info_FB{i,ii}.position=ii;
            
            %Looped cells
            if useLoop
                
                if ~isfield(temp.data,'beadID')
                    tempbeads=false(1,size(temp.data.data_sorted,1));
                    tempbeads(temp.data.beads_pos(:,3))=true;
                    temp.data.beadID.bead_pos=tempbeads;
                end
                info_FB{i,ii}.looped = temp.data.beadID.bead_pos.';
            end
            
            % XY coordinates
            
            for ROI=1:length(temp.data.s2p{1}.iscell)
                ROICoordyxtemp(ROI,:)   = temp.data.s2p{1}.stat{ROI}.med;
            end
            info_FB{i,ii}.ROICoordyx=ROICoordyxtemp(logical(temp.data.s2p{1}.iscell),:);
            
            
        end
    end
    clear temp
    % - FF
    fprintf(2,'FF animals\n')
    for i = 1:n_animals_FF
        disp(FFanimals{i})
        pos_directories = dir([FFdir(1).folder filesep FFanimals{i}]);
        match = {}; n_pos_FF(i) = 0;
        for ii = 1:length(pos_directories)
            if ~isempty(regexp(pos_directories(ii).name,'pos','ONCE'))
                match{end+1} = pos_directories(ii).name;
                n_pos_FF(i) = n_pos_FF(i) +1;
                index(i,n_pos_FF(i)) = ii;
            end
        end
        for ii = 1:n_pos_FF(i)
            path = [pos_directories(index(i,ii)).folder filesep pos_directories(index(i,ii)).name];
            disp(pos_directories(index(i,ii)).name)
            file = dir(fullfile([path filesep], '*.mat'));
            temp=load([file.folder filesep file.name]);
            
            if size(temp.data.data_sorted,4) == 1
                temp.data.data_sorted(:,:,:,2,:,:) = temp.data.data_sorted(:,:,:,1,:,:);
            end
            if size(temp.data.data_sorted,6) >20
                temp2 = temp.data.data_sorted(:,:,:,:,:,1:20);
                temp.data.data_sorted = temp2;
            end
            
            data_FF{i,ii} = temp.data.data_sorted;
            info_FF{i,ii} = temp.data.Info;
            info_FF{i,ii}.pos = file;
            info_FF{i,ii}.animal=i;
            info_FF{i,ii}.position=ii;
            
            % Looped cells
            if useLoop
                
                if ~isfield(temp.data,'beadID')
                    tempbeads=false(1,size(temp.data.data_sorted,1));
                    tempbeads(temp.data.beads_pos(:,3))=true;
                    temp.data.beadID.bead_pos=tempbeads';
                end
                
                info_FF{i,ii}.looped = temp.data.beadID.bead_pos.';
            end
            % XY coordinates
            
            for ROI=1:length(temp.data.s2p{1}.iscell)
                ROICoordyxtemp(ROI,:)   = temp.data.s2p{1}.stat{ROI}.med;
            end
            info_FF{i,ii}.ROICoordyx=ROICoordyxtemp(logical(temp.data.s2p{1}.iscell),:);
            
        end
    end
    clear temp
    
    data{1}=data_FB;
    data{2}=data_FF;
    info{1}=info_FB;
    info{2}=info_FF;
    clear data_FB; clear data_FF; clear info_FB; clear info_FF;
    
    
    n_animals_FB= size(data{1},1);n_animals_FF= size(data{2},1);
    for i=1:n_animals_FB
        n_pos_FB(i)=sum(~cellfun('isempty',data{1}(i,:)));
    end
    for i=1:n_animals_FF
        n_pos_FF(i)=sum(~cellfun('isempty',data{2}(i,:)));
    end
    n_pos{1}= n_pos_FB; n_pos{2}= n_pos_FF;
    n_ani{1}= n_animals_FB; n_ani{2}= n_animals_FF;
    
    
    % Make time matrices the same size
    for F=1:2
        for ani=1:n_ani{F}
            for pos=1:n_pos{F}(ani)
                if info{F}{ani,pos}.trial_length==4
                    data{F}{ani,pos}(:,end+1:75,:,:,:,:)=nan;
                elseif size(data{F}{ani,pos},2)>75
                    data{F}{ani,pos}=data{F}{ani,pos}(:,1:75,:,:,:,:);
                end
            end
        end
    end
    %
    
    disp('saving info ...')
    save([MainDir filesep 'info_' tv], 'info');

% Concatenate All sessions and animals
for F=1:2
    for ani=1:size(data{F},1)
        for pos=1:sum(~cellfun('isempty',data{F}(ani,:)))
            if size(data{F}{ani,pos},2)<75
                data{F}{ani,pos}(:,end+1:75,:,:,:,:)=NaN;
            end
        end
    end
end



for F=1:2
    for ani=1:size(info{F},1)
        MouseID{F}{:,ani}=info{F}{ani,1}.saveName(1,1:6); %MouseID{F}{:,ani}=info{F}{ani,1}.saveName(1,1:5);
    end
end