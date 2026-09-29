%% if need to insert a new animal
% Define the position where you want to insert the empty row
insertRow = 3;

% Create an empty row (same number of columns as C)
emptyRow = cell(1, size(DistInjSite{1}, 2));

% Insert the empty row at the desired position
DistInjSite{1} = [DistInjSite{1}(1:insertRow-1, :); emptyRow; DistInjSite{1}(insertRow:end, :)];
clear emptyRow insertRow

%% if need to insert a new position of an existing animal
% Specify insertion point
row = 5;
col = 1;

% Shift elements in row 2 to the right from column 3 onward
DistInjSiteLearning{1}(row, col+1:end+1) = DistInjSiteLearning{1}(row, col:end);  % Shift right
DistInjSiteLearning{1}{row, col} = [];  % Insert empty cell


DistInjSiteLearning{1}{11, 2} =DistInjSiteLearning{1}{17, 2};
DistInjSiteLearning{1}{17, 2} =[];


%% Add rows and columns
% Column names
%rowNames = {'FV058', 'FV060', 'FV062', 'FV068', 'FV071', 'FV102', 'FV104', 'GF022', 'GF028', 'GF029', 'GF031', 'GF032'};

rowNames = {'CAM975','CAM976','FV009','FV038', 'FV040', 'FV106', 'FV107', 'FV110','FV114','FV126', 'FV127', 'FV130', 'FV131','FV140', 'FV141', 'GF030'};

colNames = {'Pos1', 'Pos2', 'Pos3', 'Pos4', 'Pos5', 'Pos6', 'Pos7'};

% Convert to table
TDistInjSiteLearning = cell2table(DistInjSiteLearning{1}, 'VariableNames', colNames, 'RowNames', rowNames);



%% Fill cell
for i=1:4
DistInjSiteLearning{1}{8,i}    =   params;
end
clear  params

%% Save
tv = datestr(now, 'yyyy_mm_dd');
save(['E:\2P\Pop Analysis\Learning\Sept 2024 Deconv\DistInjSiteLearning_' tv],'DistInjSiteLearning')

