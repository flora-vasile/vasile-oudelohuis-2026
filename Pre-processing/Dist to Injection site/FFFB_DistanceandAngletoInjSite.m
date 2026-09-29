% Load image and params
[imageFile, imagePath] = uigetfile('*.png', 'Select Registered Image');
if isequal(imageFile,0), return; end
img = imread(fullfile(imagePath, imageFile));

try
load(fullfile(imagePath,'transform_params.mat'), 'params'); 
catch
disp('no params file in directory')
[paramFile, paramPath] = uigetfile('transform_params.mat', 'Select Params File');
if isequal(paramFile,0), return; end
load(fullfile(paramPath, paramFile), 'params');
end

% Show image
fig = figure;
imshow(img); hold on;

% Plot center line and center
plot(params.centerLine(:,1), params.centerLine(:,2), 'r-', 'LineWidth', 2);
plot(params.centerOnImage2(1), params.centerOnImage2(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);

% --- Let user select a point
disp('Click a point on the image...');
[x, y] = ginput(1);  % User selects one point

if isempty([x y])
    disp('No point selected.');
    return;
end

drawnPoint = [x, y];
center = params.centerOnImage2;

% Compute distance from center to point
distance = norm(drawnPoint - center);

% Define the custom horizontal vector from centerLine
refVec = params.centerLine(2,:) - params.centerLine(1,:);
refVec = refVec / norm(refVec);  % Normalize

% Vector from center to selected point
vecToPoint = drawnPoint - center;
vecToPoint = vecToPoint / norm(vecToPoint);  % Normalize

% Compute signed angle from refVec to vecToPoint
dotProd = dot(refVec, vecToPoint);
det = -(refVec(1)*vecToPoint(2) - refVec(2)*vecToPoint(1));  % Inverted 2D cross

angle_rad = atan2(det, dotProd);  % atan2 of det and dot gives signed angle
angle_deg = rad2deg(angle_rad);

% Convert to [0, 360] range
if angle_deg < 0
    angle_deg = angle_deg + 360;
end


% Now enforce your 0–180 (top) and 180–360 (bottom) rule
% Assuming "top half" means above the centerLine in image coordinates:
% We'll use the sign of the cross product to split top/bottom.
% The rule is already correct as written.

% Store results
params.drawnPoint = drawnPoint;
params.distanceToInjSite = distance;
params.angleToInjSite = angle_deg;

% Plot result
plot(drawnPoint(1), drawnPoint(2), 'ms', 'MarkerSize', 8, 'LineWidth', 2);
plot([center(1), drawnPoint(1)], [center(2), drawnPoint(2)], 'c--', 'LineWidth', 1.5);

% Save overlay image
saveas(fig, fullfile(imagePath, 'analysis_overlay_point.png'));

% Save updated params
save(fullfile(imagePath, 'transform_params_with_analysis_point.mat'), 'params');

disp(['Saved. Distance: ' num2str(distance, '%.2f') ...
      ', Angle: ' num2str(angle_deg, '%.2f') '°']);

clear all
