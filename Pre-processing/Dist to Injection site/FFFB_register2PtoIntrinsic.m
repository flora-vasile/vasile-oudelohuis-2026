function overlay_gui_transform_save()

    % --- Load fixed (color) image (Image 2)
    [fixedFile, fixedPath] = uigetfile({'*.tif;*.png;*.jpg'}, 'Select Reference Image (Image 2)');
    if isequal(fixedFile,0), return; end
    fixed = im2double(imread(fullfile(fixedPath, fixedFile)));
    
    % --- Load moving (grayscale) image (Image 1)
    try
    movingPath=fixedPath;
    movingFile='target_p2.tif';
    moving = im2double(imread(fullfile(movingPath, movingFile)));
    catch
    [movingFile, movingPath] = uigetfile({'*.tif;*.png;*.jpg'}, 'Select Image to Overlay (Image 1)');
    if isequal(movingFile,0), return; end
    
    moving = im2double(imread(fullfile(movingPath, movingFile)));
    end
    if size(moving,3) == 3
        moving = rgb2gray(moving);
    end

    % Normalize moving image
    moving = (moving - min(moving(:))) / (max(moving(:)) - min(moving(:)));

    % Initial parameters
    dx = 0;
    dy = 0;
    angle = 0;
    scale = 0.5;
    dragging = false;
    prevPoint = [];
  
    % Create GUI
    fig = figure('Name','Manual Registration','NumberTitle','off','Position',[100 100 1000 800]);
    
    
    % Shared variables
    dx = 0; dy = 0;
    prevPt = [];
    
    
    ax = axes('Parent', fig, 'Position', [0.05 0.2 0.9 0.75]);
    imshow(fixed, 'Parent', ax); hold on;

    % Enable dragging
    set(fig, 'WindowButtonDownFcn', @mouseDown);

    % NESTED FUNCTION: mouseDown
    function mouseDown(~, ~)
        prevPt = get(gca, 'CurrentPoint');
        set(fig, 'WindowButtonMotionFcn', @mouseDrag);
        set(fig, 'WindowButtonUpFcn', @(~,~) set(fig, 'WindowButtonMotionFcn', ''));
    end

    % NESTED FUNCTION: mouseDrag
    function mouseDrag(~, ~)
        cp = get(gca, 'CurrentPoint');
        dx = dx + (cp(1,1) - prevPt(1,1));
        dy = dy + (cp(1,2) - prevPt(1,2));
        prevPt = cp;
        updateDisplay();  % assumes this updates the image position
    end
    
    
    % Sliders
    uicontrol('Style','text','String','Rotation (deg)','Position',[20 60 100 20]);
    rotSlider = uicontrol('Style','slider','Min',-180,'Max',180,'Value',0,...
        'Position',[130 60 300 20],'Callback',@sliderCallback);
 
    rotText = uicontrol('Style','text','String','0°','Position',[440 60 50 20]);   

    uicontrol('Style','text','String','Scale','Position',[20 30 100 20]);
    scaleSlider = uicontrol('Style','slider','Min',0.01,'Max',1,'Value',0.5,...
        'Position',[130 30 300 20],'Callback',@sliderCallback);

    scaleText = uicontrol('Style','text','String','0.50','Position',[440 30 50 20]);
    
    uicontrol('Style', 'pushbutton', ...
          'String', 'Save Registered Image', ...
          'Units', 'normalized', ...
          'Position', [0.8 0.05 0.15 0.05], ...
          'Callback', @saveFinalImage);

    uicontrol('Style', 'togglebutton', 'String', 'Zoom Mode', ...
    'Units', 'normalized', 'Position', [0.8 0.9 0.15 0.05], ...
    'Callback', @toggleZoom);

    % Display overlay
    hOverlay = imshow(moving, 'Parent', ax);
    set(hOverlay, 'AlphaData', 0.5, 'ButtonDownFcn', @startDrag);

    % Mouse interaction
    set(fig, 'WindowButtonMotionFcn', @draggingFcn);
    set(fig, 'WindowButtonUpFcn', @stopDrag);

    % --- Initial Display
    updateDisplay();

    % --- CALLBACK FUNCTIONS ---
    function sliderCallback(~,~)
        angle = get(rotSlider, 'Value');
        scale = get(scaleSlider, 'Value');
        set(scaleText, 'String', sprintf('%.2f', scale));  % Update live display
        set(rotText, 'String', sprintf('%.0f°', angle));
        updateDisplay();
    end



function toggleZoom(src, ~)
    if get(src, 'Value')
        zoom on;
        set(src, 'BackgroundColor', [0.8 1 0.8]);
        set(fig, 'WindowButtonDownFcn', '');  % Disable drag while zooming
    else
        zoom off;
        set(src, 'BackgroundColor', [1 1 1]);
        set(fig, 'WindowButtonDownFcn', @mouseDown);  % Re-enable dragging
    end
end

  function updateDisplay()
    % Step 1: Scale
    scaled = imresize(moving, scale);
    sz = size(scaled);

    % Step 2: Pad to allow rotation without clipping
    diagLen = ceil(norm(sz));  % Diagonal size of the image
    padX = ceil((diagLen - sz(2)) / 2);
    padY = ceil((diagLen - sz(1)) / 2);
    padded = padarray(scaled, [padY, padX], 0, 'both');

    % Step 3: Rotation centered on padded image
    szPadded = size(padded);
    center = [szPadded(2)/2, szPadded(1)/2];
    theta = deg2rad(angle);

    % Centered rotation matrix
    T1 = [1 0 0; 0 1 0; -center 1];        % Translate to origin
    R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];  % Rotate
    T2 = [1 0 0; 0 1 0; center 1];         % Translate back
    TformMatrix = T1 * R * T2;
    tform = affine2d(TformMatrix);

    % Step 4: Rotate with padded size
    rotated = imwarp(padded, tform, 'OutputView', imref2d(szPadded));

    % Step 5: Update overlay
    set(hOverlay, 'CData', rotated);
    set(hOverlay, 'XData', [1 size(rotated,2)] + dx);
    set(hOverlay, 'YData', [1 size(rotated,1)] + dy);
    
% For grayscale:
if size(rotated, 3) == 1
    alphaMask = rotated > 0;
else
    % For RGB: use brightness as mask
    alphaMask = max(rotated, [], 3) > 0;
end

set(hOverlay, 'AlphaData', alphaMask);
  end


function saveFinalImage(~, ~)
   % Step 1: Scale
scaled = imresize(moving, scale);

% Step 2: Pad the scaled image
sz = size(scaled);
diagLen = ceil(norm(sz));
padX = ceil((diagLen - sz(2))/2);
padY = ceil((diagLen - sz(1))/2);
padded = padarray(scaled, [padY, padX], 0, 'both');

% Step 3: Rotation
szPadded = size(padded);
center = [szPadded(2)/2, szPadded(1)/2];  % [x, y]
theta = deg2rad(angle);
T1 = [1 0 0; 0 1 0; -center 1];
R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
T2 = [1 0 0; 0 1 0; center 1];
T = T1 * R * T2;
tform = affine2d(T);

% Step 4: Warp
rotated = imwarp(padded, tform, 'OutputView', imref2d(szPadded));

% Step 5: Translate
canvas = zeros(size(fixed));
xRange = round(1 + dx : size(rotated,2) + dx);
yRange = round(1 + dy : size(rotated,1) + dy);
xRange = max(1, min(xRange, size(canvas,2)));
yRange = max(1, min(yRange, size(canvas,1)));
canvas(yRange, xRange) = rotated(1:length(yRange), 1:length(xRange));

% Step 6: Blend
fusedImage = imfuse(fixed, canvas, 'blend');

% Step 7: Save image
[file, path] = uiputfile('registeredImage.png', 'Save registered image');

if ischar(file)
    imwrite(fusedImage, fullfile(path, file));
end

% Step 8: Compute center of image 1 after transform
originalCenter = size(moving) / 2;
scaledCenter = originalCenter * scale;
paddedCenter = scaledCenter + [padX, padY];
rotMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
rotatedCenter = (rotMatrix * paddedCenter')';
finalCenter = rotatedCenter + [dx, dy];  % <-- this is the correct center

    % Center line points (horizontal line)
    x1 = 1;
    x2 = size(moving, 2);
    y = originalCenter(1);  % center row

    pt1 = [x1, y];
    pt2 = [x2, y];

    % Scale
    pt1 = pt1 * scale;
    pt2 = pt2 * scale;

    % Pad
    pt1 = pt1 + [padX, padY];
    pt2 = pt2 + [padX, padY];

    % Rotate
    pt1_rot = transformPointsForward(tform, pt1);
    pt2_rot = transformPointsForward(tform, pt2);

    % Translate
    pt1_final = pt1_rot + [dx, dy];
    pt2_final = pt2_rot + [dx, dy];

%     % Optionally draw the line on the fused image
%     figure; imshow(fusedImage); hold on;
%     plot([pt1_final(1), pt2_final(1)], [pt1_final(2), pt2_final(2)], 'r-', 'LineWidth', 2);
%     plot( finalCenter(1), finalCenter(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);  % Green circle
%     title('Registered Image with Center and Center Line');

    % Save all to file
    params = struct();
    params.angle = angle;
    params.scale = scale;
    params.dx = dx;
    params.dy = dy;
    params.centerLine = [pt1_final; pt2_final];  % [start; end]
    params.centerOnImage2 = mean(params.centerLine, 1); %params.centerOnImage2 = finalCenter;

    save(fullfile(path, 'transform_params.mat'), 'params');
    
    % Optionally draw the line on the fused image
    figure; imshow(fusedImage); hold on;
    plot(params.centerOnImage2(1), params.centerOnImage2(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);  % center dot
plot(params.centerLine(:,1), params.centerLine(:,2), 'r-', 'LineWidth', 2);  % center line
    title('Registered Image with Center and Center Line');
    
end



    function startDrag(~,~)
        dragging = true;
        prevPoint = get(ax, 'CurrentPoint');
    end

    function draggingFcn(~,~)
        if dragging
            currPoint = get(ax, 'CurrentPoint');
            delta = currPoint(1,1:2) - prevPoint(1,1:2);
            dx = dx + delta(1);
            dy = dy + delta(2);
            prevPoint = currPoint;
            updateDisplay();
        end
    end

    function stopDrag(~,~)
        dragging = false;
    end




end
