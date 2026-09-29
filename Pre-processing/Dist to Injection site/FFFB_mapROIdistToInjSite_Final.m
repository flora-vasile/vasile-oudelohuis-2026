
     
    DistInjSite2{F}=cell(n_ani{F},max(n_pos{F})); 
    
     for ani=1:n_ani{F}
         for pos=1:n_pos{F}(ani)
             try
                 DistInjSite2{F}{ani,pos}=DistInjSite{F}{ani,pos};
             catch
                 DistInjSite2{F}{ani,pos}=[];
             end
         end
     end
DistInjSite{F}=DistInjSite2{F};
    
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani) 
        ROICoordyx{F}{ani,pos}=info{F}{ani,pos}.ROICoordyx;

            clear points
            if ~isempty(DistInjSite{F}{ani,pos}) & ~isempty(info{F}{ani,pos}.ROICoordyx)

% Values of 2P field size
resolution=1024;
if resolution==1024
% for 1024*1024 imaging
   sizeTarget{F}{ani,pos}=1024;
   distanceline2P = 1024;
   HorAx2P=[0 512 ; 1024 512];
   center = [512 512];     % [x0, y0]
   square_size = 1024;        % Width/height of the square
elseif resolution==512
   % for 512*512 imaging
   sizeTarget{F}{ani,pos}=512;
   distanceline2P = 512;
   HorAx2P=[0 256 ; 512 256];
   center = [256 256];     % [x0, y0]
   square_size = 512;        % Width/height of the square
end  
% convert distance from global to local
% determine scaling factor
distanceHorzAxisIntrinsic{F}{ani,pos} = sqrt((DistInjSite{F}{ani,pos}.centerLine(2,1) - DistInjSite{F}{ani,pos}.centerLine(1,1))^2 + (DistInjSite{F}{ani,pos}.centerLine(2,2) - DistInjSite{F}{ani,pos}.centerLine(1,2))^2);

ScaleFactor{F}{ani,pos} = distanceline2P ./ distanceHorzAxisIntrinsic{F}{ani,pos}; 
distanceCentInjSite{F}{ani,pos}=DistInjSite{F}{ani,pos}.distanceToInjSite.*ScaleFactor{F}{ani,pos};
points(:,1)=ROICoordyx{F}{ani,pos}(:,2); points(:,2)=ROICoordyx{F}{ani,pos}(:,1);

% --- Known values ---

distance = distanceCentInjSite{F}{ani,pos};
angle_deg = DistInjSite{F}{ani,pos}.angleToInjSite;      % angle relative to centerLine (0-360)
%angle_deg =180;
% Custom horizontal vector (rightward)
refVec = HorAx2P(2,:) - HorAx2P(1,:);
%refVec = DistInjSite{F}{ani,pos}.centerLine(2,:) - DistInjSite{F}{ani,pos}.centerLine(1,:);
refVec = refVec / norm(refVec);  % unit vector

% Rotate by angle relative to the LEFT, which is refVec rotated by 180°
angle_rad = deg2rad(angle_deg);  % shift angle from "left" to standard rotation

% Rotate refVec counterclockwise by angle_rad
R = [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)];
directionVec = (R * refVec')';  % rotated direction from center to drawn point

% Get drawn point location
drawnPoint = center + distance * directionVec;



% Assume otherPoints is an Nx2 array
if resolution==1024
pix_to_um=0.5862;
elseif resolution==512 %(2x 4 planes)
pix_to_um=0.7813;
end

for i = 1:size(points,1)
differences{F}{ani,pos}(i,:) = points(i,:) - drawnPoint;
%distancesToInjSite{F}{ani,pos}(i) = 0.5862 * (sqrt(sum(differences{F}{ani,pos}(i).^2, 2)));  % Euclidean distances + % to convert from pixels to microns
distancesToInjSite{F}{ani,pos}(i) = pix_to_um * norm(differences{F}{ani,pos}(i,:));
end


if dispfig
% --- Visualization ---
figure; hold on; axis equal;

% Plot center
plot(center(1), center(2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(center(1)+5, center(2), 'Center', 'Color', 'green');

% Plot implied drawn point
plot(drawnPoint(1), drawnPoint(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(drawnPoint(1)+5, drawnPoint(2), 'Drawn Point', 'Color', 'red');

% Plot square
rectangle('Position', [center(1)-square_size/2, center(2)-square_size/2, square_size, square_size], ...
    'EdgeColor', 'k');

% Plot other points
plot(points(:,1), points(:,2), 'bs', 'MarkerSize', 20, 'LineWidth', 1.5, 'Marker', '.');

% Lines from other points to drawn point
for i = 1:size(points, 1)
    plot([points(i,1), drawnPoint(1)], [points(i,2), drawnPoint(2)], 'c-', 'LineWidth', 0.2);
    
end

% Plot the centerLine (custom horizontal axis)
plot(HorAx2P(:,1), HorAx2P(:,2), 'm-', 'LineWidth', 2);
midPoint = mean(HorAx2P, 1);
text(midPoint(1), midPoint(2)-10, 'centerLine', 'Color', 'magenta', 'FontSize', 10, 'FontWeight', 'bold');


% Line from center to drawn point
plot([center(1), drawnPoint(1)], [center(2), drawnPoint(2)], 'r--', 'LineWidth', 3);


title(['Distance between injection site and all ROIs ' MouseID{F}{ani} ' pos ' num2str(pos)]);
close
end
            end
        end
    end
    distancesToInjSite2{F}=cell(n_ani{F},max(n_pos{F})); 
    
     for ani=1:n_ani{F}
         for pos=1:n_pos{F}(ani)
             try
                 distancesToInjSite2{F}{ani,pos}=distancesToInjSite{F}{ani,pos};
             catch
                 distancesToInjSite2{F}{ani,pos}=[];
             end
         end
     end
distancesToInjSite{F}=distancesToInjSite2{F};    
    
% logical for selection of cells far from PM
ThDist=600;
if F==1
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if ~isempty(DistInjSite{F}{ani,pos}) & ~isempty(info{F}{ani,pos}.ROICoordyx)
                SelDistInjSite{F}{ani,pos}=false(size(distancesToInjSite{F}{ani,pos},2),1);
                for ROI=1:size(distancesToInjSite{F}{ani,pos},2)
                    if   distancesToInjSite{F}{ani,pos}(ROI)>ThDist
                        SelDistInjSite{F}{ani,pos}(ROI)=true;
                    end
                end
            else
             SelDistInjSite{F}{ani,pos}=true(length(ROICoordyx{F}{ani,pos}),1);   
            end
        end
    end
else
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            SelDistInjSite{F}{ani,pos}=true(length(ROICoordyx{F}{ani,pos}),1);
        end
    end
end
