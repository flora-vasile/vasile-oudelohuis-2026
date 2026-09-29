dispfig=0;
if F==2
    if PMcoord
        for ani=1:n_ani{F}
            for pos=1:n_pos{F}(ani)
                
                ROICoordyx{F}{ani,pos}=info{F}{ani,pos}.ROICoordyx;
                % Define the perimeter of the area
                
                    
                xPerimeter = xyPM{F}{ani,pos}(:,1); 
                yPerimeter = xyPM{F}{ani,pos}(:,2);
                
               k = convhull(xPerimeter,yPerimeter);
               xBoundary = xPerimeter(k);
               yBoundary = yPerimeter(k);
                % Plot the points and convex hull

                
                % Define the points to check
                xPoints = ROICoordyx{F}{ani,pos}(:,2); 
                yPoints = ROICoordyx{F}{ani,pos}(:,1); 
                
                % Use inpolygon to determine points inside the area
                [in, on] = inpolygon(xPoints, yPoints, xBoundary, yBoundary);
                
                % Points inside, on the perimeter, and outside
                insidePointsX = xPoints(in);
                insidePointsY = yPoints(in);
                
                
                InPM{F}{ani,pos}=in|on;
                
                
                onPerimeterPointsX = xPoints(on);
                onPerimeterPointsY = yPoints(on);
                
                outsidePointsX = xPoints(~(in | on));
                outsidePointsY = yPoints(~(in | on));
                
                if dispfig
                % Plotting
                figure;
                hold on;
                
                % Plot the polygon edges
                plot(xBoundary,  yBoundary, 'b-', 'LineWidth', 2, 'DisplayName', 'Polygon Perimeter');
                
                % Plot points
                plot(insidePointsX, insidePointsY, 'go', 'MarkerSize', 8, 'DisplayName', 'Inside Points');
                plot(onPerimeterPointsX, onPerimeterPointsY, 'ro', 'MarkerSize', 8, 'DisplayName', 'On Perimeter Points');
                plot(outsidePointsX, outsidePointsY, 'kx', 'MarkerSize', 8, 'DisplayName', 'Outside Points');
                
                % Add labels and legend
                xlabel('X');
                ylabel('Y');
                title(['Polygon and Points '  MouseID{F}{ani}  ' Pos ' num2str(pos)] );
                legend('show');
                grid on;
                axis equal;
                hold off;
                set(gca, 'YDir','reverse')
                end
            end
        end
    end
else
    
        for ani=1:n_ani{F}
            for pos=1:n_pos{F}(ani)
               InPM{F}{ani,pos}=true(length(info{F}{ani,pos}.ROICoordyx),1); 
            end
        end
end