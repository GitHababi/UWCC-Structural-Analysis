close all;
clc; 

% Deprecated

%% Coord inputs
for i = 1:numel(xUnique)
    

end
  
disp("Cross Section Surface Area: " + crossSectionalArea(innerCoords,outerCoords) + " in^2");
[Centroid, MoI_canoe] = getMomentOfInertia(innerCoords,outerCoords);

disp("Moment of Inertia: " + MoI_canoe + " in^4");
disp("Centroid " + Centroid + " in from the keel");
    %% Plot

    allCoordsRight = [outerCoords; flipud(innerCoords)];
    allCoordsFull  = [allCoordsRight; -flipud(allCoordsRight(:,1)), allCoordsRight(:,2)];
    
    figure('Color','w');
    fill(allCoordsFull(:,1), allCoordsFull(:,2), [0 0 0], 'EdgeColor','k','FaceAlpha',0.5);
    hold on;
    plot(0, Centroid, 'ro', 'MarkerFaceColor','r','MarkerSize',8);
    yline(Centroid, 'r--', 'LineWidth',1.5);
    axis equal;
    xlabel('x (in)');
    ylabel('z (in)');
    title('Critical Section');
    legend('Cross Section','Centroid','Neutral Axis','Location','Best');
