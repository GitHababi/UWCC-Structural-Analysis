close all;
clc;


coordinates = readtable("HullCoordinates.csv");
x_coordinates = coordinates(:,"LongitudalPosition");
unique_x = table2array(unique(x_coordinates));
moi = zeros(size(unique_x));
centroids = zeros(size(unique_x));
for i = 1:numel(unique_x)
    current_x = unique_x(i);
    inner_rows = (coordinates.LongitudalPosition == current_x & matches(coordinates.Side,"Inner"));
    inner_coord_table = coordinates(inner_rows,:);
    innerCoords = inner_coord_table{:,["Offset","Height"]};
    
    outer_rows = (coordinates.LongitudalPosition == current_x & matches(coordinates.Side,"Outer"));
    outer_coord_table = coordinates(outer_rows,:);
    outerCoords = outer_coord_table{:,["Offset","Height"]};
    if current_x == -54
        plotCrossSection(innerCoords,outerCoords);
        disp(size(innerCoords));
        disp(size(outerCoords));
    end
    [centroids(i), moi(i), ~] = generateCrossSectionProperties(innerCoords,outerCoords);
end


figure;
hold on
plot(unique_x,moi);
plot(unique_x,centroids .* 100,"Color","r");
hold off;