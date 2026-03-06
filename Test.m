%% init

close all;
clc;


coordinates = readtable("HullCoordinates.csv");
x_coordinates = coordinates(:,"LongitudalPosition");
unique_x = table2array(unique(x_coordinates));

%% code to find total surface area on outside of canoe

% find arc length for each unique station
arc_lengths = zeros(size(unique_x));
for i = 1:numel(unique_x)
    current_x = unique_x(i);

    outer_rows = (coordinates.LongitudalPosition == current_x & matches(coordinates.Side,"Outer"));
    outer_coord_table = coordinates(outer_rows,:);
    outerCoords = outer_coord_table{:,["Offset","Height"]};
    arc_length = 0;
    [num_points, ~] = size(outerCoords);
    for j = 1:num_points-1
        arc_length = arc_length + pdist(outerCoords(j:j+1,:),'euclidean');
    end
    arc_lengths(i) = arc_length;
end

figure;
hold on
plot(unique_x,arc_lengths);
hold off;

surface_areas = cumtrapz(unique_x,arc_lengths);
disp("Surface area on outside of canoe: " + max(surface_areas) + " (in^2), " + max(surface_areas)/144 + " (ft^2)");
%% code to find all  mois and centroids
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