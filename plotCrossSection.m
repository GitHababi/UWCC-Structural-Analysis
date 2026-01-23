function plotCrossSection(innerCoords,outerCoords)
%PLOTCROSSSECTION Summary of this function goes here plots 
%   Detailed explanation goes here
    hold on 
    plot(innerCoords(:,1),innerCoords(:,2),"b*");
    plot(flipud(outerCoords(:,1)),flipud(outerCoords(:,2)),"r*");
    hold off
end

