function [innerCoords,outerCoords] = crosssectioncoords(x,CoordTable)
%CROSSSECTIONCOORDS Pulls cross section coordinates from a given x station.
%   Detailed explanation goes here
    xCoordinates = CoordTable.LongitudinalPosition;
    yCoordinates = CoordTable.Offset;
    zCoordinates = CoordTable.Height;
    side         = CoordTable.Side;

    outerIndex = (xCoordinates == x) & (side == "Outer");
    innerIndex = (xCoordinates == x) & (side == "Inner");

    [yInner, sortIndex] = sort(yCoordinates(innerIndex), 'descend');
    zInner = zCoordinates(innerIndex); zInner = zInner(sortIndex);
    [yOuter, sortIndex] = sort(yCoordinates(outerIndex), 'descend');
    zOuter = zCoordinates(outerIndex); zOuter = zOuter(sortIndex);
    
    innerCoords = [yInner zInner];
    outerCoords = [yOuter zOuter];
end

