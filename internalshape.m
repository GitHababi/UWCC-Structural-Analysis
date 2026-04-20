function volume = internalshape(outerCoords)
% INTERNALVOL Calculates the total volume inside cross section
%   Detailed explanation goes here
    y_all = [outerCoords(:,1); -flip(outerCoords(:,1))];
    z_all = [outerCoords(:,2); flip(outerCoords(:,2))];
    volume = polyshape(y_all,z_all, 'Simplify', true);
end

