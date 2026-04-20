function moi = moinertia(innerCoords,outerCoords)
%COMPUTEMOI Computes the moment of inertia (i.e. second moment of area) of a cross section about it's centroid given inner and outer
%coordinates.
%   Detailed explanation goes here
    
    outside = flipud(outerCoords);
    inside  = flipud(innerCoords);
    nOuter = size(outside,1);
    nInner = size(inside,1);
    
    [~, centroid] = csa(innerCoords,outerCoords);

    outerI = 0;

    for i = 1:nOuter-1
        y1 = abs(outside(i, 1));
        y2 = abs(outside(i+1, 1));
        z1 = outside(i, 2);
        z2 = outside(i+1, 2);
        h  = abs(z2 - z1);
        elementCentroid = (z1 + z2)/2;
        area = (y1 + y2)/2 * h;
        outerI = outerI + (y1+y2)/2*(h)^3/12 + area * (elementCentroid - centroid)^2;
    end
    
    innerI = 0;

    for i = 1:nInner-1
        y1 = abs(inside(i, 1));
        y2 = abs(inside(i+1, 1));
        z1 = inside(i, 2);
        z2 = inside(i+1, 2);
        h  = abs(z2 - z1);
        elementCentroid = (z1 + z2)/2;
        area = (y1 + y2)/2 * h;
        innerI = innerI + (y1+y2)/2*(h)^3/12 + area * (elementCentroid - centroid)^2;
    end
    
    moi = 2 * abs(innerI - outerI);
end

