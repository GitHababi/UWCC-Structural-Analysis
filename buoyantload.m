function buoyantDist = buoyantload(draft, trim, CoordTable, xUnique, ymin, zmin, padding, structDensity)
%BUOYANTLOAD Calculates the buoyant distribution for a proposed trim

N = numel(xUnique);
buoyantDist = zeros(N,1);
ymin = abs(ymin);
yLoop = [ -ymin - padding,  ymin + padding, ymin + padding, -ymin - padding];

for i = 1:N
    zWaterline = draft + tan(trim) * xUnique(i);
    crossSectionShape = internalshape(crosssectioncoords(xUnique(i),CoordTable));
    waterAreaShape = polyshape(yLoop, [zmin - padding, zmin - padding, zWaterline, zWaterline]);
    
    intersection = intersect(crossSectionShape, waterAreaShape);
    if ~isempty(intersection.Vertices)
        buoyantDist(i) = area(intersection) * structDensity; 
    end
end
end

