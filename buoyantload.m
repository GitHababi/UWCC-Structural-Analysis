function buoyantDist = buoyantload(draft, trim, internalShapes, xUnique, waterbox,  padding, structDensity)
%BUOYANTLOAD Calculates the buoyant distribution for a proposed trim

N = numel(xUnique);
buoyantDist = zeros(N,1);

for i = 1:N
    waterline = draft + tan(trim) * xUnique(i);
    crossSectionShape = internalShapes(i);
    waterAreaShape = waterbox(waterline,padding);
    intersection = intersect(crossSectionShape, waterAreaShape);
    if ~isempty(intersection.Vertices)
        buoyantDist(i) = area(intersection) * structDensity; 
    end
end
end

