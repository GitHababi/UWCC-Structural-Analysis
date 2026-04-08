function buoyantDist = buoyantload(draft, trim, internalShapes, xUnique, waterbox,  padding)
%BUOYANTLOAD Calculates the buoyant distribution for a proposed trim. all
%units in lbf, in, etc..

N = numel(xUnique);
buoyantDist = zeros(N,1);

for i = 1:N
    waterline = draft + tan(trim) * xUnique(i);
    crossSectionShape = internalShapes(i);
    waterAreaShape = waterbox(waterline,padding);
    intersection = intersect(crossSectionShape, waterAreaShape);
    if ~isempty(intersection.Vertices)
        buoyantDist(i) = area(intersection) * 62.4/1728; 
    end
end
end

