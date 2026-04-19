function [V,M] = shearmoment(xPositions,distributedLoads,pointLoads)
% xPositions -> vector of all xPositions, center around zero.
% distributedLoads -> vector of all distributedLoads corresponding to
% xPositions
% pointLoads -> 2xN matrix where first row is position, second row is
% force. 
    V = zeros(numel(xPositions),1);
    M = zeros(numel(xPositions),1);

    V = cumtrapz(xPositions,distributedLoads);

    for j = 1:size(pointLoads,2)
            for k = 1:numel(xPositions)
                if xPositions(k) > pointLoads(1,j)
                    V(k) = V(k) + pointLoads(2,j);
                end
            end
    end

    M = cumtrapz(xPositions,V);
end

