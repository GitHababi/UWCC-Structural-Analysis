function moment = momentabout(X,load,xPositions)
%MOMENTABOUT Calculates the net moment about a point X from the distributed
% load with given xPositions
%   Detailed explanation goes here
    load = load .* (xPositions - X);
    moment = trapz(xPositions, load);
end
